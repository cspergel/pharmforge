"""
PubMed Adapter - Fetches scientific literature from NCBI PubMed
Uses E-utilities API for literature search, MeSH terms, and abstracts
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
import xml.etree.ElementTree as ET
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class PubMedAdapter(AdapterProtocol):
    """
    Adapter for NCBI PubMed E-utilities API
    Searches for scientific literature related to drugs, diseases, genes, and keywords

    Capabilities:
    - Literature search with boolean logic
    - MeSH term extraction
    - Abstract retrieval
    - Citation information
    - Publication date filtering
    """

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        """
        Initialize PubMed adapter

        Args:
            email: Email for NCBI (recommended by NCBI)
            api_key: NCBI API key for higher rate limits (10 req/s vs 3 req/s)
        """
        super().__init__(
            name="pubmed",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.34 if not api_key else 0.11,  # 3 req/s or 10 req/s with key
                "timeout": 30,
                "max_results": 100,
                "email": email or "pharmforge@example.com",
                "api_key": api_key
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid search query

        Args:
            input_data: Search query string or dict with query parameters

        Returns:
            True if valid, False otherwise
        """
        if isinstance(input_data, str):
            return len(input_data.strip()) > 0
        elif isinstance(input_data, dict):
            return "query" in input_data and len(input_data["query"].strip()) > 0
        return False

    async def _search_pubmed_async(
        self,
        query: str,
        max_results: int = 100,
        sort_by: str = "relevance",
        filters: Optional[Dict[str, Any]] = None
    ) -> Optional[List[str]]:
        """
        Search PubMed and return PMIDs

        Args:
            query: Search query (supports boolean operators: AND, OR, NOT)
            max_results: Maximum number of results to return
            sort_by: Sort order ("relevance", "pub_date", "first_author")
            filters: Optional filters (date_from, date_to, article_type, etc.)

        Returns:
            List of PMIDs or None on error
        """
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "email": self.config.get("email")
        }

        # Add API key if available
        if self.config.get("api_key"):
            params["api_key"] = self.config["api_key"]

        # Add sort parameter
        sort_map = {
            "relevance": "relevance",
            "pub_date": "pub+date",
            "first_author": "first+author"
        }
        params["sort"] = sort_map.get(sort_by, "relevance")

        # Add date filters if provided
        if filters:
            if "date_from" in filters:
                params["mindate"] = filters["date_from"]
            if "date_to" in filters:
                params["maxdate"] = filters["date_to"]
            if "date_from" in filters or "date_to" in filters:
                params["datetype"] = "pdat"  # Publication date

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.ESEARCH_URL, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.error(f"PubMed search failed with status {response.status}")
                        return None

                    data = await response.json()
                    result = data.get("esearchresult", {})
                    pmids = result.get("idlist", [])
                    count = int(result.get("count", 0))

                    logger.info(f"PubMed: Found {count} articles for query: {query[:50]}...")
                    logger.info(f"PubMed: Returning {len(pmids)} PMIDs")

                    return pmids

        except asyncio.TimeoutError:
            logger.error(f"PubMed: Timeout searching for: {query}")
            return None
        except Exception as e:
            logger.error(f"PubMed: Error searching for {query}: {e}")
            return None

    async def _fetch_article_details_async(self, pmids: List[str]) -> Optional[List[Dict[str, Any]]]:
        """
        Fetch detailed article information for given PMIDs

        Args:
            pmids: List of PubMed IDs

        Returns:
            List of article details or None on error
        """
        if not pmids:
            return []

        # Fetch in batches of 200 (NCBI limit)
        batch_size = 200
        all_articles = []

        for i in range(0, len(pmids), batch_size):
            batch_pmids = pmids[i:i + batch_size]
            pmid_str = ",".join(batch_pmids)

            params = {
                "db": "pubmed",
                "id": pmid_str,
                "retmode": "xml",
                "email": self.config.get("email")
            }

            if self.config.get("api_key"):
                params["api_key"] = self.config["api_key"]

            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

            try:
                async with aiohttp.ClientSession() as session:
                    async with session.get(self.EFETCH_URL, params=params, timeout=timeout) as response:
                        if response.status != 200:
                            logger.error(f"PubMed fetch failed with status {response.status}")
                            continue

                        xml_content = await response.text()
                        articles = self._parse_pubmed_xml(xml_content)
                        all_articles.extend(articles)

                        # Rate limiting between batches
                        if i + batch_size < len(pmids):
                            await asyncio.sleep(self.config.get("rate_limit_delay", 0.34))

            except asyncio.TimeoutError:
                logger.error(f"PubMed: Timeout fetching articles for batch starting at {i}")
                continue
            except Exception as e:
                logger.error(f"PubMed: Error fetching articles: {e}")
                continue

        return all_articles if all_articles else None

    def _parse_pubmed_xml(self, xml_content: str) -> List[Dict[str, Any]]:
        """
        Parse PubMed XML response and extract article information

        Args:
            xml_content: XML string from PubMed API

        Returns:
            List of parsed article dictionaries
        """
        articles = []

        try:
            root = ET.fromstring(xml_content)

            for article_elem in root.findall(".//PubmedArticle"):
                try:
                    article_data = {}

                    # Get PMID
                    pmid_elem = article_elem.find(".//PMID")
                    article_data["pmid"] = pmid_elem.text if pmid_elem is not None else None

                    # Get Article details
                    article = article_elem.find(".//Article")
                    if article is not None:
                        # Title
                        title_elem = article.find(".//ArticleTitle")
                        article_data["title"] = title_elem.text if title_elem is not None else None

                        # Abstract
                        abstract_parts = []
                        for abstract_text in article.findall(".//AbstractText"):
                            label = abstract_text.get("Label", "")
                            text = abstract_text.text or ""
                            if label:
                                abstract_parts.append(f"{label}: {text}")
                            else:
                                abstract_parts.append(text)
                        article_data["abstract"] = " ".join(abstract_parts) if abstract_parts else None

                        # Authors
                        authors = []
                        for author in article.findall(".//Author"):
                            lastname = author.find("LastName")
                            forename = author.find("ForeName")
                            if lastname is not None:
                                name = lastname.text
                                if forename is not None:
                                    name = f"{forename.text} {name}"
                                authors.append(name)
                        article_data["authors"] = authors

                        # Journal
                        journal = article.find(".//Journal")
                        if journal is not None:
                            title = journal.find(".//Title")
                            article_data["journal"] = title.text if title is not None else None

                        # Publication date
                        pub_date = article.find(".//PubDate")
                        if pub_date is not None:
                            year = pub_date.find("Year")
                            month = pub_date.find("Month")
                            day = pub_date.find("Day")

                            date_parts = []
                            if year is not None:
                                date_parts.append(year.text)
                            if month is not None:
                                date_parts.append(month.text)
                            if day is not None:
                                date_parts.append(day.text)

                            article_data["publication_date"] = "-".join(date_parts) if date_parts else None

                    # Get MeSH terms
                    mesh_terms = []
                    for mesh in article_elem.findall(".//MeshHeading"):
                        descriptor = mesh.find("DescriptorName")
                        if descriptor is not None:
                            mesh_terms.append({
                                "term": descriptor.text,
                                "major_topic": descriptor.get("MajorTopicYN", "N") == "Y"
                            })
                    article_data["mesh_terms"] = mesh_terms

                    # Get DOI
                    for article_id in article_elem.findall(".//ArticleId"):
                        if article_id.get("IdType") == "doi":
                            article_data["doi"] = article_id.text
                            break

                    articles.append(article_data)

                except Exception as e:
                    logger.warning(f"Error parsing individual article: {e}")
                    continue

            logger.info(f"PubMed: Successfully parsed {len(articles)} articles from XML")

        except Exception as e:
            logger.error(f"PubMed: Error parsing XML: {e}")

        return articles

    def _summarize_results(self, articles: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Create summary statistics from article results

        Args:
            articles: List of article dictionaries

        Returns:
            Summary dictionary with statistics
        """
        if not articles:
            return {
                "total_articles": 0,
                "articles_with_abstracts": 0,
                "common_mesh_terms": [],
                "journals": [],
                "date_range": None
            }

        # Count articles with abstracts
        with_abstracts = sum(1 for a in articles if a.get("abstract"))

        # Extract common MeSH terms
        mesh_counts = {}
        for article in articles:
            for mesh in article.get("mesh_terms", []):
                term = mesh.get("term")
                if term:
                    mesh_counts[term] = mesh_counts.get(term, 0) + 1

        # Top 10 MeSH terms
        common_mesh = sorted(mesh_counts.items(), key=lambda x: x[1], reverse=True)[:10]

        # Extract journals
        journal_counts = {}
        for article in articles:
            journal = article.get("journal")
            if journal:
                journal_counts[journal] = journal_counts.get(journal, 0) + 1

        top_journals = sorted(journal_counts.items(), key=lambda x: x[1], reverse=True)[:5]

        # Date range
        dates = [a.get("publication_date") for a in articles if a.get("publication_date")]
        date_range = None
        if dates:
            dates_sorted = sorted(dates)
            date_range = {"earliest": dates_sorted[0], "latest": dates_sorted[-1]}

        return {
            "total_articles": len(articles),
            "articles_with_abstracts": with_abstracts,
            "common_mesh_terms": [{"term": term, "count": count} for term, count in common_mesh],
            "top_journals": [{"journal": j, "count": c} for j, c in top_journals],
            "date_range": date_range
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute PubMed literature search

        Args:
            input_data: Search query string or dict with parameters
                       String: "aspirin AND cancer"
                       Dict: {"query": "aspirin", "max_results": 50, "date_from": "2020/01/01"}
            **kwargs: Additional parameters
                     - max_results: Maximum number of results (default: 100)
                     - sort_by: Sort order ("relevance", "pub_date", "first_author")
                     - filters: Dict with date_from, date_to, etc.
                     - fetch_details: Whether to fetch full article details (default: True)

        Returns:
            AdapterResult containing article information
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid search query"
            )

        # Parse input
        if isinstance(input_data, str):
            query = input_data
            max_results = kwargs.get("max_results", self.config.get("max_results", 100))
            sort_by = kwargs.get("sort_by", "relevance")
            filters = kwargs.get("filters")
        else:
            query = input_data["query"]
            max_results = input_data.get("max_results", kwargs.get("max_results", 100))
            sort_by = input_data.get("sort_by", kwargs.get("sort_by", "relevance"))
            filters = input_data.get("filters", kwargs.get("filters"))

        fetch_details = kwargs.get("fetch_details", True)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.34)
        await asyncio.sleep(rate_delay)

        # Search PubMed
        pmids = await self._search_pubmed_async(query, max_results, sort_by, filters)

        if pmids is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search PubMed",
                metadata={"source": "pubmed", "query": query}
            )

        if not pmids:
            return AdapterResult(
                success=True,
                data={
                    "pmids": [],
                    "articles": [],
                    "summary": self._summarize_results([]),
                    "query": query
                },
                metadata={"source": "pubmed", "query": query, "total_results": 0}
            )

        # Fetch article details if requested
        articles = []
        if fetch_details:
            articles = await self._fetch_article_details_async(pmids) or []

        # Create summary
        summary = self._summarize_results(articles)

        result_data = {
            "pmids": pmids,
            "articles": articles,
            "summary": summary,
            "query": query
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "pubmed",
                "query": query,
                "adapter_version": self.version,
                "total_results": len(pmids),
                "articles_fetched": len(articles)
            }
        )
