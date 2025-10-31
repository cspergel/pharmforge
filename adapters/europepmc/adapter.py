"""
Europe PMC Adapter - Fetches full-text articles with chemical entity mentions
Uses Europe PMC RESTful Web Service API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class EuropePMCAdapter(AdapterProtocol):
    """
    Adapter for Europe PMC RESTful API

    Capabilities:
    - Full-text article search
    - Chemical entity extraction and mentions
    - MeSH term support
    - Open access article identification
    - Citation retrieval
    - Cross-references to other databases (PubMed, DOI, etc.)
    - Advanced query syntax with boolean operators
    """

    SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    FULLTEXT_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{source}/{id}/fullTextXML"
    ANNOTATIONS_URL = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds"

    def __init__(self, email: Optional[str] = None):
        """
        Initialize Europe PMC adapter

        Args:
            email: Contact email (optional but recommended)
        """
        super().__init__(
            name="europepmc",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.2,  # 5 requests/second
                "timeout": 30,
                "max_results": 100,
                "email": email or "pharmforge@example.com"
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

    def _build_query_string(
        self,
        query: str,
        filters: Optional[Dict[str, Any]] = None
    ) -> str:
        """
        Build Europe PMC query string with filters

        Args:
            query: Base search query (supports boolean operators: AND, OR, NOT)
            filters: Optional filters

        Returns:
            Complete query string
        """
        query_parts = [query]

        if filters:
            # Open access filter
            if filters.get("open_access"):
                query_parts.append("OPEN_ACCESS:Y")

            # Publication date range
            if "date_from" in filters or "date_to" in filters:
                date_from = filters.get("date_from", "*")
                date_to = filters.get("date_to", "*")
                query_parts.append(f"FIRST_PDATE:[{date_from} TO {date_to}]")

            # Source type (e.g., "MED" for PubMed, "PMC" for PMC)
            if "source" in filters:
                query_parts.append(f"SRC:{filters['source']}")

            # Publication type
            if "pub_type" in filters:
                query_parts.append(f"PUB_TYPE:\"{filters['pub_type']}\"")

            # Has chemical annotations
            if filters.get("has_chemicals"):
                query_parts.append("HAS_CHEM:Y")

            # Has full text
            if filters.get("has_fulltext"):
                query_parts.append("HAS_FT:Y")

        return " AND ".join(query_parts)

    async def _search_europepmc_async(
        self,
        query: str,
        max_results: int = 100,
        sort_by: str = "relevance",
        result_type: str = "core"
    ) -> Optional[Dict[str, Any]]:
        """
        Search Europe PMC

        Args:
            query: Search query with filters
            max_results: Maximum number of results
            sort_by: Sort order ("relevance", "date", "cited")
            result_type: "core" (metadata only) or "lite" (minimal)

        Returns:
            Search results dictionary or None on error
        """
        params = {
            "query": query,
            "pageSize": min(max_results, 1000),  # Europe PMC max is 1000
            "format": "json",
            "resultType": result_type
        }

        # Add sort parameter
        sort_map = {
            "relevance": "relevance",
            "date": "P_PDATE_D desc",
            "cited": "CITED desc"
        }
        if sort_by in sort_map:
            params["sort"] = sort_map[sort_by]

        # Add email if available
        if self.config.get("email"):
            params["email"] = self.config["email"]

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.SEARCH_URL, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        logger.error(f"Europe PMC search failed: {response.status} - {error_text}")
                        return None

                    data = await response.json()
                    result_list = data.get("resultList", {})
                    results = result_list.get("result", [])
                    hit_count = int(data.get("hitCount", 0))

                    logger.info(f"Europe PMC: Found {hit_count} articles for query: {query[:50]}...")
                    logger.info(f"Europe PMC: Returning {len(results)} articles")

                    return {
                        "results": results,
                        "hit_count": hit_count,
                        "total_returned": len(results)
                    }

        except asyncio.TimeoutError:
            logger.error(f"Europe PMC: Timeout searching for: {query}")
            return None
        except Exception as e:
            logger.error(f"Europe PMC: Error searching for {query}: {e}")
            return None

    async def _fetch_annotations_async(
        self,
        article_ids: List[str],
        source: str = "MED"
    ) -> Optional[Dict[str, Any]]:
        """
        Fetch chemical and other annotations for articles

        Args:
            article_ids: List of article IDs (PMIDs or PMCIDs)
            source: Source database ("MED" for PubMed, "PMC" for PMC)

        Returns:
            Annotations dictionary or None on error
        """
        if not article_ids:
            return {}

        # Build article ID string
        id_params = [f"{source}:{aid}" for aid in article_ids[:100]]  # Limit to 100 per request

        params = {
            "articleIds": ",".join(id_params),
            "type": "Chemicals",  # Can also be "Diseases", "Organisms", "Gene_Proteins"
            "format": "JSON"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    self.ANNOTATIONS_URL,
                    data=params,
                    timeout=timeout
                ) as response:
                    if response.status != 200:
                        logger.warning(f"Europe PMC annotations fetch failed: {response.status}")
                        return None

                    data = await response.json()
                    logger.info(f"Europe PMC: Fetched annotations for {len(article_ids)} articles")

                    return data

        except asyncio.TimeoutError:
            logger.error("Europe PMC: Timeout fetching annotations")
            return None
        except Exception as e:
            logger.error(f"Europe PMC: Error fetching annotations: {e}")
            return None

    async def _fetch_fulltext_async(
        self,
        article_id: str,
        source: str = "MED"
    ) -> Optional[str]:
        """
        Fetch full-text XML for an article

        Args:
            article_id: Article ID (PMID or PMCID)
            source: Source database ("MED" or "PMC")

        Returns:
            Full-text XML string or None if not available
        """
        url = self.FULLTEXT_URL.format(source=source, id=article_id)

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=timeout) as response:
                    if response.status == 200:
                        fulltext = await response.text()
                        logger.info(f"Europe PMC: Retrieved full text for {article_id}")
                        return fulltext
                    elif response.status == 404:
                        logger.info(f"Europe PMC: Full text not available for {article_id}")
                        return None
                    else:
                        logger.warning(f"Europe PMC: Full text fetch failed for {article_id}: {response.status}")
                        return None

        except asyncio.TimeoutError:
            logger.error(f"Europe PMC: Timeout fetching full text for {article_id}")
            return None
        except Exception as e:
            logger.error(f"Europe PMC: Error fetching full text for {article_id}: {e}")
            return None

    def _parse_article_results(self, results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Parse and normalize article results

        Args:
            results: Raw article results from Europe PMC

        Returns:
            List of normalized article dictionaries
        """
        articles = []

        for result in results:
            try:
                article = {
                    "id": result.get("id"),
                    "source": result.get("source"),
                    "pmid": result.get("pmid"),
                    "pmcid": result.get("pmcid"),
                    "doi": result.get("doi"),
                    "title": result.get("title"),
                    "abstract": result.get("abstractText"),
                    "authors": self._parse_authors(result.get("authorString", "")),
                    "journal": {
                        "title": result.get("journalTitle"),
                        "volume": result.get("journalVolume"),
                        "issue": result.get("issue"),
                        "pages": result.get("pageInfo")
                    },
                    "publication_date": result.get("firstPublicationDate"),
                    "pub_year": result.get("pubYear"),
                    "pub_type": result.get("pubType"),
                    "is_open_access": result.get("isOpenAccess") == "Y",
                    "has_fulltext": result.get("hasTextMinedTerms") == "Y",
                    "citations": {
                        "cited_by_count": int(result.get("citedByCount", 0))
                    },
                    "keywords": result.get("keywordList", {}).get("keyword", []) if result.get("keywordList") else [],
                    "mesh_terms": self._extract_mesh_terms(result),
                    "chemicals": [],  # Will be populated by annotations
                    "url": f"https://europepmc.org/article/{result.get('source')}/{result.get('id')}"
                }

                articles.append(article)

            except Exception as e:
                logger.warning(f"Error parsing article result: {e}")
                continue

        return articles

    def _parse_authors(self, author_string: str) -> List[str]:
        """
        Parse author string into list

        Args:
            author_string: Semicolon or comma-separated author names

        Returns:
            List of author names
        """
        if not author_string:
            return []

        # Try semicolon first, then comma
        if ";" in author_string:
            return [a.strip() for a in author_string.split(";") if a.strip()]
        else:
            return [a.strip() for a in author_string.split(",") if a.strip()]

    def _extract_mesh_terms(self, result: Dict[str, Any]) -> List[Dict[str, str]]:
        """
        Extract MeSH terms from article result

        Args:
            result: Article result dictionary

        Returns:
            List of MeSH term dictionaries
        """
        mesh_terms = []

        mesh_heading_list = result.get("meshHeadingList", {})
        if mesh_heading_list:
            for heading in mesh_heading_list.get("meshHeading", []):
                if isinstance(heading, dict):
                    descriptor = heading.get("descriptorName")
                    if descriptor:
                        mesh_terms.append({
                            "term": descriptor,
                            "major_topic": heading.get("majorTopic_YN") == "Y"
                        })

        return mesh_terms

    def _merge_annotations(
        self,
        articles: List[Dict[str, Any]],
        annotations: Optional[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Merge chemical annotations into article data

        Args:
            articles: List of article dictionaries
            annotations: Annotations data from API

        Returns:
            Articles with merged annotations
        """
        if not annotations:
            return articles

        # Build annotation lookup by article ID
        annotation_lookup = {}
        for annotation_group in annotations:
            if isinstance(annotation_group, dict):
                ext_id = annotation_group.get("extId")
                annotations_list = annotation_group.get("annotations", [])
                if ext_id:
                    annotation_lookup[ext_id] = annotations_list

        # Merge annotations into articles
        for article in articles:
            article_id = article.get("id")
            if article_id in annotation_lookup:
                chemicals = []
                for ann in annotation_lookup[article_id]:
                    if isinstance(ann, dict):
                        chemicals.append({
                            "name": ann.get("exact"),
                            "type": ann.get("type"),
                            "database": ann.get("provider"),
                            "id": ann.get("id")
                        })
                article["chemicals"] = chemicals

        return articles

    def _summarize_results(self, articles: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Create summary statistics from article results

        Args:
            articles: List of article dictionaries

        Returns:
            Summary dictionary
        """
        if not articles:
            return {
                "total_articles": 0,
                "open_access_count": 0,
                "with_fulltext_count": 0,
                "common_chemicals": [],
                "top_journals": [],
                "date_range": None
            }

        # Count open access and full text
        open_access_count = sum(1 for a in articles if a.get("is_open_access"))
        fulltext_count = sum(1 for a in articles if a.get("has_fulltext"))

        # Extract common chemicals
        chemical_counts = {}
        for article in articles:
            for chemical in article.get("chemicals", []):
                name = chemical.get("name")
                if name:
                    chemical_counts[name] = chemical_counts.get(name, 0) + 1

        common_chemicals = sorted(chemical_counts.items(), key=lambda x: x[1], reverse=True)[:20]

        # Extract top journals
        journal_counts = {}
        for article in articles:
            journal = article.get("journal", {}).get("title")
            if journal:
                journal_counts[journal] = journal_counts.get(journal, 0) + 1

        top_journals = sorted(journal_counts.items(), key=lambda x: x[1], reverse=True)[:5]

        # Date range
        dates = [a.get("publication_date") for a in articles if a.get("publication_date")]
        date_range = None
        if dates:
            dates_sorted = sorted(dates)
            date_range = {"earliest": dates_sorted[0], "latest": dates_sorted[-1]}

        # Citation statistics
        citation_counts = [a.get("citations", {}).get("cited_by_count", 0) for a in articles]

        return {
            "total_articles": len(articles),
            "open_access_count": open_access_count,
            "open_access_percentage": (open_access_count / len(articles) * 100) if articles else 0,
            "with_fulltext_count": fulltext_count,
            "common_chemicals": [{"name": c, "count": count} for c, count in common_chemicals],
            "top_journals": [{"journal": j, "count": c} for j, c in top_journals],
            "date_range": date_range,
            "avg_citations": sum(citation_counts) / len(citation_counts) if citation_counts else 0,
            "max_citations": max(citation_counts) if citation_counts else 0
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Europe PMC article search

        Args:
            input_data: Search query string or dict with parameters
                       String: "aspirin AND cancer"
                       Dict: {
                           "query": "aspirin",
                           "max_results": 50,
                           "filters": {
                               "open_access": True,
                               "has_chemicals": True,
                               "date_from": "2020-01-01"
                           }
                       }
            **kwargs: Additional parameters
                     - max_results: Maximum number of results (default: 100)
                     - sort_by: Sort order ("relevance", "date", "cited")
                     - filters: Dict with various filters
                     - fetch_annotations: Whether to fetch chemical annotations (default: True)
                     - fetch_fulltext: Whether to fetch full text for results (default: False)

        Returns:
            AdapterResult containing article information with chemical entities
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
            base_query = input_data
            max_results = kwargs.get("max_results", self.config.get("max_results", 100))
            sort_by = kwargs.get("sort_by", "relevance")
            filters = kwargs.get("filters")
        else:
            base_query = input_data["query"]
            max_results = input_data.get("max_results", kwargs.get("max_results", 100))
            sort_by = input_data.get("sort_by", kwargs.get("sort_by", "relevance"))
            filters = input_data.get("filters", kwargs.get("filters"))

        fetch_annotations = kwargs.get("fetch_annotations", True)
        fetch_fulltext = kwargs.get("fetch_fulltext", False)

        # Build complete query string
        query = self._build_query_string(base_query, filters)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.2)
        await asyncio.sleep(rate_delay)

        # Search Europe PMC
        search_results = await self._search_europepmc_async(query, max_results, sort_by)

        if search_results is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search Europe PMC",
                metadata={"source": "europepmc", "query": query}
            )

        # Parse article results
        articles = self._parse_article_results(search_results["results"])

        # Fetch chemical annotations if requested
        if fetch_annotations and articles:
            await asyncio.sleep(rate_delay)
            article_ids = [a["id"] for a in articles if a.get("id")]
            # Determine source (prefer PMC for better annotations)
            source = "PMC" if any(a.get("pmcid") for a in articles) else "MED"
            annotations = await self._fetch_annotations_async(article_ids, source)
            if annotations:
                articles = self._merge_annotations(articles, annotations)

        # Fetch full text for selected articles if requested
        if fetch_fulltext and articles:
            for article in articles[:5]:  # Limit to first 5 to avoid excessive requests
                await asyncio.sleep(rate_delay)
                article_id = article.get("pmcid") or article.get("pmid")
                source = "PMC" if article.get("pmcid") else "MED"
                if article_id:
                    fulltext = await self._fetch_fulltext_async(article_id, source)
                    if fulltext:
                        article["fulltext_xml"] = fulltext

        # Create summary
        summary = self._summarize_results(articles)

        result_data = {
            "articles": articles,
            "summary": summary,
            "query": query,
            "hit_count": search_results["hit_count"],
            "total_returned": search_results["total_returned"]
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "europepmc",
                "query": query,
                "adapter_version": self.version,
                "total_results": search_results["hit_count"],
                "articles_returned": len(articles)
            }
        )
