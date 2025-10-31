"""
Lens.org Adapter - Unified search for patents and scholarly works
Uses Lens.org Scholarly API for patent and paper search with citation analysis
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from urllib.parse import quote

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class LensAdapter(AdapterProtocol):
    """
    Adapter for Lens.org Scholarly API

    Capabilities:
    - Patent search with chemistry filters
    - Scholarly article search
    - Citation network analysis
    - Patent family information
    - Chemical structure search (when available)
    - Cross-referencing between patents and papers

    Requires API token from: https://www.lens.org/lens/user/subscriptions#scholar
    Free academic access available.
    """

    PATENT_SEARCH_URL = "https://api.lens.org/patent/search"
    SCHOLARLY_SEARCH_URL = "https://api.lens.org/scholarly/search"

    def __init__(self, api_token: Optional[str] = None):
        """
        Initialize Lens.org adapter

        Args:
            api_token: Lens.org API token (required for API access)
        """
        super().__init__(
            name="lens",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.5,  # Conservative rate limiting
                "timeout": 30,
                "max_results": 100,
                "api_token": api_token
            }
        )
        self.version = "1.0.0"

        if not api_token:
            logger.warning("Lens.org: No API token provided. Set via config or environment variable.")

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

    def _build_lens_query(
        self,
        query: str,
        search_type: str = "patent",
        filters: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Build Lens.org API query structure

        Args:
            query: Search query text
            search_type: "patent" or "scholarly"
            filters: Optional filters (date_published, jurisdiction, etc.)

        Returns:
            Query dictionary for Lens.org API
        """
        # Build the query structure
        lens_query = {
            "query": {
                "match": {
                    "_all": query
                }
            },
            "size": self.config.get("max_results", 100),
            "sort": [{"relevance": "desc"}]
        }

        # Add filters if provided
        if filters:
            filter_clauses = []

            # Date filters
            if "date_from" in filters or "date_to" in filters:
                date_field = "date_published" if search_type == "patent" else "year_published"
                date_range = {}
                if "date_from" in filters:
                    date_range["gte"] = filters["date_from"]
                if "date_to" in filters:
                    date_range["lte"] = filters["date_to"]

                filter_clauses.append({
                    "range": {date_field: date_range}
                })

            # Patent-specific filters
            if search_type == "patent":
                if "jurisdiction" in filters:
                    filter_clauses.append({
                        "term": {"jurisdiction": filters["jurisdiction"]}
                    })

                if "patent_type" in filters:
                    filter_clauses.append({
                        "term": {"type": filters["patent_type"]}
                    })

                if "applicant" in filters:
                    filter_clauses.append({
                        "match": {"applicant": filters["applicant"]}
                    })

            # Scholarly-specific filters
            if search_type == "scholarly":
                if "open_access" in filters and filters["open_access"]:
                    filter_clauses.append({
                        "term": {"open_access": True}
                    })

                if "source_type" in filters:
                    filter_clauses.append({
                        "term": {"source.type": filters["source_type"]}
                    })

            if filter_clauses:
                lens_query["query"] = {
                    "bool": {
                        "must": lens_query["query"],
                        "filter": filter_clauses
                    }
                }

        return lens_query

    async def _search_patents_async(
        self,
        query: str,
        filters: Optional[Dict[str, Any]] = None
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Search Lens.org patents

        Args:
            query: Search query
            filters: Optional filters

        Returns:
            List of patent records or None on error
        """
        if not self.config.get("api_token"):
            logger.error("Lens.org: API token required for patent search")
            return None

        lens_query = self._build_lens_query(query, "patent", filters)

        headers = {
            "Authorization": f"Bearer {self.config['api_token']}",
            "Content-Type": "application/json",
            "User-Agent": "PharmForge/1.0"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    self.PATENT_SEARCH_URL,
                    json=lens_query,
                    headers=headers,
                    timeout=timeout
                ) as response:
                    if response.status == 401:
                        logger.error("Lens.org: Invalid or missing API token")
                        return None
                    elif response.status == 429:
                        logger.error("Lens.org: Rate limit exceeded")
                        return None
                    elif response.status != 200:
                        error_text = await response.text()
                        logger.error(f"Lens.org patent search failed: {response.status} - {error_text}")
                        return None

                    data = await response.json()
                    patents = data.get("data", [])
                    total = data.get("total", 0)

                    logger.info(f"Lens.org: Found {total} patents, returning {len(patents)}")

                    return self._parse_patent_results(patents)

        except asyncio.TimeoutError:
            logger.error(f"Lens.org: Timeout searching patents for: {query}")
            return None
        except Exception as e:
            logger.error(f"Lens.org: Error searching patents: {e}")
            return None

    async def _search_scholarly_async(
        self,
        query: str,
        filters: Optional[Dict[str, Any]] = None
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Search Lens.org scholarly works

        Args:
            query: Search query
            filters: Optional filters

        Returns:
            List of scholarly records or None on error
        """
        if not self.config.get("api_token"):
            logger.error("Lens.org: API token required for scholarly search")
            return None

        lens_query = self._build_lens_query(query, "scholarly", filters)

        headers = {
            "Authorization": f"Bearer {self.config['api_token']}",
            "Content-Type": "application/json",
            "User-Agent": "PharmForge/1.0"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.post(
                    self.SCHOLARLY_SEARCH_URL,
                    json=lens_query,
                    headers=headers,
                    timeout=timeout
                ) as response:
                    if response.status == 401:
                        logger.error("Lens.org: Invalid or missing API token")
                        return None
                    elif response.status == 429:
                        logger.error("Lens.org: Rate limit exceeded")
                        return None
                    elif response.status != 200:
                        error_text = await response.text()
                        logger.error(f"Lens.org scholarly search failed: {response.status} - {error_text}")
                        return None

                    data = await response.json()
                    works = data.get("data", [])
                    total = data.get("total", 0)

                    logger.info(f"Lens.org: Found {total} scholarly works, returning {len(works)}")

                    return self._parse_scholarly_results(works)

        except asyncio.TimeoutError:
            logger.error(f"Lens.org: Timeout searching scholarly works for: {query}")
            return None
        except Exception as e:
            logger.error(f"Lens.org: Error searching scholarly works: {e}")
            return None

    def _parse_patent_results(self, patents: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Parse and normalize patent results from Lens.org

        Args:
            patents: Raw patent data from API

        Returns:
            List of normalized patent dictionaries
        """
        parsed_patents = []

        for patent in patents:
            try:
                parsed = {
                    "lens_id": patent.get("lens_id"),
                    "patent_number": patent.get("publication_number"),
                    "title": patent.get("title"),
                    "abstract": patent.get("abstract"),
                    "applicants": patent.get("applicant", []),
                    "inventors": patent.get("inventor", []),
                    "publication_date": patent.get("date_published"),
                    "filing_date": patent.get("filing_date"),
                    "jurisdiction": patent.get("jurisdiction"),
                    "patent_type": patent.get("type"),
                    "family_id": patent.get("simple_family_id"),
                    "citations": {
                        "cited_by_count": patent.get("cited_by_patent_count", 0),
                        "citing_count": patent.get("citing_patent_count", 0),
                        "backward_citations": patent.get("backward_citation_count", 0)
                    },
                    "classifications": {
                        "cpc": patent.get("cpc_classifications", []),
                        "ipc": patent.get("ipc_classifications", [])
                    },
                    "url": f"https://www.lens.org/lens/patent/{patent.get('lens_id')}"
                }

                # Extract chemical entities if available
                if "chemical_entities" in patent:
                    parsed["chemicals"] = patent["chemical_entities"]

                parsed_patents.append(parsed)

            except Exception as e:
                logger.warning(f"Error parsing patent result: {e}")
                continue

        return parsed_patents

    def _parse_scholarly_results(self, works: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Parse and normalize scholarly results from Lens.org

        Args:
            works: Raw scholarly data from API

        Returns:
            List of normalized scholarly dictionaries
        """
        parsed_works = []

        for work in works:
            try:
                # Extract authors
                authors = []
                for author in work.get("authors", []):
                    if isinstance(author, dict):
                        name_parts = []
                        if "first_name" in author:
                            name_parts.append(author["first_name"])
                        if "last_name" in author:
                            name_parts.append(author["last_name"])
                        if name_parts:
                            authors.append(" ".join(name_parts))

                parsed = {
                    "lens_id": work.get("lens_id"),
                    "title": work.get("title"),
                    "abstract": work.get("abstract"),
                    "authors": authors,
                    "year_published": work.get("year_published"),
                    "date_published": work.get("date_published"),
                    "source": {
                        "title": work.get("source", {}).get("title"),
                        "type": work.get("source", {}).get("type"),
                        "publisher": work.get("source", {}).get("publisher")
                    },
                    "citations": {
                        "cited_by_count": work.get("scholarly_citations_count", 0),
                        "references_count": work.get("references_count", 0)
                    },
                    "identifiers": {
                        "doi": work.get("doi"),
                        "pmid": work.get("pubmed_id"),
                        "issn": work.get("issn")
                    },
                    "open_access": work.get("open_access", False),
                    "url": f"https://www.lens.org/lens/scholar/{work.get('lens_id')}"
                }

                # Extract fields of study/keywords
                if "fields_of_study" in work:
                    parsed["fields_of_study"] = work["fields_of_study"]

                # Extract chemical mentions if available
                if "chemical_entities" in work:
                    parsed["chemicals"] = work["chemical_entities"]

                parsed_works.append(parsed)

            except Exception as e:
                logger.warning(f"Error parsing scholarly result: {e}")
                continue

        return parsed_works

    def _build_citation_network(
        self,
        items: List[Dict[str, Any]],
        item_type: str = "patent"
    ) -> Dict[str, Any]:
        """
        Build citation network from results

        Args:
            items: List of patents or scholarly works
            item_type: "patent" or "scholarly"

        Returns:
            Citation network dictionary
        """
        nodes = []
        edges = []

        for item in items:
            node = {
                "id": item.get("lens_id"),
                "title": item.get("title"),
                "type": item_type,
                "cited_by_count": item.get("citations", {}).get("cited_by_count", 0)
            }
            nodes.append(node)

            # Add citation edges (if citation details are available)
            # Note: Full citation network requires additional API calls
            # This is a simplified version

        return {
            "nodes": nodes,
            "edges": edges,
            "stats": {
                "total_nodes": len(nodes),
                "total_edges": len(edges),
                "avg_citations": sum(n["cited_by_count"] for n in nodes) / len(nodes) if nodes else 0
            }
        }

    def _summarize_results(
        self,
        patents: Optional[List[Dict[str, Any]]] = None,
        scholarly: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, Any]:
        """
        Create summary statistics from results

        Args:
            patents: List of patent results
            scholarly: List of scholarly results

        Returns:
            Summary dictionary
        """
        summary = {}

        if patents:
            # Top applicants
            applicant_counts = {}
            for patent in patents:
                for applicant in patent.get("applicants", []):
                    name = applicant if isinstance(applicant, str) else applicant.get("name", "")
                    if name:
                        applicant_counts[name] = applicant_counts.get(name, 0) + 1

            top_applicants = sorted(applicant_counts.items(), key=lambda x: x[1], reverse=True)[:5]

            # Citation statistics
            citation_counts = [p.get("citations", {}).get("cited_by_count", 0) for p in patents]

            summary["patents"] = {
                "total": len(patents),
                "top_applicants": [{"name": a, "count": c} for a, c in top_applicants],
                "avg_citations": sum(citation_counts) / len(citation_counts) if citation_counts else 0,
                "max_citations": max(citation_counts) if citation_counts else 0
            }

        if scholarly:
            # Top sources
            source_counts = {}
            for work in scholarly:
                source = work.get("source", {}).get("title")
                if source:
                    source_counts[source] = source_counts.get(source, 0) + 1

            top_sources = sorted(source_counts.items(), key=lambda x: x[1], reverse=True)[:5]

            # Citation statistics
            citation_counts = [w.get("citations", {}).get("cited_by_count", 0) for w in scholarly]

            # Open access count
            open_access_count = sum(1 for w in scholarly if w.get("open_access"))

            summary["scholarly"] = {
                "total": len(scholarly),
                "top_sources": [{"source": s, "count": c} for s, c in top_sources],
                "avg_citations": sum(citation_counts) / len(citation_counts) if citation_counts else 0,
                "max_citations": max(citation_counts) if citation_counts else 0,
                "open_access_count": open_access_count,
                "open_access_percentage": (open_access_count / len(scholarly) * 100) if scholarly else 0
            }

        return summary

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute Lens.org search

        Args:
            input_data: Search query string or dict with parameters
                       String: "aspirin pharmaceutical"
                       Dict: {
                           "query": "aspirin",
                           "search_type": "both",  # "patent", "scholarly", or "both"
                           "filters": {
                               "date_from": "2020-01-01",
                               "jurisdiction": "US"
                           }
                       }
            **kwargs: Additional parameters
                     - search_type: "patent", "scholarly", or "both" (default: "both")
                     - filters: Dict with various filters
                     - build_citation_network: Whether to build citation network (default: False)

        Returns:
            AdapterResult containing patent and/or scholarly information
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid search query"
            )

        # Check for API token
        if not self.config.get("api_token"):
            return AdapterResult(
                success=False,
                data=None,
                error="Lens.org API token required. Get one at https://www.lens.org/lens/user/subscriptions#scholar"
            )

        # Parse input
        if isinstance(input_data, str):
            query = input_data
            search_type = kwargs.get("search_type", "both")
            filters = kwargs.get("filters")
        else:
            query = input_data["query"]
            search_type = input_data.get("search_type", kwargs.get("search_type", "both"))
            filters = input_data.get("filters", kwargs.get("filters"))

        build_citation_network = kwargs.get("build_citation_network", False)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.5)
        await asyncio.sleep(rate_delay)

        # Execute searches based on type
        patents = None
        scholarly = None

        if search_type in ["patent", "both"]:
            patents = await self._search_patents_async(query, filters)
            if patents is None and search_type == "patent":
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to search Lens.org patents",
                    metadata={"source": "lens", "query": query}
                )

        if search_type in ["scholarly", "both"]:
            if search_type == "both":
                await asyncio.sleep(rate_delay)  # Rate limit between searches
            scholarly = await self._search_scholarly_async(query, filters)
            if scholarly is None and search_type == "scholarly":
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to search Lens.org scholarly works",
                    metadata={"source": "lens", "query": query}
                )

        # Build citation network if requested
        citation_network = None
        if build_citation_network:
            if patents:
                patent_network = self._build_citation_network(patents, "patent")
                citation_network = {"patent_network": patent_network}
            if scholarly:
                scholarly_network = self._build_citation_network(scholarly, "scholarly")
                if citation_network:
                    citation_network["scholarly_network"] = scholarly_network
                else:
                    citation_network = {"scholarly_network": scholarly_network}

        # Create summary
        summary = self._summarize_results(patents, scholarly)

        result_data = {
            "patents": patents or [],
            "scholarly": scholarly or [],
            "summary": summary,
            "citation_network": citation_network,
            "query": query,
            "search_type": search_type
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "lens",
                "query": query,
                "adapter_version": self.version,
                "patent_count": len(patents) if patents else 0,
                "scholarly_count": len(scholarly) if scholarly else 0
            }
        )
