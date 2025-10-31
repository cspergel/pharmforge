"""
GEO Adapter - Fetches gene expression datasets from NCBI Gene Expression Omnibus
Uses E-utilities API for dataset search, retrieval, and metadata access
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
import xml.etree.ElementTree as ET
from urllib.parse import quote
import re

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class GEOAdapter(AdapterProtocol):
    """
    Adapter for NCBI Gene Expression Omnibus (GEO) E-utilities API
    Searches and retrieves gene expression datasets

    Capabilities:
    - Search datasets by GEO accession (GSE, GDS)
    - Search by gene name/symbol
    - Search by disease/condition/keyword
    - Retrieve dataset metadata
    - Access dataset summaries
    - Link to downloadable data files

    GEO databases:
    - gds: GEO DataSets (curated datasets)
    - geoprofiles: GEO Profiles (gene expression profiles)
    """

    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"

    # FTP base URL for downloading full datasets
    FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo"

    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        """
        Initialize GEO adapter

        Args:
            email: Email for NCBI (recommended by NCBI guidelines)
            api_key: NCBI API key for higher rate limits (10 req/s vs 3 req/s)
        """
        super().__init__(
            name="geo",
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
        Validate that input is a valid search query or GEO accession

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

    def _is_geo_accession(self, query: str) -> Optional[str]:
        """
        Check if query is a GEO accession and return its type

        Args:
            query: Query string

        Returns:
            Accession type ("GSE", "GDS", "GPL", "GSM") or None
        """
        query = query.strip().upper()
        # Check for GEO accession patterns
        if re.match(r'^GSE\d+$', query):
            return "GSE"  # Series
        elif re.match(r'^GDS\d+$', query):
            return "GDS"  # Dataset
        elif re.match(r'^GPL\d+$', query):
            return "GPL"  # Platform
        elif re.match(r'^GSM\d+$', query):
            return "GSM"  # Sample
        return None

    def _build_ftp_url(self, accession: str, file_type: str = "soft") -> str:
        """
        Build FTP URL for downloading GEO data

        Args:
            accession: GEO accession (e.g., GSE1234)
            file_type: File type ("soft", "matrix", "suppl")

        Returns:
            FTP URL string
        """
        acc_type = self._is_geo_accession(accession)
        if not acc_type:
            return None

        # Extract number and create nnn format
        acc_num = re.search(r'\d+', accession).group()
        acc_base = accession[:3]  # GSE, GDS, GPL, GSM

        # Create nnn format (e.g., 1234 -> 1nnn)
        if len(acc_num) >= 3:
            nnn_format = acc_num[:-3] + "nnn"
        else:
            nnn_format = "nnn"

        if acc_type == "GSE":
            if file_type == "soft":
                return f"{self.FTP_BASE}/series/{acc_base}{nnn_format}/{accession}/soft/{accession}_family.soft.gz"
            elif file_type == "matrix":
                return f"{self.FTP_BASE}/series/{acc_base}{nnn_format}/{accession}/matrix/{accession}_series_matrix.txt.gz"
            elif file_type == "suppl":
                return f"{self.FTP_BASE}/series/{acc_base}{nnn_format}/{accession}/suppl/"
        elif acc_type == "GDS":
            return f"{self.FTP_BASE}/datasets/{acc_base}{nnn_format}/{accession}/soft/{accession}.soft.gz"

        return None

    async def _search_geo_async(
        self,
        query: str,
        database: str = "gds",
        max_results: int = 100,
        filters: Optional[Dict[str, Any]] = None
    ) -> Optional[Dict[str, Any]]:
        """
        Search GEO and return UIDs and metadata

        Args:
            query: Search query (supports boolean operators: AND, OR, NOT)
            database: Database to search ("gds" for datasets, "geoprofiles" for profiles)
            max_results: Maximum number of results to return
            filters: Optional filters (organism, entry_type, date_from, date_to)

        Returns:
            Dict with UIDs and search info or None on error
        """
        # Build query with filters
        search_terms = [query]

        if filters:
            if "organism" in filters:
                search_terms.append(f"{filters['organism']}[ORGN]")
            if "entry_type" in filters:
                search_terms.append(f"{filters['entry_type']}[ETYP]")
            if "platform" in filters:
                search_terms.append(f"{filters['platform']}[Platform]")

        final_query = " AND ".join(search_terms)

        params = {
            "db": database,
            "term": final_query,
            "retmax": max_results,
            "retmode": "json",
            "email": self.config.get("email")
        }

        # Add API key if available
        if self.config.get("api_key"):
            params["api_key"] = self.config["api_key"]

        # Add date filters if provided
        if filters:
            if "date_from" in filters:
                params["mindate"] = filters["date_from"]
            if "date_to" in filters:
                params["maxdate"] = filters["date_to"]
            if "date_from" in filters or "date_to" in filters:
                params["datetype"] = "pdat"

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.ESEARCH_URL, params=params, timeout=timeout) as response:
                    if response.status != 200:
                        logger.error(f"GEO search failed with status {response.status}")
                        return None

                    data = await response.json()
                    result = data.get("esearchresult", {})
                    uids = result.get("idlist", [])
                    count = int(result.get("count", 0))

                    logger.info(f"GEO: Found {count} datasets for query: {query[:50]}...")
                    logger.info(f"GEO: Returning {len(uids)} UIDs")

                    return {
                        "uids": uids,
                        "count": count,
                        "query": final_query,
                        "database": database
                    }

        except asyncio.TimeoutError:
            logger.error(f"GEO: Timeout searching for: {query}")
            return None
        except Exception as e:
            logger.error(f"GEO: Error searching for {query}: {e}")
            return None

    async def _fetch_summaries_async(
        self,
        uids: List[str],
        database: str = "gds"
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Fetch dataset summaries for given UIDs

        Args:
            uids: List of GEO UIDs
            database: Database type ("gds" or "geoprofiles")

        Returns:
            List of dataset summaries or None on error
        """
        if not uids:
            return []

        # Fetch in batches of 200 (NCBI limit)
        batch_size = 200
        all_summaries = []

        for i in range(0, len(uids), batch_size):
            batch_uids = uids[i:i + batch_size]
            uid_str = ",".join(batch_uids)

            params = {
                "db": database,
                "id": uid_str,
                "retmode": "json",
                "email": self.config.get("email")
            }

            if self.config.get("api_key"):
                params["api_key"] = self.config["api_key"]

            timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 30))

            try:
                async with aiohttp.ClientSession() as session:
                    async with session.get(self.ESUMMARY_URL, params=params, timeout=timeout) as response:
                        if response.status != 200:
                            logger.error(f"GEO summary fetch failed with status {response.status}")
                            continue

                        data = await response.json()
                        result = data.get("result", {})

                        # Extract summaries for each UID
                        for uid in batch_uids:
                            if uid in result:
                                summary = self._parse_summary(result[uid], database)
                                if summary:
                                    all_summaries.append(summary)

                        # Rate limiting between batches
                        if i + batch_size < len(uids):
                            await asyncio.sleep(self.config.get("rate_limit_delay", 0.34))

            except asyncio.TimeoutError:
                logger.error(f"GEO: Timeout fetching summaries for batch starting at {i}")
                continue
            except Exception as e:
                logger.error(f"GEO: Error fetching summaries: {e}")
                continue

        return all_summaries if all_summaries else None

    def _parse_summary(self, summary_data: Dict[str, Any], database: str) -> Optional[Dict[str, Any]]:
        """
        Parse summary data from eSummary response

        Args:
            summary_data: Raw summary data from NCBI
            database: Database type

        Returns:
            Parsed summary dictionary
        """
        try:
            parsed = {
                "uid": summary_data.get("uid"),
                "database": database
            }

            if database == "gds":
                # GEO DataSets specific fields
                parsed.update({
                    "accession": summary_data.get("accession"),
                    "title": summary_data.get("title"),
                    "summary": summary_data.get("summary"),
                    "gpl": summary_data.get("gpl"),  # Platform ID
                    "gse": summary_data.get("gse"),  # Series ID
                    "taxon": summary_data.get("taxon"),
                    "entrytype": summary_data.get("entrytype"),
                    "gdstype": summary_data.get("gdstype"),
                    "ptechtype": summary_data.get("ptechtype"),  # Platform technology type
                    "valtype": summary_data.get("valtype"),  # Value type
                    "n_samples": summary_data.get("n_samples"),
                    "pubmedids": summary_data.get("pubmedids", []),
                    "projects": summary_data.get("projects", []),
                    "ftplink": summary_data.get("ftplink"),
                })
            elif database == "geoprofiles":
                # GEO Profiles specific fields
                # Note: geoprofiles may have different field names in API response
                parsed.update({
                    "accession": summary_data.get("accession") or summary_data.get("name"),
                    "title": summary_data.get("title") or summary_data.get("description"),
                    "organism": summary_data.get("organism") or summary_data.get("taxon"),
                    "genesymbol": summary_data.get("genesymbol") or summary_data.get("symbol"),
                    "genefullname": summary_data.get("genefullname") or summary_data.get("name"),
                })

            return parsed

        except Exception as e:
            logger.warning(f"Error parsing summary: {e}")
            return None

    async def _fetch_by_accession_async(
        self,
        accession: str
    ) -> Optional[Dict[str, Any]]:
        """
        Fetch dataset by GEO accession directly

        Args:
            accession: GEO accession (GSE, GDS, etc.)

        Returns:
            Dataset information or None
        """
        acc_type = self._is_geo_accession(accession)
        if not acc_type:
            logger.warning(f"Invalid GEO accession format: {accession}")
            return None

        # Determine which database to use
        if acc_type in ["GSE", "GDS"]:
            database = "gds"
            query = f"{accession}[ACCN]"
        else:
            # For GPL and GSM, we might need different handling
            return {"error": f"Direct fetch for {acc_type} not fully implemented. Use search instead."}

        # Search for the specific accession
        search_result = await self._search_geo_async(query, database=database, max_results=1)

        if not search_result or not search_result.get("uids"):
            return None

        # Fetch summary
        summaries = await self._fetch_summaries_async(search_result["uids"], database=database)

        if summaries and len(summaries) > 0:
            dataset = summaries[0]

            # Add download URLs
            dataset["download_urls"] = {
                "soft": self._build_ftp_url(accession, "soft"),
                "matrix": self._build_ftp_url(accession, "matrix"),
                "supplementary": self._build_ftp_url(accession, "suppl")
            }

            return dataset

        return None

    def _summarize_results(self, datasets: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Create summary statistics from dataset results

        Args:
            datasets: List of dataset dictionaries

        Returns:
            Summary dictionary with statistics
        """
        if not datasets:
            return {
                "total_datasets": 0,
                "organisms": [],
                "platform_types": [],
                "entry_types": []
            }

        # Count organisms
        organism_counts = {}
        for ds in datasets:
            taxon = ds.get("taxon")
            if taxon:
                organism_counts[taxon] = organism_counts.get(taxon, 0) + 1

        # Count platform types
        platform_counts = {}
        for ds in datasets:
            ptech = ds.get("ptechtype")
            if ptech:
                platform_counts[ptech] = platform_counts.get(ptech, 0) + 1

        # Count entry types
        entry_counts = {}
        for ds in datasets:
            etype = ds.get("entrytype")
            if etype:
                entry_counts[etype] = entry_counts.get(etype, 0) + 1

        # Count datasets with PubMed links
        with_pubmed = sum(1 for ds in datasets if ds.get("pubmedids"))

        return {
            "total_datasets": len(datasets),
            "datasets_with_pubmed": with_pubmed,
            "organisms": [{"organism": org, "count": count}
                         for org, count in sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)],
            "platform_types": [{"type": ptype, "count": count}
                              for ptype, count in sorted(platform_counts.items(), key=lambda x: x[1], reverse=True)],
            "entry_types": [{"type": etype, "count": count}
                           for etype, count in sorted(entry_counts.items(), key=lambda x: x[1], reverse=True)]
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute GEO dataset search or retrieval

        Args:
            input_data: Search query string or dict with parameters
                       String examples:
                       - "breast cancer" (keyword search)
                       - "TP53" (gene search)
                       - "GSE1234" (accession search)
                       Dict: {"query": "diabetes", "organism": "Homo sapiens", "max_results": 50}
            **kwargs: Additional parameters
                     - database: "gds" (datasets) or "geoprofiles" (gene profiles)
                     - max_results: Maximum number of results (default: 100)
                     - filters: Dict with organism, entry_type, platform, date_from, date_to
                     - fetch_details: Whether to fetch full summaries (default: True)

        Returns:
            AdapterResult containing dataset information
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
            database = kwargs.get("database", "gds")
            max_results = kwargs.get("max_results", self.config.get("max_results", 100))
            filters = kwargs.get("filters")
        else:
            query = input_data["query"]
            database = input_data.get("database", kwargs.get("database", "gds"))
            max_results = input_data.get("max_results", kwargs.get("max_results", 100))
            filters = input_data.get("filters", kwargs.get("filters"))

        fetch_details = kwargs.get("fetch_details", True)

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.34)
        await asyncio.sleep(rate_delay)

        # Check if query is a GEO accession
        acc_type = self._is_geo_accession(query)
        if acc_type:
            logger.info(f"GEO: Detected accession query: {query} (type: {acc_type})")
            dataset = await self._fetch_by_accession_async(query)

            if not dataset:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Dataset not found for accession: {query}",
                    metadata={"source": "geo", "query": query}
                )

            return AdapterResult(
                success=True,
                data={
                    "query": query,
                    "query_type": "accession",
                    "accession_type": acc_type,
                    "dataset": dataset
                },
                cache_hit=False,
                metadata={
                    "source": "geo",
                    "query": query,
                    "adapter_version": self.version,
                    "accession": query
                }
            )

        # Keyword/gene search
        search_result = await self._search_geo_async(query, database, max_results, filters)

        if search_result is None:
            return AdapterResult(
                success=False,
                data=None,
                error="Failed to search GEO",
                metadata={"source": "geo", "query": query}
            )

        uids = search_result.get("uids", [])
        total_count = search_result.get("count", 0)

        if not uids:
            return AdapterResult(
                success=True,
                data={
                    "uids": [],
                    "datasets": [],
                    "summary": self._summarize_results([]),
                    "query": query,
                    "total_count": 0
                },
                metadata={"source": "geo", "query": query, "total_results": 0}
            )

        # Fetch dataset details if requested
        datasets = []
        if fetch_details:
            datasets = await self._fetch_summaries_async(uids, database) or []

        # Create summary
        summary = self._summarize_results(datasets)

        result_data = {
            "uids": uids,
            "datasets": datasets,
            "summary": summary,
            "query": query,
            "database": database,
            "total_count": total_count,
            "returned_count": len(datasets)
        }

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "geo",
                "query": query,
                "adapter_version": self.version,
                "database": database,
                "total_results": total_count,
                "datasets_fetched": len(datasets)
            }
        )
