"""
DisGeNET Adapter - Fetches gene-disease and variant-disease associations from DisGeNET API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class DisGeNETAdapter(AdapterProtocol):
    """
    Adapter for DisGeNET REST API
    Fetches gene-disease associations and variant-disease associations

    Note: DisGeNET API requires authentication for most endpoints.
    This adapter uses the public endpoints which have limited functionality.
    For full access, obtain an API key from https://www.disgenet.org/signup
    """

    BASE_URL = "https://www.disgenet.org/api"

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize DisGeNET adapter

        Args:
            api_key: Optional API key for authenticated requests
        """
        super().__init__(
            name="disgenet",
            adapter_type="api",
            config={
                "rate_limit_delay": 1.0,  # 1 request/second for public API
                "timeout": 60,
                "max_results": 100
            }
        )
        self.version = "1.0.0"
        self.api_key = api_key
        self._session = None

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid gene symbol, disease name, or variant ID

        Args:
            input_data: String containing gene symbol, disease name, or variant ID

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        return True

    def _get_headers(self) -> Dict[str, str]:
        """
        Get HTTP headers including API key if available

        Returns:
            Dictionary of HTTP headers
        """
        headers = {
            "Accept": "application/json",
            "Content-Type": "application/json"
        }

        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"

        return headers

    async def _get_gene_disease_associations(
        self,
        gene_symbol: str,
        max_results: int = 100
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get disease associations for a gene

        Args:
            gene_symbol: Gene symbol (e.g., "BRCA2")
            max_results: Maximum number of results

        Returns:
            List of disease associations or None
        """
        # DisGeNET public API endpoint for gene-disease associations
        url = f"{self.BASE_URL}/gda/gene/{gene_symbol}"

        params = {
            "source": "ALL",  # All data sources
            "limit": max_results
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(
                    url,
                    headers=self._get_headers(),
                    params=params,
                    timeout=timeout
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        # DisGeNET returns data in different formats depending on version
                        if isinstance(data, list):
                            return data
                        elif isinstance(data, dict):
                            return data.get("results", [])
                        return []
                    elif response.status == 401:
                        logger.error("DisGeNET: Authentication required. Please provide API key.")
                        return None
                    elif response.status == 404:
                        logger.warning(f"DisGeNET: Gene not found: {gene_symbol}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"DisGeNET API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"DisGeNET: Timeout fetching data for {gene_symbol}")
            return None
        except Exception as e:
            logger.error(f"DisGeNET: Error fetching data for {gene_symbol}: {e}")
            return None

    async def _get_disease_gene_associations(
        self,
        disease_name: str,
        max_results: int = 100
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get gene associations for a disease

        Args:
            disease_name: Disease name (e.g., "breast cancer")
            max_results: Maximum number of results

        Returns:
            List of gene associations or None
        """
        # Note: This endpoint may require authentication
        url = f"{self.BASE_URL}/gda/disease/{disease_name}"

        params = {
            "source": "ALL",
            "limit": max_results
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(
                    url,
                    headers=self._get_headers(),
                    params=params,
                    timeout=timeout
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        if isinstance(data, list):
                            return data
                        elif isinstance(data, dict):
                            return data.get("results", [])
                        return []
                    elif response.status == 401:
                        logger.error("DisGeNET: Authentication required for disease queries")
                        return None
                    elif response.status == 404:
                        logger.warning(f"DisGeNET: Disease not found: {disease_name}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"DisGeNET API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"DisGeNET: Timeout fetching data for {disease_name}")
            return None
        except Exception as e:
            logger.error(f"DisGeNET: Error fetching data for {disease_name}: {e}")
            return None

    async def _get_variant_disease_associations(
        self,
        variant_id: str,
        max_results: int = 100
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get disease associations for a variant

        Args:
            variant_id: Variant ID (e.g., "rs123456")
            max_results: Maximum number of results

        Returns:
            List of disease associations or None
        """
        url = f"{self.BASE_URL}/vda/variant/{variant_id}"

        params = {
            "source": "ALL",
            "limit": max_results
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(
                    url,
                    headers=self._get_headers(),
                    params=params,
                    timeout=timeout
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        if isinstance(data, list):
                            return data
                        elif isinstance(data, dict):
                            return data.get("results", [])
                        return []
                    elif response.status == 401:
                        logger.error("DisGeNET: Authentication required for variant queries")
                        return None
                    elif response.status == 404:
                        logger.warning(f"DisGeNET: Variant not found: {variant_id}")
                        return []
                    else:
                        error_text = await response.text()
                        logger.error(f"DisGeNET API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"DisGeNET: Timeout fetching data for {variant_id}")
            return None
        except Exception as e:
            logger.error(f"DisGeNET: Error fetching data for {variant_id}: {e}")
            return None

    def _summarize_associations(self, associations: List[Dict[str, Any]], assoc_type: str) -> Dict[str, Any]:
        """
        Summarize association data into useful metrics

        Args:
            associations: List of association records
            assoc_type: Type of association ("gda" or "vda")

        Returns:
            Dictionary of summarized metrics
        """
        if not associations:
            return {
                "num_associations": 0,
                "sources": [],
                "evidence_types": [],
                "avg_score": None,
                "max_score": None
            }

        sources = set()
        evidence_types = set()
        scores = []

        for assoc in associations:
            # Collect sources
            source = assoc.get("source")
            if source:
                sources.add(source)

            # Collect evidence types
            evidence = assoc.get("evidence_type") or assoc.get("evidenceType")
            if evidence:
                evidence_types.add(evidence)

            # Collect scores
            score = assoc.get("score") or assoc.get("geneDisease_Score") or assoc.get("variantDisease_Score")
            if score:
                try:
                    scores.append(float(score))
                except (ValueError, TypeError):
                    pass

        return {
            "num_associations": len(associations),
            "sources": list(sources),
            "evidence_types": list(evidence_types),
            "avg_score": sum(scores) / len(scores) if scores else None,
            "max_score": max(scores) if scores else None,
            "min_score": min(scores) if scores else None
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the DisGeNET query

        Args:
            input_data: Gene symbol, disease name, or variant ID
            **kwargs: Additional parameters
                - query_type: "gene", "disease", or "variant" (default: "gene")
                - max_results: Maximum number of results (default: 100)

        Returns:
            AdapterResult containing DisGeNET association data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a non-empty string"
            )

        query_type = kwargs.get("query_type", "gene")
        max_results = kwargs.get("max_results", self.config.get("max_results", 100))

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 1.0)
        await asyncio.sleep(rate_delay)

        # Execute query based on type
        associations = None

        try:
            if query_type == "gene":
                associations = await self._get_gene_disease_associations(input_data, max_results)
            elif query_type == "disease":
                associations = await self._get_disease_gene_associations(input_data, max_results)
            elif query_type == "variant":
                associations = await self._get_variant_disease_associations(input_data, max_results)
            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Invalid query_type: {query_type}. Must be 'gene', 'disease', or 'variant'"
                )

            if associations is None:
                return AdapterResult(
                    success=False,
                    data=None,
                    error="Failed to fetch data from DisGeNET (may require API key)",
                    metadata={
                        "source": "disgenet",
                        "query": input_data,
                        "query_type": query_type
                    }
                )

            # Summarize results
            summary = self._summarize_associations(associations, query_type)

            # Format the detailed associations
            formatted_associations = []
            for assoc in associations[:max_results]:
                formatted = {
                    "gene_symbol": assoc.get("gene_symbol") or assoc.get("geneSymbol"),
                    "disease_name": assoc.get("disease_name") or assoc.get("diseaseName"),
                    "disease_id": assoc.get("disease_id") or assoc.get("diseaseId"),
                    "score": assoc.get("score") or assoc.get("geneDisease_Score") or assoc.get("variantDisease_Score"),
                    "evidence_type": assoc.get("evidence_type") or assoc.get("evidenceType"),
                    "source": assoc.get("source"),
                    "pmid": assoc.get("pmid"),  # PubMed IDs
                }

                # Add variant-specific fields if present
                if query_type == "variant":
                    formatted["variant_id"] = assoc.get("variant_id") or assoc.get("variantId")
                    formatted["chromosome"] = assoc.get("chromosome")
                    formatted["position"] = assoc.get("position")

                formatted_associations.append(formatted)

            result_data = {
                "query": input_data,
                "query_type": query_type,
                "summary": summary,
                "associations": formatted_associations
            }

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "disgenet",
                    "query": input_data,
                    "query_type": query_type,
                    "adapter_version": self.version,
                    "authenticated": self.api_key is not None
                }
            )

        except Exception as e:
            logger.error(f"DisGeNET: Error processing results: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error processing DisGeNET data: {str(e)}",
                metadata={
                    "source": "disgenet",
                    "query": input_data,
                    "query_type": query_type
                }
            )
