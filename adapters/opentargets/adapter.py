"""
OpenTargets Adapter - Fetches target-disease associations from OpenTargets Platform GraphQL API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging
from aiohttp import TCPConnector, ClientTimeout

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class OpenTargetsAdapter(AdapterProtocol):
    """
    Adapter for OpenTargets Platform GraphQL API
    Fetches target-disease associations, drug mechanisms, and genetic evidence
    """

    # Correct endpoint - the old platform-api.opentargets.io is deprecated
    BASE_URL = "https://api.platform.opentargets.org/api/v4/graphql"
    FALLBACK_URLS = [
        "https://api.platform.opentargets.org/api/v4/graphql",
        "http://api.platform.opentargets.org/api/v4/graphql"  # HTTP fallback if HTTPS fails
    ]

    def __init__(self, config: Optional[Dict[str, Any]] = None):
        """
        Initialize OpenTargets adapter

        Args:
            config: Optional configuration dictionary with keys:
                - rate_limit_delay: Delay between requests in seconds (default: 0.1)
                - timeout: Request timeout in seconds (default: 60)
                - max_results: Maximum number of results to return (default: 50)
                - base_url: Custom API endpoint URL (default: uses BASE_URL)
                - max_retries: Maximum number of retry attempts (default: 3)
                - retry_delay: Initial retry delay in seconds (default: 1.0)
        """
        default_config = {
            "rate_limit_delay": 0.1,  # 10 requests/second
            "timeout": 60,
            "max_results": 50,
            "max_retries": 3,
            "retry_delay": 1.0
        }
        merged_config = {**default_config, **(config or {})}

        super().__init__(
            name="opentargets",
            adapter_type="api",
            config=merged_config
        )
        self.version = "1.2.0"

        # Use custom URL if provided, otherwise use default
        self.api_url = self.config.get("base_url", self.BASE_URL)

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid gene ID, disease ID, or drug name

        Args:
            input_data: String containing gene ID, disease ID, or drug name

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        return True

    async def _query_graphql(self, query: str, variables: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Execute a GraphQL query against OpenTargets API with retry logic

        Args:
            query: GraphQL query string
            variables: Variables for the query

        Returns:
            Query result data or None on error
        """
        timeout = ClientTimeout(
            total=self.config.get("timeout", 60),
            connect=30,  # Connection timeout
            sock_read=30  # Socket read timeout
        )

        payload = {
            "query": query,
            "variables": variables
        }

        max_retries = self.config.get("max_retries", 3)
        retry_delay = self.config.get("retry_delay", 1.0)

        # Try multiple URLs if the primary fails
        urls_to_try = [self.api_url] + [url for url in self.FALLBACK_URLS if url != self.api_url]

        last_error = None

        for url_index, url in enumerate(urls_to_try):
            for attempt in range(max_retries):
                try:
                    # Create a new session with proper DNS and connection settings
                    connector = TCPConnector(
                        family=0,  # Use both IPv4 and IPv6
                        ssl=False if url.startswith("http://") else None,  # Disable SSL verification for HTTP
                        force_close=True,  # Force close connections to avoid stale connections
                        enable_cleanup_closed=True
                    )

                    async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
                        logger.debug(f"OpenTargets: Attempting query to {url} (attempt {attempt + 1}/{max_retries})")

                        async with session.post(
                            url,
                            json=payload,
                            headers={"Content-Type": "application/json"}
                        ) as response:
                            if response.status == 200:
                                data = await response.json()
                                if "errors" in data:
                                    error_messages = [err.get("message", str(err)) for err in data["errors"]]
                                    logger.error(f"OpenTargets GraphQL errors: {error_messages}")
                                    return None
                                logger.info(f"OpenTargets: Successfully queried {url}")
                                return data.get("data")
                            else:
                                error_text = await response.text()
                                last_error = f"HTTP {response.status}: {error_text}"
                                logger.warning(f"OpenTargets API error {response.status} on {url}: {error_text}")
                                # Don't retry on client errors (4xx)
                                if 400 <= response.status < 500:
                                    return None
                                # Retry on server errors (5xx)
                                if attempt < max_retries - 1:
                                    await asyncio.sleep(retry_delay * (2 ** attempt))  # Exponential backoff
                                    continue

                except asyncio.TimeoutError:
                    last_error = f"Timeout after {self.config.get('timeout')}s"
                    logger.warning(f"OpenTargets: Timeout on {url} (attempt {attempt + 1}/{max_retries})")
                    if attempt < max_retries - 1:
                        await asyncio.sleep(retry_delay * (2 ** attempt))
                        continue

                except aiohttp.ClientConnectorError as e:
                    last_error = f"Connection error: {str(e)}"
                    logger.warning(f"OpenTargets: Connection error on {url}: {e}")
                    if attempt < max_retries - 1:
                        await asyncio.sleep(retry_delay * (2 ** attempt))
                        continue
                    # If this is the last attempt with this URL, try the next URL
                    break

                except aiohttp.ClientError as e:
                    last_error = f"Client error: {str(e)}"
                    logger.warning(f"OpenTargets: Client error on {url}: {e}")
                    if attempt < max_retries - 1:
                        await asyncio.sleep(retry_delay * (2 ** attempt))
                        continue

                except Exception as e:
                    last_error = f"Unexpected error: {str(e)}"
                    logger.warning(f"OpenTargets: Unexpected error on {url}: {e}")
                    if attempt < max_retries - 1:
                        await asyncio.sleep(retry_delay * (2 ** attempt))
                        continue

            # If we exhausted retries for this URL, the loop will continue to the next URL

        # If we get here, all URLs and retries failed
        logger.error(
            f"OpenTargets: All connection attempts failed. "
            f"Tried URLs: {urls_to_try}. Last error: {last_error}. "
            f"Please check your network connection and DNS settings."
        )
        return None

    async def _get_target_info(self, gene_id: str) -> Optional[Dict[str, Any]]:
        """
        Get information about a target gene

        Args:
            gene_id: Ensembl gene ID (e.g., "ENSG00000139618" for BRCA2)

        Returns:
            Target information or None
        """
        query = """
        query TargetInfo($ensemblId: String!) {
          target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            approvedName
            biotype
            functionDescriptions
            pathways {
              pathway
              pathwayId
            }
            tractability {
              label
              modality
              value
            }
            associatedDiseases(page: {index: 0, size: 10}) {
              count
              rows {
                disease {
                  id
                  name
                }
                score
                datatypeScores {
                  id
                  score
                }
              }
            }
          }
        }
        """

        variables = {"ensemblId": gene_id}
        result = await self._query_graphql(query, variables)

        if result and "target" in result:
            return result["target"]
        return None

    async def _get_disease_associations(self, gene_id: str, max_results: int = 50) -> Optional[List[Dict[str, Any]]]:
        """
        Get disease associations for a target gene

        Args:
            gene_id: Ensembl gene ID
            max_results: Maximum number of results to return

        Returns:
            List of disease associations or None
        """
        query = """
        query DiseaseAssociations($ensemblId: String!, $size: Int!) {
          target(ensemblId: $ensemblId) {
            associatedDiseases(page: {index: 0, size: $size}) {
              rows {
                disease {
                  id
                  name
                  description
                  therapeuticAreas {
                    id
                    name
                  }
                }
                score
                datatypeScores {
                  id
                  score
                }
              }
            }
          }
        }
        """

        variables = {"ensemblId": gene_id, "size": max_results}
        result = await self._query_graphql(query, variables)

        if result and "target" in result and result["target"]:
            associations = result["target"].get("associatedDiseases", {}).get("rows", [])
            return associations
        return None

    async def _get_drug_mechanisms(self, gene_id: str, max_results: int = 50) -> Optional[List[Dict[str, Any]]]:
        """
        Get drug mechanisms for a target gene

        Args:
            gene_id: Ensembl gene ID
            max_results: Maximum number of results to return

        Returns:
            List of drug mechanisms or None
        """
        query = """
        query DrugMechanisms($ensemblId: String!, $size: Int!) {
          target(ensemblId: $ensemblId) {
            knownDrugs(size: $size) {
              rows {
                drug {
                  id
                  name
                  drugType
                  maximumClinicalTrialPhase
                }
                mechanismOfAction
                phase
                status
                urls {
                  url
                  name
                }
              }
            }
          }
        }
        """

        variables = {"ensemblId": gene_id, "size": max_results}
        result = await self._query_graphql(query, variables)

        if result and "target" in result and result["target"]:
            drugs = result["target"].get("knownDrugs", {}).get("rows", [])
            return drugs
        return None

    async def _search_target(self, query_string: str) -> Optional[str]:
        """
        Search for a target by gene symbol or name

        Args:
            query_string: Gene symbol or name to search for

        Returns:
            Ensembl gene ID or None
        """
        query = """
        query SearchTarget($queryString: String!) {
          search(queryString: $queryString, entityNames: ["target"], page: {index: 0, size: 1}) {
            hits {
              id
              name
              entity
            }
          }
        }
        """

        variables = {"queryString": query_string}
        result = await self._query_graphql(query, variables)

        if result and "search" in result:
            hits = result["search"].get("hits", [])
            if hits:
                return hits[0].get("id")
        return None

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the OpenTargets query

        Args:
            input_data: Gene ID (ENSG...), gene symbol (e.g., "BRCA2"), or disease ID
            **kwargs: Additional parameters
                - query_type: "target_info", "disease_associations", "drug_mechanisms" (default: "all")
                - max_results: Maximum number of results (default: 50)

        Returns:
            AdapterResult containing OpenTargets data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a non-empty string"
            )

        query_type = kwargs.get("query_type", "all")
        max_results = kwargs.get("max_results", self.config.get("max_results", 50))

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.1)
        await asyncio.sleep(rate_delay)

        gene_id = input_data

        # If input doesn't look like an Ensembl ID, try to search for it
        if not gene_id.startswith("ENSG"):
            logger.info(f"OpenTargets: Searching for gene: {gene_id}")
            found_id = await self._search_target(gene_id)
            if not found_id:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Gene not found: {gene_id}",
                    metadata={"source": "opentargets", "query": input_data}
                )
            gene_id = found_id
            logger.info(f"OpenTargets: Found gene ID: {gene_id}")

        # Execute requested queries
        result_data = {
            "gene_id": gene_id,
            "query_type": query_type
        }

        try:
            if query_type in ["all", "target_info"]:
                target_info = await self._get_target_info(gene_id)
                if target_info:
                    result_data["target_info"] = {
                        "symbol": target_info.get("approvedSymbol"),
                        "name": target_info.get("approvedName"),
                        "biotype": target_info.get("biotype"),
                        "function": target_info.get("functionDescriptions", [])[:3],  # Top 3 functions
                        "go_terms": [
                            go.get("term", {}).get("name")
                            for go in target_info.get("go", [])[:10]  # Top 10 GO terms
                        ]
                    }

            if query_type in ["all", "disease_associations"]:
                associations = await self._get_disease_associations(gene_id, max_results)
                if associations is not None:
                    result_data["disease_associations"] = [
                        {
                            "disease_id": assoc.get("disease", {}).get("id"),
                            "disease_name": assoc.get("disease", {}).get("name"),
                            "overall_score": assoc.get("score"),
                            "datatype_scores": {
                                dt.get("id"): dt.get("score")
                                for dt in assoc.get("datatypeScores", [])
                            },
                            "therapeutic_areas": [
                                ta.get("name")
                                for ta in assoc.get("disease", {}).get("therapeuticAreas", [])
                            ]
                        }
                        for assoc in associations[:max_results]
                    ]

            if query_type in ["all", "drug_mechanisms"]:
                drugs = await self._get_drug_mechanisms(gene_id, max_results)
                if drugs is not None:
                    result_data["drug_mechanisms"] = [
                        {
                            "drug_id": drug.get("drug", {}).get("id"),
                            "drug_name": drug.get("drug", {}).get("name"),
                            "drug_type": drug.get("drug", {}).get("drugType"),
                            "mechanism": drug.get("mechanismOfAction"),
                            "phase": drug.get("phase"),
                            "max_clinical_phase": drug.get("drug", {}).get("maximumClinicalTrialPhase"),
                            "status": drug.get("status")
                        }
                        for drug in drugs[:max_results]
                    ]

        except Exception as e:
            logger.error(f"OpenTargets: Error processing results: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error processing OpenTargets data: {str(e)}",
                metadata={"source": "opentargets", "gene_id": gene_id}
            )

        # Check if we got any data - note that empty lists are still valid data
        has_data = any(key in result_data for key in ["target_info", "disease_associations", "drug_mechanisms"])

        if not has_data:
            return AdapterResult(
                success=False,
                data=None,
                error=f"No data found for gene: {gene_id}. The gene may not exist or the API may be experiencing issues.",
                metadata={"source": "opentargets", "gene_id": gene_id}
            )

        return AdapterResult(
            success=True,
            data=result_data,
            cache_hit=False,
            metadata={
                "source": "opentargets",
                "gene_id": gene_id,
                "original_query": input_data,
                "adapter_version": self.version
            }
        )
