"""
GTEx Adapter - Fetches gene expression data from GTEx Portal API
"""
from typing import Any, Dict, Optional, List
import aiohttp
import asyncio
import logging

from backend.core.adapters.protocol import AdapterProtocol, AdapterResult

logger = logging.getLogger(__name__)


class GTExAdapter(AdapterProtocol):
    """
    Adapter for GTEx (Genotype-Tissue Expression) Portal API
    Fetches gene expression levels by tissue and eQTL data

    Note: GTEx API has rate limits and may require special handling
    """

    BASE_URL = "https://gtexportal.org/api/v2"

    def __init__(self):
        super().__init__(
            name="gtex",
            adapter_type="api",
            config={
                "rate_limit_delay": 0.2,  # 5 requests/second
                "timeout": 60,
                "max_tissues": 100
            }
        )
        self.version = "1.0.0"

    def validate_input(self, input_data: Any) -> bool:
        """
        Validate that input is a valid gene symbol or ID

        Args:
            input_data: String containing gene symbol or Ensembl ID

        Returns:
            True if valid, False otherwise
        """
        if not isinstance(input_data, str):
            return False
        if len(input_data) == 0:
            return False
        return True

    async def _get_gene_expression(
        self,
        gene_symbol: str,
        dataset: str = "gtex_v8"
    ) -> Optional[Dict[str, Any]]:
        """
        Get gene expression data across tissues

        Args:
            gene_symbol: Gene symbol (e.g., "BRCA2")
            dataset: GTEx dataset version

        Returns:
            Dictionary of expression data or None
        """
        url = f"{self.BASE_URL}/expression/medianGeneExpression"

        params = {
            "gencodeId": gene_symbol,
            "datasetId": dataset,
            "format": "json"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    elif response.status == 404:
                        logger.warning(f"GTEx: Gene not found: {gene_symbol}")
                        return None
                    else:
                        error_text = await response.text()
                        logger.error(f"GTEx API error {response.status}: {error_text}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"GTEx: Timeout fetching expression for {gene_symbol}")
            return None
        except Exception as e:
            logger.error(f"GTEx: Error fetching expression for {gene_symbol}: {e}")
            return None

    async def _get_gene_info(self, gene_symbol: str) -> Optional[Dict[str, Any]]:
        """
        Get gene information from GTEx

        Args:
            gene_symbol: Gene symbol

        Returns:
            Gene information or None
        """
        url = f"{self.BASE_URL}/reference/gene"

        params = {
            "gencodeId": gene_symbol,
            "format": "json"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    else:
                        logger.warning(f"GTEx: Gene info not found for {gene_symbol}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"GTEx: Timeout fetching gene info for {gene_symbol}")
            return None
        except Exception as e:
            logger.error(f"GTEx: Error fetching gene info for {gene_symbol}: {e}")
            return None

    async def _get_tissue_expression(
        self,
        gene_symbol: str,
        tissue: str,
        dataset: str = "gtex_v8"
    ) -> Optional[Dict[str, Any]]:
        """
        Get gene expression in a specific tissue

        Args:
            gene_symbol: Gene symbol
            tissue: Tissue name
            dataset: GTEx dataset version

        Returns:
            Expression data for tissue or None
        """
        url = f"{self.BASE_URL}/expression/geneExpression"

        params = {
            "gencodeId": gene_symbol,
            "tissueSiteDetailId": tissue,
            "datasetId": dataset,
            "format": "json"
        }

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data
                    else:
                        logger.warning(f"GTEx: Expression not found for {gene_symbol} in {tissue}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"GTEx: Timeout fetching tissue expression")
            return None
        except Exception as e:
            logger.error(f"GTEx: Error fetching tissue expression: {e}")
            return None

    async def _get_eqtls(
        self,
        gene_symbol: str,
        tissue: Optional[str] = None,
        dataset: str = "gtex_v8"
    ) -> Optional[List[Dict[str, Any]]]:
        """
        Get eQTL (expression quantitative trait loci) data

        Args:
            gene_symbol: Gene symbol
            tissue: Optional tissue name
            dataset: GTEx dataset version

        Returns:
            List of eQTL data or None
        """
        url = f"{self.BASE_URL}/association/singleTissueEqtl"

        params = {
            "gencodeId": gene_symbol,
            "datasetId": dataset,
            "format": "json"
        }

        if tissue:
            params["tissueSiteDetailId"] = tissue

        timeout = aiohttp.ClientTimeout(total=self.config.get("timeout", 60))

        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(url, params=params, timeout=timeout) as response:
                    if response.status == 200:
                        data = await response.json()
                        if isinstance(data, dict):
                            return data.get("singleTissueEqtl", [])
                        elif isinstance(data, list):
                            return data
                        return []
                    else:
                        logger.warning(f"GTEx: eQTLs not found for {gene_symbol}")
                        return None
        except asyncio.TimeoutError:
            logger.error(f"GTEx: Timeout fetching eQTLs")
            return None
        except Exception as e:
            logger.error(f"GTEx: Error fetching eQTLs: {e}")
            return None

    def _format_expression_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Format expression data into standardized structure

        Args:
            data: Raw GTEx expression data

        Returns:
            Formatted expression data
        """
        if not data:
            return {
                "tissues": [],
                "median_expression": {},
                "max_expression_tissue": None,
                "max_expression_value": None,
                "num_tissues_expressed": 0
            }

        # GTEx returns data in different formats depending on endpoint
        # Handle both possible formats
        tissues = []
        expression_values = {}

        if "medianGeneExpression" in data:
            for entry in data["medianGeneExpression"]:
                tissue = entry.get("tissueSiteDetailId", "")
                median = entry.get("median", 0)
                tissues.append(tissue)
                expression_values[tissue] = median

        elif "geneExpression" in data:
            for entry in data["geneExpression"]:
                tissue = entry.get("tissueSiteDetailId", "")
                expression = entry.get("expression", 0)
                tissues.append(tissue)
                expression_values[tissue] = expression

        # Find tissue with highest expression
        max_tissue = None
        max_value = None
        if expression_values:
            max_tissue = max(expression_values, key=expression_values.get)
            max_value = expression_values[max_tissue]

        # Count tissues where gene is expressed (threshold: TPM > 1)
        num_expressed = sum(1 for val in expression_values.values() if val > 1)

        return {
            "tissues": tissues,
            "median_expression": expression_values,
            "max_expression_tissue": max_tissue,
            "max_expression_value": max_value,
            "num_tissues_expressed": num_expressed
        }

    async def execute(self, input_data: Any, **kwargs) -> AdapterResult:
        """
        Execute the GTEx query

        Args:
            input_data: Gene symbol or Ensembl gene ID
            **kwargs: Additional parameters
                - query_type: "expression", "tissue", or "eqtl" (default: "expression")
                - tissue: Specific tissue name (for "tissue" or "eqtl" queries)
                - dataset: GTEx dataset version (default: "gtex_v8")

        Returns:
            AdapterResult containing GTEx expression data
        """
        # Validate input
        if not self.validate_input(input_data):
            return AdapterResult(
                success=False,
                data=None,
                error="Invalid input: must be a non-empty string"
            )

        query_type = kwargs.get("query_type", "expression")
        tissue = kwargs.get("tissue")
        dataset = kwargs.get("dataset", "gtex_v8")

        # Rate limiting
        rate_delay = self.config.get("rate_limit_delay", 0.2)
        await asyncio.sleep(rate_delay)

        gene_symbol = input_data

        try:
            result_data = {
                "gene_symbol": gene_symbol,
                "query_type": query_type,
                "dataset": dataset
            }

            if query_type == "expression":
                # Get median expression across all tissues
                expression_data = await self._get_gene_expression(gene_symbol, dataset)

                if not expression_data:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"No expression data found for gene: {gene_symbol}",
                        metadata={"source": "gtex", "gene": gene_symbol}
                    )

                formatted = self._format_expression_data(expression_data)
                result_data.update(formatted)

                # Also try to get gene info
                gene_info = await self._get_gene_info(gene_symbol)
                if gene_info:
                    result_data["gene_info"] = gene_info

            elif query_type == "tissue":
                if not tissue:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error="tissue parameter required for 'tissue' query type"
                    )

                # Get expression in specific tissue
                tissue_data = await self._get_tissue_expression(gene_symbol, tissue, dataset)

                if not tissue_data:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"No expression data found for {gene_symbol} in {tissue}",
                        metadata={"source": "gtex", "gene": gene_symbol, "tissue": tissue}
                    )

                result_data["tissue"] = tissue
                result_data["expression_data"] = tissue_data

            elif query_type == "eqtl":
                # Get eQTL data
                eqtl_data = await self._get_eqtls(gene_symbol, tissue, dataset)

                if eqtl_data is None:
                    return AdapterResult(
                        success=False,
                        data=None,
                        error=f"No eQTL data found for gene: {gene_symbol}",
                        metadata={"source": "gtex", "gene": gene_symbol}
                    )

                result_data["eqtls"] = eqtl_data
                result_data["num_eqtls"] = len(eqtl_data)

                if tissue:
                    result_data["tissue"] = tissue

            else:
                return AdapterResult(
                    success=False,
                    data=None,
                    error=f"Invalid query_type: {query_type}. Must be 'expression', 'tissue', or 'eqtl'"
                )

            return AdapterResult(
                success=True,
                data=result_data,
                cache_hit=False,
                metadata={
                    "source": "gtex",
                    "gene": gene_symbol,
                    "query_type": query_type,
                    "dataset": dataset,
                    "adapter_version": self.version
                }
            )

        except Exception as e:
            logger.error(f"GTEx: Error processing results: {e}")
            return AdapterResult(
                success=False,
                data=None,
                error=f"Error processing GTEx data: {str(e)}",
                metadata={"source": "gtex", "gene": gene_symbol}
            )
