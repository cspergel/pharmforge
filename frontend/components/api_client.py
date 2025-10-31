"""
API Client for PharmForge Frontend
Handles all backend communication with proper error handling and caching
"""
import os
import requests
from typing import Optional, Dict, Any, List
import logging
from dataclasses import dataclass
from datetime import datetime

logger = logging.getLogger(__name__)


@dataclass
class APIResponse:
    """Standardized API response wrapper"""
    success: bool
    data: Optional[Dict[str, Any]] = None
    error: Optional[str] = None
    status_code: Optional[int] = None


class PharmForgeAPIClient:
    """
    Client for communicating with PharmForge backend API

    Reads API URL from environment variable PHARMFORGE_API_URL
    Default: http://backend:8000
    """

    def __init__(self, base_url: Optional[str] = None, timeout: int = 10):
        """
        Initialize API client

        Args:
            base_url: Optional override for API URL
            timeout: Default timeout for requests in seconds
        """
        self.base_url = base_url or os.getenv("PHARMFORGE_API_URL", "http://localhost:8000")
        self.timeout = timeout
        logger.info(f"PharmForge API Client initialized with base_url: {self.base_url}")

    def _make_request(
        self,
        method: str,
        endpoint: str,
        data: Optional[Dict] = None,
        params: Optional[Dict] = None,
        timeout: Optional[int] = None
    ) -> APIResponse:
        """
        Make HTTP request with error handling

        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint path
            data: Request body (for POST/PUT)
            params: Query parameters
            timeout: Request timeout override

        Returns:
            APIResponse with success status and data/error
        """
        url = f"{self.base_url}{endpoint}"
        timeout = timeout or self.timeout

        try:
            logger.debug(f"{method} {url} with data={data}, params={params}")

            response = requests.request(
                method=method,
                url=url,
                json=data,
                params=params,
                timeout=timeout
            )

            # Check for HTTP errors
            if response.status_code >= 400:
                error_msg = f"HTTP {response.status_code}: {response.text}"
                logger.error(error_msg)
                return APIResponse(
                    success=False,
                    error=error_msg,
                    status_code=response.status_code
                )

            # Parse JSON response
            try:
                response_data = response.json()
            except ValueError:
                response_data = {"raw": response.text}

            return APIResponse(
                success=True,
                data=response_data,
                status_code=response.status_code
            )

        except requests.exceptions.Timeout:
            error_msg = f"Request timeout after {timeout}s"
            logger.error(error_msg)
            return APIResponse(success=False, error=error_msg)

        except requests.exceptions.ConnectionError as e:
            error_msg = f"Connection error: {str(e)}"
            logger.error(error_msg)
            return APIResponse(success=False, error=error_msg)

        except Exception as e:
            error_msg = f"Unexpected error: {str(e)}"
            logger.error(error_msg, exc_info=True)
            return APIResponse(success=False, error=error_msg)

    # ==================== Health ====================

    def get_health(self) -> APIResponse:
        """
        Check system health

        Returns:
            APIResponse with health status
        """
        return self._make_request("GET", "/health", timeout=2)

    # ==================== Adapters ====================

    def list_adapters(self) -> APIResponse:
        """
        List all registered adapters

        Returns:
            APIResponse with adapter list
        """
        return self._make_request("GET", "/api/adapters/list")

    def get_adapter(self, adapter_name: str) -> APIResponse:
        """
        Get details for a specific adapter

        Args:
            adapter_name: Name of the adapter

        Returns:
            APIResponse with adapter metadata
        """
        return self._make_request("GET", f"/api/adapters/{adapter_name}/info")

    def test_adapter(self, adapter_name: str, smiles: str, protein_data: Optional[Dict] = None) -> APIResponse:
        """
        Test an adapter with a SMILES string

        Args:
            adapter_name: Name of the adapter to test
            smiles: SMILES string to test with
            protein_data: Optional protein data for docking adapters

        Returns:
            APIResponse with test results
        """
        data = {"input_data": {"smiles": smiles}, "use_cache": True}

        if protein_data:
            data["input_data"]["protein_data"] = protein_data

        # Use longer timeout for computationally intensive adapters
        # DiffDock and other ML models can take several minutes
        timeout = 300  # 5 minutes for docking/ML adapters
        if adapter_name in ['diffdock', 'gnina', 'vina', 'openmm', 'aizynthfinder']:
            timeout = 300
        else:
            timeout = 60  # 1 minute for other adapters

        return self._make_request(
            "POST",
            f"/api/adapters/{adapter_name}/execute",
            data=data,
            timeout=timeout
        )

    # ==================== Runs ====================

    def create_run(
        self,
        name: str,
        input_type: str,
        nl_query: Optional[str] = None,
        smiles_list: Optional[List[str]] = None,
        pipeline_config: Optional[Dict] = None,
        target_protein: Optional[str] = None,
        n_candidates: int = 10
    ) -> APIResponse:
        """
        Create a new pipeline run

        Args:
            name: Human-readable run name
            input_type: Type of input ("nl", "batch", "evolution")
            nl_query: Natural language query (for NL mode)
            smiles_list: List of SMILES strings (for batch mode)
            pipeline_config: Optional pipeline configuration
            target_protein: Target protein identifier
            n_candidates: Number of top candidates to return

        Returns:
            APIResponse with run details
        """
        # Build request payload
        payload = {
            "name": name,
            "input_type": input_type,
        }

        if nl_query:
            payload["nl_query"] = nl_query

        # For now, use a default SMILES list for NL queries
        # In production, this will be generated by Arcana
        if not smiles_list:
            smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]  # Default test molecules

        payload["input_smiles"] = smiles_list

        if pipeline_config:
            payload["pipeline_config"] = pipeline_config

        if target_protein:
            payload["target_protein"] = target_protein

        payload["n_candidates"] = n_candidates

        return self._make_request("POST", "/api/v1/runs", data=payload, timeout=30)

    def get_run(self, run_id: str) -> APIResponse:
        """
        Get status and details for a pipeline run

        Args:
            run_id: Unique run identifier

        Returns:
            APIResponse with run details
        """
        return self._make_request("GET", f"/api/v1/runs/{run_id}")

    def list_runs(
        self,
        page: int = 1,
        page_size: int = 20,
        status: Optional[str] = None,
        user_id: Optional[str] = None
    ) -> APIResponse:
        """
        List all pipeline runs with pagination

        Args:
            page: Page number (1-indexed)
            page_size: Number of items per page
            status: Optional status filter
            user_id: Optional user filter

        Returns:
            APIResponse with paginated runs
        """
        params = {
            "page": page,
            "page_size": page_size
        }

        if status:
            params["status"] = status
        if user_id:
            params["user_id"] = user_id

        return self._make_request("GET", "/api/v1/runs", params=params)

    def get_run_results(self, run_id: str, top_n: int = 10) -> APIResponse:
        """
        Get ranked results for a completed run

        Args:
            run_id: Unique run identifier
            top_n: Number of top results to return

        Returns:
            APIResponse with ranked compounds
        """
        return self._make_request(
            "GET",
            f"/api/v1/runs/{run_id}/results",
            params={"top_n": top_n}
        )

    def delete_run(self, run_id: str) -> APIResponse:
        """
        Delete a pipeline run and all associated results

        Args:
            run_id: Unique run identifier

        Returns:
            APIResponse with deletion status
        """
        return self._make_request("DELETE", f"/api/v1/runs/{run_id}")


# Global singleton instance
_api_client = None


def get_api_client() -> PharmForgeAPIClient:
    """
    Get or create global API client instance

    Returns:
        PharmForgeAPIClient singleton
    """
    global _api_client
    if _api_client is None:
        _api_client = PharmForgeAPIClient()
    return _api_client
