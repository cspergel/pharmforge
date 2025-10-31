"""
Background task utilities for PharmForge

Provides wrapper functions to execute async coroutines in background tasks
compatible with FastAPI's BackgroundTasks.
"""
import asyncio


def run_async_in_background(coro_fn, *args, **kwargs):
    """
    Wrap an async coroutine function for execution as a background task.

    FastAPI's BackgroundTasks expects synchronous functions, but many of our
    pipeline operations are async. This wrapper allows async functions to be
    executed properly as background tasks.

    Args:
        coro_fn: Async function (coroutine function) to execute
        *args: Positional arguments to pass to the coroutine function
        **kwargs: Keyword arguments to pass to the coroutine function

    Returns:
        Synchronous wrapper function that can be added to BackgroundTasks

    Usage:
        >>> async def async_pipeline(smiles: str, run_id: str):
        ...     await process_pipeline(smiles, run_id)
        ...
        >>> # In FastAPI endpoint:
        >>> background_tasks.add_task(
        ...     run_async_in_background(async_pipeline, "CCO", "run_123")
        ... )

    Example:
        >>> from fastapi import BackgroundTasks
        >>> async def my_async_task(value: int):
        ...     await asyncio.sleep(0.1)
        ...     print(f"Processed: {value}")
        ...
        >>> background_tasks = BackgroundTasks()
        >>> background_tasks.add_task(
        ...     run_async_in_background(my_async_task, 42)
        ... )
    """
    def _runner():
        asyncio.run(coro_fn(*args, **kwargs))
    return _runner
