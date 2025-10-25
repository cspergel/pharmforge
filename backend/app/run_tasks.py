
import asyncio
def run_async_in_background(coro_fn, *args, **kwargs):
    def _runner():
        asyncio.run(coro_fn(*args, **kwargs))
    return _runner
