
from fastapi import APIRouter, Response
import json

router = APIRouter()

@router.get('/health')
def health(db_ok: bool = True, redis_ok: bool = True):
    checks = {'db': db_ok, 'redis': redis_ok}
    healthy = all(checks.values())
    code = 200 if healthy else 503
    return Response(
        content=json.dumps({'status':'ok' if healthy else 'degraded', **checks}),
        media_type='application/json',
        status_code=code
    )
