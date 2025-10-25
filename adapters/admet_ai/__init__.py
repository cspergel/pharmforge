"""
ADMET-AI Adapter

Provides ADMET property predictions using the ADMET-AI library
(Chemprop-based models trained on TDC datasets).

Reference: https://github.com/swansonk14/admet_ai
"""

from .adapter import ADMETaiAdapter

__all__ = ["ADMETaiAdapter"]
