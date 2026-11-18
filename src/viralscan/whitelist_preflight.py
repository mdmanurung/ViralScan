"""Compatibility re-export: the barcode-whitelist preflight lives in :mod:`viralscan.chemistry`.

Kept because ``tests/_frozen_oracles`` (frozen, not edited) imports this path. New code
imports from :mod:`viralscan.chemistry` (PLAN WHITELIST-FOLD-01).
"""

from viralscan.chemistry import (
    DEFAULT_MIN_MATCH_RATE,
    DEFAULT_SAMPLE,
    WhitelistCheck,
    check_whitelist,
    load_whitelist,
    whitelist_match_rate,
)

__all__ = [
    "DEFAULT_MIN_MATCH_RATE",
    "DEFAULT_SAMPLE",
    "WhitelistCheck",
    "check_whitelist",
    "load_whitelist",
    "whitelist_match_rate",
]
