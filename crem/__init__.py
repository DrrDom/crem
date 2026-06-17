#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 14-08-2019
# copyright       : Pavel Polishchuk 2019
# license         :
#==============================================================================


def _resolve_version():
    try:
        from pathlib import Path
        from setuptools_scm import get_version
        return get_version(root=Path(__file__).resolve().parent.parent)
    except Exception:
        pass
    try:
        from ._version import version as _v
        return _v
    except ImportError:
        return "0.0.0+unknown"


__version__ = _resolve_version()
del _resolve_version
