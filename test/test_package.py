from __future__ import annotations

import importlib.metadata

import archicrop as m


def test_version():
    assert importlib.metadata.version("archicrop") == m.__version__
