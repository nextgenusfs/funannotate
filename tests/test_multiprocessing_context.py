"""The aux_scripts call runMultiProgress() from unguarded top-level code.

macOS starts workers with "spawn", which re-imports and re-runs the calling
script in every worker. runMultiProgress must ask for the "fork" context.
"""
import multiprocessing

import funannotate.library as lib


def _square(x):
    return x * x


def test_run_multi_progress_uses_the_fork_context(monkeypatch):
    requested = []
    real = multiprocessing.get_context

    def recording(method=None):
        requested.append(method)
        return real(method)

    monkeypatch.setattr(multiprocessing, "get_context", recording)
    lib.runMultiProgress(_square, [1, 2, 3], 2, progress=False)
    assert requested == ["fork"]


def test_run_multi_progress_still_runs_the_jobs():
    # Not a result check of the helper's return value; it ends without raising
    # and without leaving workers behind.
    lib.runMultiProgress(_square, list(range(6)), 2, progress=False)
    assert multiprocessing.active_children() == []
