def __getattr__(name):
    if name == "__version__":
        from .__version__ import __version__
        globals()["__version__"] = __version__
        return __version__
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))
