__all__ = ["SiteSelector", "SiteSelectorManager", "PolarBinRenderer"]

def __getattr__(name: str):
    if name in __all__:
        from .matplotlib import selectors, renderers
        mapping = {
            "SiteSelector": selectors.SiteSelector,
            "SiteSelectorManager": selectors.SiteSelectorManager,
            "PolarBinRenderer": renderers.PolarBinRenderer,
        }
        return mapping[name]
    raise AttributeError(name)
