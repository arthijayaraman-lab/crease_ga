class Plugins:
    """mBuild recipe object."""

    pass


plugins = {}
from pkg_resources import iter_entry_points

import sys
#if sys.version_info < (3, 8):
#        from importlib_metadata import entry_points
#else:
#        from importlib.metadata import entry_points

#plugins = entry_points(group="crease_ga.plugins")


for ep in iter_entry_points(group="crease_ga.plugins",name=None):
    plugins[ep.name] = ep
