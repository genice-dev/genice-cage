#!/usr/bin/env python
import sys
from genice2_cage.formats.cage import __doc__ as doc
import distutils.core
import jinja2 as jj

genice_link="[GenIce2](https://github.com/vitroid/GenIce)"

setup = distutils.core.run_setup("setup.py")

d = {
    "usage"   : "\n".join(doc.splitlines()[2:]),
    "summary" : doc.splitlines()[1].replace("GenIce2", genice_link),
    "version" : setup.get_version(),
    "package" : setup.get_name(),
    "base"    : setup.get_name().replace("-", "_"),
    "url"     : setup.get_url(),
    "genice"  : genice_link,
    "requires": setup.install_requires,
}


print(jj.Template(sys.stdin.read()).render(d))
