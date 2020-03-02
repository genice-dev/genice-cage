#!/usr/bin/env python
import sys
from genice_cage.formats.cage import __doc__ as doc
import distutils.core
import jinja2 as jj

setup = distutils.core.run_setup("setup.py")

d = {
    "usage"   : "\n".join(doc.splitlines()[2:]),
    "version" : setup.get_version(),
    "package" : setup.get_name(),
    "url"     : setup.get_url(),
    "genice"  : "[GenIce](https://github.com/vitroid/GenIce)",
    "requires": setup.install_requires,
}


print(jj.Template(sys.stdin.read()).render(d))
