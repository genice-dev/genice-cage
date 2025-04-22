#!/usr/bin/env python
from jinja2 import Template, BaseLoader, Environment, FileSystemLoader
import toml
import genice2_cage.formats.cage

import sys

project = toml.load("pyproject.toml")

project |= {
    "usage": genice2_cage.formats.cage.desc["usage"],
    "version": project["tool"]["poetry"]["version"],
}

t = Environment(loader=FileSystemLoader(searchpath=".")).get_template(sys.argv[1])
markdown_en = t.render(**project)
print(markdown_en)
