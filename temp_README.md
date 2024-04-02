![Logo]({{tool.genice.urls.logo}})

# [{{project.name}}]({{project.urls.Homepage}})

A format plugin for [GenIce2]({{tool.genice.urls.repository}}) to detect cages.

version {{version}}

## Requirements

{% for item in tool.poetry.dependencies %}* {{item}}{{tool.poetry.dependencies[item]}}
{% endfor %}

## Installation from PyPI

```shell
% pip install {{project.name}}
```

## Usage

{%- filter indent %}
    {{usage}}
{%- endfilter %}


## Test in place

```shell
% make test
```

## Algorithms

* M. Matsumoto, A. Baba, and I. Ohmine, Topological building blocks of hydrogen bond network in water, J. Chem. Phys. 127, 134504 (2007); [doi:10.1063/1.2772627](http://dx.doi.org/doi:10.1063/1.2772627)
