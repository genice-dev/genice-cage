# [genice-cage](https://github.com/vitroid/genice-cage/)

A [GenIce](https://github.com/vitroid/GenIce) plugin to detect cage-like topologies.

version 0.1

## Requirements

* countrings>=0.1.6
* genice>=0.25

## Installation from PyPI

    % pip install genice_cage

## Manual Installation

### System-wide installation

    % make install

### Private installation

Copy the files in genice_cage/formats/ into your local formats/ folder.

## Usage

    
    Usage: 
        % genice CS1 -r 2 2 2 -f cage[12,14-16:maxring=6] 
        % genice CRN1 -f cage[3-10:json] 
        % genice CRN1 -f cage[3-10:yaplot] 
    
    It may not work with a small structure. (In the example above, the unit cell of CS1 is extended to 2x2x2 so as to avoid detecting cell-spanning wierd cages.)
    
    Options:
        Cage sizes to be listed, separated by commas and ranged with hyphens. (default is 3 to 8)
        maxring=x  Specify the maximum ring size (default=8).
        json       Output in [JSON](https://www.json.org/) format.
        yaplot     Visualize with [Yaplot](https://github.com/vitroid/Yaplot/). Cages are drawn in different layers according to the number of faces, and faces are colored according to the number of vertices.

## Test in place

    % make test

## Algorithms

* M. Matsumoto, A. Baba, and I. Ohmine, Topological building blocks of hydrogen bond network in water, J. Chem. Phys. 127, 134504 (2007); [doi:10.1063/1.2772627](http://dx.doi.org/doi:10.1063/1.2772627)
