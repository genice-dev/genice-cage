![Logo](https://raw.githubusercontent.com/vitroid/GenIce/develop/logo/genice-v0.png)

# [genice2-cage](https://github.com/vitroid/genice-cage)

A format plugin for [GenIce2](https://github.com/vitroid/GenIce) to detect cages.

version 2.2

## Requirements

* python^3.9
* GenIce2^2.2.5.2
* graphstat^0.2.1
* cycless^0.4.2
* yaplotlib^0.1.2
* numpy^1.26.4
* networkx^3.2.1


## Installation from PyPI

```shell
% pip install genice2-cage
```

## Usage
        
    A GenIce2 format plugin to detect cage-like topologies.

    Usage:
        % genice2 CS1 -r 2 2 2 -f cage[12,14-16:ring=-6]
        % genice2 CRN1 -f cage[sizes=3-10:json]
        % genice2 CRN1 -f cage[sizes=3-10:yaplot]
        % genice2 CS2 -w tip4p -f cage[gromacs:sizes=-16:ring=5,6]
        % analice2 traj.gro -O OW -H HW[12] -w tip4p -f cage[quad]
        % analice2 traj.gro -O OW -H HW[12] -w tip4p -f cage[quad:json]
        % genice2 FAU -r 2 2 2 -f cage[-26:maxring=12:json2]

    It may not work with a small structure. (In the example above, the unit cell of CS1 is extended to 2x2x2 so as to avoid detecting cell-spanning wierd cages.)

    Options:
        Cage sizes to be listed, separated by commas and ranged with hyphens. (e.g. -4,6,8-10,16-) (default is 3-16)
        ring=3,5-6 Specify the ring sizes that cages are built of (default is 3-8, maximum is 8).
        json       Output values in [JSON](https://www.json.org/) format.
        json2      Output values in [JSON](https://www.json.org/) format (Assess cage locations based on HB network topology by labeling them).
        yaplot     Visualize cages with [Yaplot](https://github.com/vitroid/Yaplot/). Cages are drawn in different layers according to the number of faces, and faces are colored according to the number of vertices.
        gromacs    Output individual cages in Gromacs format. (EXPERIMENTAL)
        quad       Quadcage order parameter to identify the Frank-Kasper-type crystal structures.[JMM2011] Cages sizes and maximum ring size are set appropriately automatically.
        python     Output cage types in python format convenient for GenIce lattice modules.
    * [JMM2011] Jacobson, L. C., Matsumoto, M. & Molinero, V. Order parameters for the multistep crystallization of clathrate hydrates. J. Chem. Phys. 135, 074501 (2011).[doi:10.1063/1.3613667](https://doi.org/10.1063/1.3613667)



## Test in place

```shell
% make test
```

## Algorithms

* M. Matsumoto, A. Baba, and I. Ohmine, Topological building blocks of hydrogen bond network in water, J. Chem. Phys. 127, 134504 (2007); [doi:10.1063/1.2772627](http://dx.doi.org/doi:10.1063/1.2772627)
