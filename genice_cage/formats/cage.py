#!/usr/bin/env python
# coding: utf-8
"""
A GenIce format plugin to detect cage-like topologies.

Usage: 
    % genice CS1 -r 2 2 2 -f cage[12,14-16:maxring=6] 
    % genice CRN1 -f cage[3-10:json] 
    % genice CRN1 -f cage[3-10:yaplot] 
    % genice CS2 -w tip4p -f cage[gromacs:-16:maxring=6]
    % analice traj.gro -O OW -H HW[12] -w tip4p -f cage[quad]
    % analice traj.gro -O OW -H HW[12] -w tip4p -f cage[quad:json]

It may not work with a small structure. (In the example above, the unit cell of CS1 is extended to 2x2x2 so as to avoid detecting cell-spanning wierd cages.)

Options:
    Cage sizes to be listed, separated by commas and ranged with hyphens. (e.g. -4,6,8-10,16-) (default is 3-8)
    maxring=x  Specify the maximum ring size (default=8).
    json       Output values in [JSON](https://www.json.org/) format.
    yaplot     Visualize cages with [Yaplot](https://github.com/vitroid/Yaplot/). Cages are drawn in different layers according to the number of faces, and faces are colored according to the number of vertices.
    gromacs    Output individual cages in Gromacs format. (EXPERIMENTAL)
    quad       Quadcage order parameter to identify the Frank-Kasper-type crystal structures.[JMM2011] Cages sizes and maximum ring size are set appropriately automatically.

* [JMM2011] Jacobson, L. C., Matsumoto, M. & Molinero, V. Order parameters for the multistep crystallization of clathrate hydrates. J. Chem. Phys. 135, 074501 (2011).[doi:10.1063/1.3613667](https://doi.org/10.1063/1.3613667)
"""

# standard modules
import json
from logging import getLogger
from collections import defaultdict
from math import log2

# external modules
import networkx as nx
import numpy as np
from attrdict import AttrDict

# public modules developed by myself
from countrings import countrings_nx as cr
from genice_cage.polyhed import Polyhed
import yaplotlib as yp

def centerOfMass(members, rpos):
    logger = getLogger()
    dsum = np.zeros(3)
    for member in members:
        d = rpos[member] - rpos[members[0]]
        d -= np.floor(d+0.5)
        dsum += d
    com = rpos[members[0]] + dsum / len(members)
    com -= np.floor(com)
    return com

def is_spanning(ring, coord):
    logger = getLogger()
    dsum = 0
    for i in range(len(ring)):
        j,k = ring[i-1], ring[i]
        d = coord[j] - coord[k]
        d -= np.floor(d+0.5)
        dsum += d
    return np.linalg.norm(dsum) > 1e-4

def hook2(lattice):
    global options

    logger = getLogger()
    logger.info("Hook2: Cages and vitrites")

    cell = lattice.repcell.mat
    positions = lattice.reppositions
    graph = nx.Graph(lattice.graph) #undirected
    maxringsize = options.maxring
    ringlist = [[int(x) for x in ring] for ring in cr.CountRings(graph).rings_iter(maxringsize)]
    for ring in ringlist:
        assert not is_spanning(ring, positions), "Some ring is spanning the cell."
    ringpos = [centerOfMass(ringnodes, positions) for ringnodes in ringlist]
    logger.info("  Rings: {0}".format(len(ringlist)))
    maxcagesize = max(options.sizes)
    cages = []
    for cage in Polyhed(ringlist, maxcagesize):
        if len(cage) in options.sizes:
            cages.append(list(cage))
    logger.info("  Cages: {0}".format(len(cages)))
    cagepos = np.array([centerOfMass(cage, ringpos) for cage in cages])
    if options.gromacs:
        options.rings = ringlist
        options.cages = cages
        logger.debug("##2 {0}".format(cages))
        logger.info("Hook2: end.")
        return

    if options.quad:
        oncage = defaultdict(list)
        for cage in cages:
            nodes = set()
            for ringid in cage:
                nodes |= set(ringlist[ringid])
            for node in nodes:
                oncage[node].append(len(cage))
        op = dict()
        for node in oncage:
            count = [0 for i in range(17)]
            for v in oncage[node]:
                count[v] += 1
            v = "{0}{1}{2}{3}".format(count[12], count[14], count[15], count[16])
            op[node] = v
            
        stat = defaultdict(int)
        for node in sorted(op):
            v = op[node]
            stat[v] += 1

        if options.json:
            output = dict()
            N = positions.shape[0]
            output["op"] = {str(k):v for k,v in op.items()}
            output["stat"] = {k:v/N for k,v in stat.items()}
            print(json.dumps(output, indent=2, sort_keys=True))
        else:
            for node in sorted(op):
                print(node, op[node])
            print("# Statistics")
            for v in sorted(stat):
                print("{0} {1} {2}/{3}".format(v, stat[v]/positions.shape[0], stat[v], positions.shape[0]))
        #ideal = {"CS2": {"2002": 0.7058823529411765,
        #                 "3001": 0.23529411764705882,
        #                 "4000": 0.058823529411764705},
        #         "CS1": {"0400": 0.13043478260869565,
        #                 "1300": 0.8695652173913043},
        #}
        #for ref in ideal:
        #    dKL = 0
        #    for v in sorted(stat):
        #        if v in ideal[ref]:
        #            dKL += ideal[ref][v]*(log2(ideal[ref][v]) - log2(stat[v]/positions.shape[0]))
        #    print("{0} dKL={1}".format(ref, dKL))
    elif options.json:
        output = dict()
        output["rings"] = ringlist
        output["cages"] = cages
        output["ringpos"] = [[x,y,z] for x,y,z in ringpos]
        output["cagepos"] = [[x,y,z] for x,y,z in cagepos]
        print(json.dumps(output, indent=2, sort_keys=True))
    elif options.yaplot:
        s = ""
        for c, cage in enumerate(cages):
            nodes = dict()
            cagesize = len(cage)
            for ringid in cage:
                ns = ringlist[ringid]
                for node in ns:
                    if node not in nodes:
                        # relative pos of the node
                        nodepos = positions[node] - cagepos[c]
                        nodepos -= np.floor( nodepos + 0.5 )
                        # shrink a little
                        nodes[node] = nodepos * 0.9
                s += yp.Color(len(ns))
                s += yp.Layer(cagesize)
                polygon = (np.array([nodes[node] for node in ns]) + cagepos[c]) @ cell
                s += yp.Polygon(polygon)
        print(s + "\n")
    else:
        # human-friendly redundant format
        for cageid, cage in enumerate(cages):
            print("Cage {0}: ({1}, {2}, {3}) {4} hedron".format(cageid, *cagepos[cageid], len(cage)))
            for ringid in sorted(cage):
                print("  Ring {0}: ({1}, {2}, {3}) {4} gon".format(ringid, *ringpos[ringid], len(ringlist[ringid])))
                print("    Nodes: {0}".format(ringlist[ringid]))
    logger.info("Hook2: end.")
    return True # terminate


def hook6(lattice):
    global options
    lattice.logger.info("Hook6: Output in Gromacs format.")
    lattice.logger.info("  Total number of atoms: {0}".format(len(lattice.atoms)))
    cellmat = lattice.repcell.mat
    s = ""
    mols = defaultdict(list)
    for atom in lattice.atoms:
        resno, resname, atomname, position, order = atom
        logger.debug(atom)
        mols[order].append(atom)
    for cage in options.cages:
        logger.debug(cage)
        s += "Generated by GenIce https://github.com/vitroid/GenIce \n"
        f = ""
        cagemols = set()
        atomcount = 0
        for ring in cage:
            cagemols |= set(options.rings[ring])
        for mol in cagemols:
            for atom in mols[mol]:
                resno, resname, atomname, position, order = atom
                f += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(order,resname, atomname, atomcount+1,position[0],position[1],position[2])
                atomcount += 1
        s += "{0}\n".format(atomcount) + f
        if cellmat[1,0] == 0 and cellmat[2,0] == 0 and cellmat[2,1] == 0:
            s += "    {0} {1} {2}\n".format(cellmat[0,0],cellmat[1,1],cellmat[2,2])
        else:
            assert cellmat[0,1] == 0 and cellmat[0,2] == 0 and cellmat[1,2] == 0
            s += "    {0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(cellmat[0,0],
                                                                    cellmat[1,1],
                                                                    cellmat[2,2],
                                                                    cellmat[0,1],
                                                                    cellmat[0,2],
                                                                    cellmat[1,0],
                                                                    cellmat[1,2],
                                                                    cellmat[2,0],
                                                                    cellmat[2,1],
            )
    print(s,end="")
    lattice.logger.info("Hook6: end.")
    


def argparser(lattice, arg):
    global options

    logger = getLogger()
    logger.info("Hook0: Parse options for cage plugin.")

    options=AttrDict({"sizes":set(),
                      "maxring":8,
                      "json":False,
                      "gromacs":False,
                      "yaplot":False,
                      "quad":False,
    })

    if arg != "":
        for a in arg.split(":"):
            decl = a.split("=")
            if len(decl) == 2:
                if decl[0] == "maxring":
                    options.maxring = int(decl[1])
                else:
                    assert False, "Wrong declaration."
            elif a in ("json", "JSON"):
                options.json = True
            elif a in ("gromacs",):
                options.gromacs = True
            elif a in ("yaplot",):
                options.yaplot = True
            elif a in ("quad",):
                options.quad = True
                options.sizes = set([12,14,15,16])
                options.maxring = 6
            else:
                # value list for cage sizes
                for v in a.split(","):
                    w = v.split("-")
                    if len(w) == 2:
                        if w[0] == "":
                            w[0] = "1"
                        if w[1] == "":
                            w[1] = "20"
                        for x in range(int(w[0]), int(w[1])+1):
                            options.sizes.add(x)
                    else:
                        options.sizes.add(int(v))

    if len(options.sizes) == 0:
        options.sizes = {3,4,5,6,7,8}

    logger.info("  Max ring size: {0}".format(options.maxring))
    logger.info("  Cage sizes:    {0}".format(options.sizes))

    logger.info("Hook0: end.")

hooks = {0:argparser, 2:hook2, 6:hook6}
