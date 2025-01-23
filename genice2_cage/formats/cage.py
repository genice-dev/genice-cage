#!/usr/bin/env python
# coding: utf-8

"""
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
"""

desc = {
    "ref": {
        "JMM2011": "Jacobson, L. C., Matsumoto, M. & Molinero, V. Order parameters for the multistep crystallization of clathrate hydrates. J. Chem. Phys. 135, 074501 (2011).[doi:10.1063/1.3613667](https://doi.org/10.1063/1.3613667)"
    },
    "brief": "Cage analysis.",
    "usage": __doc__,
}


# standard modules
import json
import string
from logging import getLogger
from collections import defaultdict, Counter

# from math import log2

# external modules
import networkx as nx
import numpy as np

# old and not python3.11 compat
# from attrdict import AttrDict

# public modules developed by myself
from cycless.polyhed import polyhedra_iter, cage_to_graph
from cycless.cycles import centerOfMass, cycles_iter
import genice2.formats
import yaplotlib as yp
from genice2.molecules import serialize
from graphstat import GraphStat


def assign_unused_label(basename, labels):
    enum = 0
    label = f"A{basename}"
    while label in labels:
        char = string.ascii_lowercase[enum]
        label = f"A{basename}{char}"
        enum += 1
    return label


def make_cage_expression(ring_ids, ringlist):
    # list the sizes of rings in the ringlist
    ringsizes = np.array([len(ring) for ring in ringlist])
    # count the occurrences
    values, counts = np.unique(ringsizes, return_counts=True)
    # sort and stringify
    index = " ".join([f"{ringsize}^{counts[i]}" for i, ringsize in enumerate(values)])
    return index


def rangeparser(s, min=1, max=20):
    # value list for cage sizes
    values = set()
    for v in s.split(","):
        w = v.split("-")
        if len(w) == 2:
            if w[0] == "":
                w[0] = min
            else:
                w[0] = int(w[0])
            if w[1] == "":
                w[1] = max
            else:
                w[1] = int(w[1])
            for x in range(w[0], w[1] + 1):
                values.add(x)
        else:
            values.add(int(v))
    return values


class Format(genice2.formats.Format):
    def __init__(self, **kwargs):
        logger = getLogger()
        options = {
            "sizes": set(),
            "ring": None,
            "json": False,
            "json2": False,
            "gromacs": False,
            "yaplot": False,
            "quad": False,
            "python": False,
        }
        unknown = dict()
        for k, v in kwargs.items():
            if k == "maxring":
                options["ring"] = [x for x in range(3, int(v) + 1)]
            elif k == "ring":
                options["ring"] = rangeparser(v, min=3, max=8)
            elif k == "sizes":
                options["sizes"] = rangeparser(v, min=3, max=20)
            elif k in ("json", "JSON"):
                options["json"] = v
            elif k in ("gromacs",):
                options["gromacs"] = v
            elif k in ("yaplot",):
                options["yaplot"] = v
            elif k in ("python",):
                options["python"] = v
            elif k in ("quad",):
                options["quad"] = v
                options["sizes"] = set([12, 14, 15, 16])
                options["ring"] = set([5, 6])
            elif k in ("json2",):
                options["json2"] = v
            else:
                # value list for cage sizes
                options["sizes"] = rangeparser(k, min=3)
        super().__init__(**kwargs)

        if len(options["sizes"]) == 0:
            options["sizes"] = set([x for x in range(3, 17)])
        if options["ring"] is None:
            options["ring"] = set([x for x in range(3, 8)])

        logger.info("  Ring sizes: {0}".format(options["ring"]))
        logger.info("  Cage sizes: {0}".format(options["sizes"]))
        self.options = options

    def hooks(self):
        return {2: self.Hook2, 6: self.Hook6}

    def Hook2(self, ice):
        logger = getLogger()
        logger.info("Hook2: Cages and vitrites")

        db = GraphStat()
        labels = set()
        g_id2label = dict()
        cagetypes = []

        cell = ice.repcell.mat
        positions = ice.reppositions
        graph = nx.Graph(ice.graph)  # undirected
        ringsize = self.options["ring"]
        ringlist = [
            [int(x) for x in ring]
            for ring in cycles_iter(graph, max(ringsize), pos=positions)
        ]
        ringpos = np.array(
            [centerOfMass(ringnodes, positions) for ringnodes in ringlist]
        )
        logger.info("  Rings: {0}".format(len(ringlist)))
        maxcagesize = max(self.options["sizes"])
        cages = []
        for cage in polyhedra_iter(ringlist, maxcagesize):
            if len(cage) in self.options["sizes"]:
                valid = True
                for ringid in cage:
                    if len(ringlist[ringid]) not in ringsize:
                        valid = False
                if valid:
                    cages.append(list(cage))
        logger.info("  Cages: {0}".format(len(cages)))
        cagepos = np.array([centerOfMass(cage, ringpos) for cage in cages])

        if self.options["gromacs"]:
            self.options["rings"] = ringlist
            self.options["cages"] = cages
            logger.debug("##2 {0}".format(cages))
            logger.info("Hook2: end.")
            return

        if self.options["quad"]:
            oncage = defaultdict(list)
            for cage in cages:
                nodes = set()
                for ringid in cage:
                    nodes |= set(ringlist[ringid])
                for node in nodes:
                    oncage[node].append(len(cage))
            op = dict()
            for node in oncage:
                count = [0 for _ in range(17)]
                for v in oncage[node]:
                    count[v] += 1
                v = "{0}{1}{2}{3}".format(count[12], count[14], count[15], count[16])
                op[node] = v

            stat = defaultdict(int)
            for node in sorted(op):
                v = op[node]
                stat[v] += 1

            if self.options["json"]:
                output = dict()
                N = positions.shape[0]
                output["op"] = {str(k): v for k, v in op.items()}
                output["stat"] = {k: v / N for k, v in stat.items()}
                print(json.dumps(output, indent=2, sort_keys=True))
            else:
                for node in sorted(op):
                    print(node, op[node])
                print("# Statistics")
                for v in sorted(stat):
                    print(
                        "{0} {1} {2}/{3}".format(
                            v, stat[v] / positions.shape[0], stat[v], positions.shape[0]
                        )
                    )

            # ideal = {"CS2": {"2002": 0.7058823529411765,
            #                  "3001": 0.23529411764705882,
            #                  "4000": 0.058823529411764705},
            #          "CS1": {"0400": 0.13043478260869565,
            #                  "1300": 0.8695652173913043},
            # }
            # for ref in ideal:
            #    dKL = 0
            #    for v in sorted(stat):
            #        if v in ideal[ref]:
            #            dKL += ideal[ref][v]*(log2(ideal[ref][v]) - log2(stat[v]/positions.shape[0]))
            #    print("{0} dKL={1}".format(ref, dKL))

        elif self.options["json"]:
            output = dict()
            output["rings"] = ringlist
            output["cages"] = cages
            output["ringpos"] = [[x, y, z] for x, y, z in ringpos]
            output["cagepos"] = [[x, y, z] for x, y, z in cagepos]
            print(json.dumps(output, indent=2, sort_keys=True))
        elif self.options["json2"]:
            for cage in cages:
                g = cage_to_graph(cage, ringlist)
                cagesize = len(cage)
                g_id = db.query_id(g)
                # if it is a new cage type
                if g_id < 0:
                    # new type!
                    # register the last query
                    g_id = db.register()
                    # prepare a new label
                    label = assign_unused_label(cagesize, labels)
                    g_id2label[g_id] = label
                    labels.add(label)
                    # cage expression
                    index = make_cage_expression(cage, ringlist)
                    logger.info(f"  Cage type: {label} ({index})")
                else:
                    label = g_id2label[g_id]
                cagetypes.append(label)
            if len(cagepos) == 0:
                logger.info("  No cages detected.")
                print("No cages detected.")
            output = dict(
                Rings=len(ringlist), Cages=len(cages), details=dict(Counter(cagetypes))
            )
            print(json.dumps(output, indent=2, sort_keys=True))
        elif self.options["yaplot"]:
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
                            nodepos -= np.floor(nodepos + 0.5)
                            # shrink a little
                            nodes[node] = nodepos * 0.9
                    s += yp.Color(len(ns))
                    s += yp.Layer(cagesize)
                    polygon = (
                        np.array([nodes[node] for node in ns]) + cagepos[c]
                    ) @ cell
                    s += yp.Polygon(polygon)
            print(s)
        elif self.options["python"]:
            import graphstat as gs

            db = gs.GraphStat()
            labels = set()
            g_id2label = dict()
            print('cages="""')
            for c, cage in enumerate(cages):
                g = cage_to_graph(cage, ringlist)
                cagesize = len(cage)
                g_id = db.query_id(g)
                if g_id < 0:
                    g_id = db.register()
                    enum = 0
                    label = "{0}".format(cagesize, enum)
                    while label in labels:
                        enum += 1
                        label = "{0}_{1}".format(cagesize, enum)
                    g_id2label[g_id] = label
                    labels.add(label)
                else:
                    label = g_id2label[g_id]
                print("{0:10s} {1:.4f} {2:.4f} {3:.4f}".format(label, *cagepos[c]))
            print('"""')
        else:
            # human-friendly redundant format
            for cageid, cage in enumerate(cages):
                print(
                    "Cage {0}: ({1}, {2}, {3}) {4} hedron".format(
                        cageid, *cagepos[cageid], len(cage)
                    )
                )
                for ringid in sorted(cage):
                    print(
                        "  Ring {0}: ({1}, {2}, {3}) {4} gon".format(
                            ringid, *ringpos[ringid], len(ringlist[ringid])
                        )
                    )
                    print("    Nodes: {0}".format(ringlist[ringid]))
        logger.info("Hook2: end.")
        return True  # terminate

    def Hook6(self, ice):
        logger = getLogger()
        logger.info("Hook6: Output in Gromacs format.")
        atoms = []
        for mols in ice.universe:
            atoms += serialize(mols)
        cellmat = ice.repcell.mat
        s = ""
        mols = defaultdict(list)
        for atom in atoms:
            resno, resname, atomname, position, order = atom
            logger.debug(atom)
            mols[order].append(atom)
        for cage in self.options["cages"]:
            logger.debug(cage)
            s += "Generated by GenIce https://github.com/vitroid/GenIce \n"
            f = ""
            cagemols = set()
            atomcount = 0
            for ring in cage:
                cagemols |= set(self.options["rings"][ring])
            for mol in cagemols:
                for atom in mols[mol]:
                    resno, resname, atomname, position, order = atom
                    f += "{0:5d}{1:5s}{2:>5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}\n".format(
                        order,
                        resname,
                        atomname,
                        atomcount + 1,
                        position[0],
                        position[1],
                        position[2],
                    )
                    atomcount += 1
            s += "{0}\n".format(atomcount) + f
            if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
                s += "    {0} {1} {2}\n".format(
                    cellmat[0, 0], cellmat[1, 1], cellmat[2, 2]
                )
            else:
                assert cellmat[0, 1] == 0 and cellmat[0, 2] == 0 and cellmat[1, 2] == 0
                s += "    {0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(
                    cellmat[0, 0],
                    cellmat[1, 1],
                    cellmat[2, 2],
                    cellmat[0, 1],
                    cellmat[0, 2],
                    cellmat[1, 0],
                    cellmat[1, 2],
                    cellmat[2, 0],
                    cellmat[2, 1],
                )
        print(s, end="")
        logger.info("Hook6: end.")
