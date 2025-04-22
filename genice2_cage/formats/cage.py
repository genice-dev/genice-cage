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

# 標準ライブラリ
import json
import string
from logging import getLogger
from collections import defaultdict, Counter

# 外部ライブラリ
import networkx as nx
import numpy as np

# カスタムモジュール
from cycless.polyhed import polyhedra_iter, cage_to_graph
from cycless.cycles import centerOfMass, cycles_iter
import genice2.formats
import yaplotlib as yp
from genice2.molecules import serialize
from graphstat import GraphStat

# グローバル変数
desc = {
    "ref": {
        "JMM2011": "Jacobson, L. C., Matsumoto, M. & Molinero, V. Order parameters for the multistep crystallization of clathrate hydrates. J. Chem. Phys. 135, 074501 (2011).[doi:10.1063/1.3613667](https://doi.org/10.1063/1.3613667)"
    },
    "brief": "Cage analysis.",
    "usage": __doc__,
}


class CageAnalyzer:
    """ケージ分析のためのユーティリティクラス"""

    @staticmethod
    def assign_unused_label(basename, labels):
        """未使用のラベルを生成する"""
        enum = 0
        label = f"A{basename}"
        while label in labels:
            char = string.ascii_lowercase[enum]
            label = f"A{basename}{char}"
            enum += 1
        return label

    @staticmethod
    def make_cage_expression(ring_ids, ringlist):
        """ケージの表現を生成する"""
        ringsizes = np.array([len(ring) for ring in ringlist])
        values, counts = np.unique(ringsizes, return_counts=True)
        return " ".join(
            [f"{ringsize}^{counts[i]}" for i, ringsize in enumerate(values)]
        )

    @staticmethod
    def parse_range(s, min=1, max=20):
        """範囲指定をパースする"""
        values = set()
        for v in s.split(","):
            w = v.split("-")
            if len(w) == 2:
                start = min if w[0] == "" else int(w[0])
                end = max if w[1] == "" else int(w[1])
                values.update(range(start, end + 1))
            else:
                values.add(int(v))
        return values


class Format(genice2.formats.Format):
    """GenIce2のフォーマットプラグイン"""

    def __init__(self, **kwargs):
        """初期化"""
        logger = getLogger()
        self.options = self._initialize_options(kwargs)
        super().__init__(**kwargs)
        self._validate_options()
        logger.info(f"  Ring sizes: {self.options['ring']}")
        logger.info(f"  Cage sizes: {self.options['sizes']}")

    def _initialize_options(self, kwargs):
        """オプションの初期化"""
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

        for k, v in kwargs.items():
            if k == "maxring":
                options["ring"] = list(range(3, int(v) + 1))
            elif k == "ring":
                options["ring"] = CageAnalyzer.parse_range(v, min=3, max=8)
            elif k == "sizes":
                options["sizes"] = CageAnalyzer.parse_range(v, min=3, max=20)
            elif k in ("json", "JSON"):
                options["json"] = v
            elif k == "gromacs":
                options["gromacs"] = v
            elif k == "yaplot":
                options["yaplot"] = v
            elif k == "python":
                options["python"] = v
            elif k == "quad":
                options["quad"] = v
                options["sizes"] = {12, 14, 15, 16}
                options["ring"] = {5, 6}
            elif k == "json2":
                options["json2"] = v
            else:
                options["sizes"] = CageAnalyzer.parse_range(k, min=3)

        return options

    def _validate_options(self):
        """オプションの検証"""
        if not self.options["sizes"]:
            self.options["sizes"] = set(range(3, 17))
        if self.options["ring"] is None:
            self.options["ring"] = set(range(3, 8))

    def hooks(self):
        """フックの定義"""
        return {2: self.Hook2, 6: self.Hook6}

    def Hook2(self, ice):
        """ケージとビトリットの分析を行うフック"""
        logger = getLogger()
        logger.info("Hook2: Cages and vitrites analysis started")

        try:
            db = GraphStat()
            labels = set()
            g_id2label = dict()
            cagetypes = []

            # 基本データの取得
            cell = ice.repcell.mat
            positions = ice.reppositions
            graph = nx.Graph(ice.graph)  # undirected
            ringsize = self.options["ring"]

            # リングの検出
            ringlist = [
                [int(x) for x in ring]
                for ring in cycles_iter(graph, max(ringsize), pos=positions)
            ]
            if not ringlist:
                logger.warning("No rings detected in the structure")
                return False

            ringpos = np.array(
                [centerOfMass(ringnodes, positions) for ringnodes in ringlist]
            )
            logger.info(f"  Rings detected: {len(ringlist)}")

            # ケージの検出
            maxcagesize = max(self.options["sizes"])
            cages = []
            for cage in polyhedra_iter(ringlist, maxcagesize):
                if len(cage) in self.options["sizes"]:
                    valid = True
                    for ringid in cage:
                        if len(ringlist[ringid]) not in ringsize:
                            valid = False
                            break
                    if valid:
                        cages.append(list(cage))

            if not cages:
                logger.warning("No valid cages detected")
                return False

            logger.info(f"  Cages detected: {len(cages)}")
            cagepos = np.array([centerOfMass(cage, ringpos) for cage in cages])

            # 出力形式に応じた処理
            if self.options["gromacs"]:
                self.options["rings"] = ringlist
                self.options["cages"] = cages
                logger.debug(f"Gromacs output prepared: {len(cages)} cages")
                return True

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
                    v = "{0}{1}{2}{3}".format(
                        count[12], count[14], count[15], count[16]
                    )
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
                                v,
                                stat[v] / positions.shape[0],
                                stat[v],
                                positions.shape[0],
                            )
                        )

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
                        label = CageAnalyzer.assign_unused_label(cagesize, labels)
                        g_id2label[g_id] = label
                        labels.add(label)
                        # cage expression
                        index = CageAnalyzer.make_cage_expression(cage, ringlist)
                        logger.info(f"  Cage type: {label} ({index})")
                    else:
                        label = g_id2label[g_id]
                    cagetypes.append(label)
                if len(cagepos) == 0:
                    logger.info("  No cages detected.")
                    print("No cages detected.")
                output = dict(
                    Rings=len(ringlist),
                    Cages=len(cages),
                    details=dict(Counter(cagetypes)),
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

        except Exception as e:
            logger.error(f"Error in Hook2: {str(e)}")
            return False

        logger.info("Hook2: Analysis completed successfully")
        return True

    def Hook6(self, ice):
        """Gromacs形式での出力を行うフック"""
        logger = getLogger()
        logger.info("Hook6: Gromacs format output started")

        try:
            if not hasattr(self.options, "cages") or not hasattr(self.options, "rings"):
                logger.error("Required cage data not found")
                return False

            # 原子データの取得
            atoms = []
            for mols in ice.universe:
                atoms.extend(serialize(mols))

            if not atoms:
                logger.error("No atoms found in the structure")
                return False

            cellmat = ice.repcell.mat
            output = []
            mols = defaultdict(list)

            # 分子データの整理
            for atom in atoms:
                resno, resname, atomname, position, order = atom
                mols[order].append(atom)

            # ケージごとの出力
            for cage in self.options["cages"]:
                cage_output = []
                cage_output.append(
                    "Generated by GenIce https://github.com/vitroid/GenIce\n"
                )

                cagemols = set()
                atomcount = 0
                for ring in cage:
                    cagemols.update(self.options["rings"][ring])

                # 原子データの出力
                for mol in cagemols:
                    for atom in mols[mol]:
                        resno, resname, atomname, position, order = atom
                        cage_output.append(
                            f"{order:5d}{resname:5s}{atomname:>5s}{atomcount + 1:5d}"
                            f"{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}\n"
                        )
                        atomcount += 1

                # ヘッダーとセル情報の追加
                cage_output.insert(1, f"{atomcount}\n")

                if cellmat[1, 0] == 0 and cellmat[2, 0] == 0 and cellmat[2, 1] == 0:
                    cage_output.append(
                        f"    {cellmat[0, 0]} {cellmat[1, 1]} {cellmat[2, 2]}\n"
                    )
                else:
                    assert (
                        cellmat[0, 1] == 0 and cellmat[0, 2] == 0 and cellmat[1, 2] == 0
                    )
                    cage_output.append(
                        f"    {cellmat[0, 0]} {cellmat[1, 1]} {cellmat[2, 2]} "
                        f"{cellmat[0, 1]} {cellmat[0, 2]} {cellmat[1, 0]} "
                        f"{cellmat[1, 2]} {cellmat[2, 0]} {cellmat[2, 1]}\n"
                    )

                output.extend(cage_output)

            print("".join(output), end="")
            logger.info("Hook6: Gromacs format output completed successfully")
            return True

        except Exception as e:
            logger.error(f"Error in Hook6: {str(e)}")
            return False
