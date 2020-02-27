#!/usr/bin/env python3

import json
import sys
import networkx as nx

# read cages in JSON format.
# e.g.
#   % genice CS1 -f cage[json:-16:ring=5-6] > cs1.json
#   % python3 this_program.py cs1.json
filename = sys.argv[1]
with open(filename) as f:
    cageinfo = json.load(f)
cages = cageinfo["cages"]
rings = cageinfo["rings"]

# rebuild a graph from cages and rings

# appearance frequency
g_counts = dict()
for cage in cages:
    # prepare an empty graph
    g = nx.Graph()
    for ringid in cage:
        # add a cycle to the graph
        nx.add_cycle(g, rings[ringid])
    found = None
    for g_known in g_counts:
        # if the graph is already in the dict,
        if nx.is_isomorphic(g_known, g):
            # escape the loop
            found = g_known
            break
    if found is None:
        # g is a new shape
        g_counts[g] = 1
    else:
        # count up the known one.
        g_counts[found] += 1

print(g_counts)
