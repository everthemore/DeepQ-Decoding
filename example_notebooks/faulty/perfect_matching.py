import random
import os
import csv
import subprocess
import time
import copy

import blossom5.pyMatch as pm

## Functions to perform minimum weight matching on the toric and planar topological codes
## There are 4 variants of this here to carry out the matching in the 2D and 3D (imperfect
## measurement) cases, and for the toric and planar codes. Each takes a list of anyon_positions,
## constructs the corresponding graph problem, and interfaces with the Blossom V algorithm (Kologomorov)
## to perform minimum weight matching.


def match_planar_3D(lattice_size,stabilizer_type,anyon_positions,time_space_weights=[1,1],boundary_weight = -1 ,print_graph=False):

    """ Finds a matching to fix the errors in a 3D planar code given the positions of '-1' stabilizer outcomes

    Parameters:
    -----------
    lattice_size -- The dimension of the code
    stabilizer_type -- defines the stabilizer basis, can take the value "star" or "plaquette"
    anyon_positions -- A list of the locations of all '-1' value stabilizers in the 3D parity lattice. [[x0,y0,t0],[x1,y1,t1],...]
    time_space_weights -- The multiplicative weighting that should be assigned to graph edges in the [space,time] dimensions. Default: [1,1]
    boundary_weight -- multiplicative weight to be assigned to edges matching to the boundary. if no boundary_weight specified, set boundary_weight = space_weight
    print_graph -- Set to True to print the constructed graph. Default: False.

    Returns:
    --------
    A list containing all the input anyon positions grouped into pairs. [[[x0,y0,t0],[x1,y1,t1]],[[x2,y2,t2],...

    """

    max_time_separation = 15  # This determines the maximum time separation of edges that are added to the graph
    [wS,wT]=time_space_weights
    wB = wS if boundary_weight == -1 else boundary_weight #if boundary weight not specifiedm, let wB=wS

    total_time=len(anyon_positions)
    nodes_list=[item for sublist in anyon_positions for item in sublist]
    n_nodes=len(nodes_list)

    # exclude edge case where no anyons exist.
    if n_nodes==0:
        return []

    node_index=[]
    count=0

    for x in anyon_positions:
        node_index+=[[count+i for i in range(len(x))]]
        count+=len(x)

    b_node_index=[[index+n_nodes for index in t] for t in node_index]

    all_boundary_nodes=[]
    all_boundary_nodes2=[]


    ## LOOKUP TABLES

    m = 2*lattice_size +1


    weight_lookup={}
    for p in range(-1,m+1):
        weight_lookup[p]={}
        for q in range(-1,m+1):
            weight_lookup[p][q]=wS*abs(p-q)


    ## CONSTRUCT GRAPH
    ## create a graph containing all possible matchings between pairs of anyons (given constraints)
    ## This is represented as three lists:

    nodes1 = []
    nodes2 = []
    weights = []

    ## PART 1: Complete graph between all real nodes

    for i in range(n_nodes -1):
        (pt,p0,p1)=nodes_list[i]

        for j in range(i+1,n_nodes):
            (qt,q0,q1)=nodes_list[j]

            wt=(qt-pt)
            if wt>=max_time_separation: break

            weight = weight_lookup[q0][p0]+weight_lookup[q1][p1]+wt*wT

            nodes1 +=[i]
            nodes2 +=[j]
            weights+=[weight]


    ## PART 2: Generate list of boundary nodes linked to each real node

    boundary_nodes_list = []

    for i in range(n_nodes):

        (pt,p0,p1)=nodes_list[i]

        if stabilizer_type =="star":
            (bt,b0,b1)=(pt,p0,-1 if p1<lattice_size else m)
        elif stabilizer_type=="plaquette":
            (bt,b0,b1)=(pt,-1 if p0<lattice_size else m,p1)
        else:
            print("stabilizer_type must be either *star* or *plaquette*")
            sys.exit(0)

        weight = weight_lookup[p0][b0]+weight_lookup[p1][b1]

        nodes1+=[i]
        nodes2+=[i+n_nodes]
        weights+=[int(weight*wB/wS)]

        boundary_nodes_list+=[(bt,b0,b1)]



 ## PART 3: Complete graph between all boundary nodes

    for i in range(n_nodes -1):
        (pt,p0,p1)=boundary_nodes_list[i]

        for j in range(i+1,n_nodes):
            (qt,q0,q1)=boundary_nodes_list[j]
            wt=(qt-pt)
            if wt>=5: break

            nodes1 +=[n_nodes+i]
            nodes2 +=[n_nodes+j]
            weights+=[0]

    n_edges=len(nodes1)


    ## MAKE MATCHING.
    ## Call the blossom5 perfect matching algorithm to return a matching.
    ## The form of the returned variable <matching> is a list of pairs of node numbers.

    matching = pm.getMatching_fast(2*n_nodes,nodes1,nodes2,weights)


    ## REFORMAT MATCHING PAIRS
    ## Take <matching> and turn it into a list of paired anyon positions.

    matching_pairs=[[i,matching[i]] for i in range(2*n_nodes) if matching[i]>i]

    all_positions=nodes_list+boundary_nodes_list

    points=[] if len(matching_pairs)==0 else [[all_positions[i] for i in x] for x in matching_pairs]

    return points
