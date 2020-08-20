# TSP SOLVER USING LINEAR PROGRAMMING
# gurobi tsp solver is based on gurobi.github.io/modeling-examples/traveling_salesman/tsp.html
# By Kirill Sechkar
# v0.1.1, 22.7.20

import numpy as np
import time

import gurobipy as gp
from gurobipy import GRB
from itertools import combinations,product


#-------------------TSP SOLVER (NON-CAPACITATED PROBLEM)--------------------
# (for infinite capacities)
#solver
def tsp_lp_gurobi(D):
    # PART 0: technical GUROBI works
    # set output flag to 0 to prevent GUROBI from printing out logs
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()


    # PART 1: initial preparations

    # PART 1.1: get the auxiliary array of node indices
    global nodes
    nodes = []
    for i in range(0, len(D)):
        nodes.append(i)

    # PART 1.2: convert distance matrix into the dictionary format used by gurobi
    dist = {(i, j): D[i][j] for i, j in product(nodes, nodes) if i != j}

    # PART 1.3: create the gurobi model
    m = gp.Model(env=env)

    # PART 1.4: create a boolean matrix indicating the trip (commonly known as X in literature)
    vars = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='e')


    # PART 2: add initial constraints
    m.addConstrs(vars.sum(c, '*') == 1 for c in nodes)  # each node has 1 incoming edge
    m.addConstrs(vars.sum('*', c) == 1 for c in nodes)  # each node has 1 outgoing edge


    # PART 3: optimise the model
    m._vars = vars
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)  # optimise while adding lazy subtour-eliminating constraints


    # PART 4: reconstruct tour from m.vars
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)  # get the edges selected as tour
    tour = subtour(selected)  # get tour from selected edges

    assert len(tour) == len(nodes)  # sanity check that the tour is actually complete

    return tour


# callback function, uses lazy constraints to eliminate sub-tours
def subtourelim(model, where):
    if where == GRB.Callback.MIPSOL:
        # make a list of edges selected in the solution
        vals = model.cbGetSolution(model._vars)
        selected = gp.tuplelist((i, j) for i, j in model._vars.keys() if vals[i, j] > 0.5)
        # find the shortest cycle (subtour) in the selected edge list
        tour = subtour(selected)
        if len(tour) < len(nodes):
            # add subtour elimination constraint for the subtour
            model.cbLazy(
                gp.quicksum(model._vars[i, j] + model._vars[j, i] for i, j in combinations(tour, 2)) <= len(tour) - 0.5)


# find the shortest subtour among edges selected by the solution
# for the final solution, this is the tour that has to be found
def subtour(edges):
    unvisited = nodes[:]
    cycle = nodes[:]  # dummy - will be replaced
    while (len(unvisited)!=0):
        thiscycle = []
        neighbours = unvisited
        while (len(neighbours)!=0):
            current = neighbours[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbours = [j for i, j in edges.select(current, '*')
                         if j in unvisited]
        if len(thiscycle) <= len(cycle):
            cycle = thiscycle  # new shortest subtour
    return cycle


#-----------------LP PROBLEM SOLVER (CAPACITATED PROBLEM)-------------------
#solver
def lp_cap(D, cap, maxtime):
    # PART 0: technical GUROBI works

    # PART 0.1: set output flag to 0 to stop gurobi from printing out the log
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()

    # PART 0.2: set maximum optimisation time
    if (maxtime == None):
        # default value:
        mxt = 1  # optimisation time of more than 1s never lead to improved results, so is unnecessary
    else:
        mxt = maxtime

    # PART 1: initial preparations

    # PART 1.1: get the auxiliary array of node indices
    # note the 0 node is skipped - it's an auxiliary only used by the TSP solver
    global nodes
    nodes = []
    for i in range(1, len(D)):
        nodes.append(i)

    # PART 1.2: convert distance matrix into the dictionary format used by GUROBI
    dist = {(i, j): D[i][j] for i, j in product(nodes, nodes) if i != j}

    # PART 1.3: create the GUROBI model
    m = gp.Model(env=env)

    # PART 1.4: create a boolean matrix indicating the trip (commonly known as X in literature)
    # Note: obj set to -1 as we need to MAXIMISE the number of Xij=1, i.e. minimise the sum of -Xij
    vars = m.addVars(dist.keys(), vtype=GRB.BINARY, name='e')

    # PART 1.5: copy capacity into a global variable to let other functions use it
    global gurcap
    gurcap = cap


    # PART 2: add constraints

    # PART 2.1: Limit number of incoming and outgoing edges for each node
    m.addConstrs(vars.sum(c, '*') <= 1 for c in nodes)  # each node has at most 1 incoming edge
    m.addConstrs(vars.sum('*', c) <= 1 for c in nodes)  # each node has at most 1 outgoing edge

    # PART 2.2: eliminate cycles and chains that are too long
    # create dummy variables needed for constraints
    u = {}
    for i in nodes:
        u[i] = m.addVar(lb=0, ub=gurcap, vtype="C", name="u(%s)" % i)

    # add the constraints
    if (gurcap > 1.5):
        m.addConstrs(u[j] - u[i] >= 1 - gurcap * (1 - vars[i, j]) for i, j in combinations(nodes, 2))
        m.addConstrs(u[i] - u[j] >= 1 - gurcap * (1 - vars[j, i]) for i, j in combinations(nodes, 2))

    # PART 2.3: 2-hop loop elimination
    m.addConstr(
        gp.quicksum(vars[i, j] * dist[(i, j)] + vars[j, i] * dist[(j, i)] for i, j in combinations(nodes, 2)) == 0)


    # PART 3: optimise the model
    m.Params.TIME_LIMIT = mxt # set optimisation time limit
    m.setObjective(gp.quicksum(vars[i, j] + vars[j, i] for i, j in combinations(nodes, 2)), GRB.MAXIMIZE) #set objective
    m.optimize() #optimise


    # PART 4: reconstruct the chains from m.vars (matrix X in literature)
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)  # get the edges selected as tour
    chains, iscycle = recover(selected)  # get tour from selected edges

    # sanity check that there are no cycles
    for i in range(0, len(chains)):
        if (iscycle[i]):
            print('Error! Cycle detected!')
            break

    return chains


# recover the chains
def recover(edges):
    # PART 1: initial preparations
    iscycle = []  # which of chains are cycles (True for a cycle)
    chains = []  # all the chains covering the cycle
    unvisited = np.ones((len(nodes)+1),dtype=bool)  #tells if a node was visited, position[0] is not needed
    unvisited[0]=False  # first False just ensures compatibility with early versions of program)

    # PART 2: get the chains
    while True:
        # see if no more unvisited nodes are left
        noneleft=True
        for i in range(1,len(unvisited)):
            if (unvisited[i]):
                noneleft=False
                break

        # proceed to return statement if indeed no unvisited nodes are left
        if(noneleft):
            break

        # otherwise, start reconstructing the chain which the unvisited node we found belongs to
        thischain = [nodes[i - 1]]
        unvisited[i] = False # and NOW, the node has been visited

        # now we need to find the nodes that go before and after it in the chain
        # first, get the nodes after
        while True:
            # get the next node in chain
            next = edges.select(thischain[-1], '*')

            # if the chain's end was reached, this isn't a cycle
            if(next==[]):
                iscycle.append(False)
                break
            # if the node we began from is encountered again, it's a cycle and we have walked it all
            elif(next[0][1]==thischain[0]):
                iscycle.append(True)
                break
            # if neither, just add the next node in the chain
            else:
                thischain.append(next[0][1])
                unvisited[next[0][1]]=False

        # now get the previous nodes (if it wasn't a cycle)
        if not (iscycle[-1]):
            # get previous node in the chain (noted as 'next' to draw a parallel with the previous part)
            next = edges.select('*',thischain[0])
            # keep getting previous nodes until there are none
            while (next!=[]):
                thischain.insert(0,next[0][0])
                unvisited[next[0][0]] = False
                next = edges.select('*', thischain[0])

        # record the obtained chain
        chains.append(thischain.copy())

    return chains, iscycle


#------------------------------MAIN (TESTING ONLY)--------------------------
def main():
    A=[[0,0,0,0],
       [0,0,0,1],
       [0,0,0,0],
       [0,0,0,0]]
    surelymore = 0
    for i in range(0, len(A)):
        for j in A[i]:
            surelymore += j
    for i in range(0, len(A)):
        A[i][i] = surelymore
    # t1=time.time()
    print(lp_cap(A,2,0.1))
    # print(time.time()-t1)


if __name__ == "__main__":
        main()