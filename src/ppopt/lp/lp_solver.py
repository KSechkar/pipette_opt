# LP SOLVER USING LINEAR PROGRAMMING
# the tsp solver function is based on gurobi.github.io/modeling-examples/traveling_salesman/tsp.html
# By Kirill Sechkar
# v0.1.10, 30.5.21

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
    m.addConstrs(vars.sum(c,'*') == 1 for c in nodes)  # each node has 1 incoming edge
    m.addConstrs(vars.sum('*',c) == 1 for c in nodes)  # each node has 1 outgoing edge


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

    # PART 1: initial preparations
    # PART 1.1: get the auxiliary array of node indices
    nodes = []
    for i in range(0, len(D)):
        nodes.append(i)
    wellnodes = nodes[1:len(nodes)] # record only nodes matched to wells (i.e. except 0)
    
    # PART 1.2: convert distance matrix into the dictionary format used by GUROBI
    dist = {(i, j): D[i][j] for i, j in product(nodes, nodes) if i != j}

    # PART 1.3: create the GUROBI model
    m = gp.Model(env=env)

    # PART 1.4: create a boolean matrix indicating the trip (commonly known as X in literature)
    # Note: the objective will be to minimise number of chosen edges leaving the 'depot', i.e. sum X_0,c for all c
    vars = m.addVars(dist.keys(), vtype=GRB.BINARY, name='e')

    # PART 1.5: copy capacity into a global variable to let other functions use it
    global gurcap
    gurcap = cap


    # PART 2: add constraints

    # PART 2.1: Limit number of incoming and outgoing edges for each node
    m.addConstrs(vars.sum(c, '*') == 1 for c in wellnodes)  # each well node has at most 1 incoming edge
    m.addConstrs(vars.sum('*', c) == 1 for c in wellnodes)  # each well node has at most 1 outgoing edge

    # PART 2.2: eliminate cycles and chains that are too long
    # create dummy variables needed for constraints
    u = {}
    for i in wellnodes:
        u[i] = m.addVar(lb=0, ub=gurcap, vtype="C", name="u(%s)" % i)

    # add the constraints
    if (gurcap > 1.5):
        m.addConstrs(u[j] - u[i] >= 1 - gurcap * (1 - vars[i, j]) for i, j in combinations(wellnodes, 2))
        m.addConstrs(u[i] - u[j] >= 1 - gurcap * (1 - vars[j, i]) for i, j in combinations(wellnodes, 2))

    # PART 2.3: all edges selected into chains must be zero-length
    m.addConstr(
        gp.quicksum(vars[i, j] * dist[(i, j)] + vars[j, i] * dist[(j, i)] for i, j in combinations(nodes, 2)) == 0)

    # PART 3: optimise the model
    m.Params.TIME_LIMIT = maxtime # set optimisation time limit
    m.setObjective(gp.quicksum(vars[0, c] for c in wellnodes), GRB.MINIMIZE) #set objective
    m.optimize() #optimise


    # PART 4: reconstruct the chains from m.vars (matrix X in literature)
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)  # get the edges selected as tour

    return recover(selected)  # get tour from selected edges


# recover the chains
def recover(edges):
    chains=[] # make an empty list of 'chains'

    # go through all routes starting at the 0 ('depot') node
    for edgeout in edges.select(0,'*'):
        curnode = edgeout[1] # current node is the starting node of the route
        curchain=[curnode] # record starting node into chain

        # get the rest of the chain
        while True:
            curnode = edges.select(curnode,'*')[0][1] # current node is now the next node
            if(curnode!=0): # if it's a well node, record into chain
                curchain.append(curnode)
            else: # if returning to the 'depot', the whole route has been recorded
                break

        # record the obtained chain into the output list
        chains.append(curchain.copy())
    return chains


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
    t1=time.time()
    print(lp_cap(A,2,0.1))
    print(time.time()-t1)


if __name__ == "__main__":
    import time
    main()
