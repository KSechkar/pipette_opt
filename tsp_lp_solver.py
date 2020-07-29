# TSP SOLVER USING LINEAR PROGRAMMING
# cvxpy parts based on gist.github.com/AshNguyen/3bbaef5203f574cc5a2bad255ca59069
# pure gurobi based on gurobi.github.io/modeling-examples/traveling_salesman/tsp.html
# By Kirill Sechkar
# v0.1.1, 22.7.20

import cvxpy as cvx
import numpy as np
import time

import gurobipy as gp
from gurobipy import GRB
from itertools import combinations,product

"""
#------------------------USING CVXPY --------------------
# recovers the tour from the boolean matrix X encoding it
def recover(Xval):
    node = 0
    tour = [node]
    while True:
        for j in range(len(Xval)):
            if Xval[j, tour[-1]] > 0.99:
                if j not in tour:
                    tour.append(j)
                else:
                    return tour


# solver function (uses Miller–Tucker–Zemlin formulation)
def tsp_lp_mtz(D):
    # get number of nodes
    numnodes = len(D)
    # matrix indicating the trip
    X = cvx.Variable((numnodes, numnodes), boolean=True)

    # objective
    obj = cvx.Minimize(sum(sum(cvx.multiply(D, X))))

    # basic condition
    Ones=np.ones(numnodes,dtype=int)
    constraints = [(cvx.sum(X, axis=0) == Ones), (cvx.sum(X, axis=1) == Ones)]

    # subtour elimination
    u = cvx.Variable(numnodes,integer=True)
    for i in range(1, numnodes):
        for j in range(1, numnodes):
            if i != j:
                constraints.append(u[i] - u[j] + (numnodes-1) * X[i, j] <= (numnodes - 2))

    for i in range(numnodes):
        constraints.append(u[i] >= 0)
        constraints.append(u[i] <= (numnodes-1))

    # solving the problem (prob.solve returns cost)
    prob = cvx.Problem(obj, constraints)
    opt = prob.solve(solver=cvx.GUROBI)
    return recover(X.value)


# solver function (uses lazy subtour elimination)
def tsp_lp_lazy(D):
    # get number of nodes
    numnodes = len(D)
    # matrix indicating the trip
    X = cvx.Variable((numnodes, numnodes), boolean=True)

    # objective
    obj = cvx.Minimize(sum(sum(cvx.multiply(D, X))))

    # basic condition
    Ones = np.ones(numnodes)
    constraints = [(cvx.sum(X, axis=0) == Ones), (cvx.sum(X, axis=1) == Ones)]

    # preliminary solution, which might involve subtours
    prob = cvx.Problem(obj, constraints)
    time1 = time.time()
    opt = prob.solve(solver=cvx.GUROBI)

    while True:
        subtour = recover(X.value)
        if (len(subtour) == numnodes):
            #print("Minimal time: ", opt)
            #print("Optimal tour: ", subtour)
            #print("Converge time: ", time.time() - time1)
            return subtour

        else:
            #print("Try: ", subtour)
            nots = [j for j in range(numnodes) if j not in subtour]
            constraints.append(sum(X[i, j] for i in subtour for j in nots) >= 1)
            prob = cvx.Problem(obj, constraints)
            opt = prob.solve(solver=cvx.MOSEK)
"""


#-----------------PURE GUROBI------------
#solver
def tsp_lp_gurobi(D):
    # PART 1: initial preparations
    # the array of node indices, is auxiliary
    global nodes
    nodes = []
    for i in range(0, len(D)):
        nodes.append(i)

    # convert distance matrix into the dictionary format gurobi uses
    dist = {(i, j): D[i][j] for i, j in product(nodes, nodes) if i != j}

    # set output falg to 0 to stop gurobi from printing
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()

    # create gurobi model
    m = gp.Model(env=env)

    # create matrix indicating the trip (commonly known as X in literature)
    vars = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='e')

    # PART 2: add initial constraints
    m.addConstrs(vars.sum(c, '*') == 1 for c in nodes)  # each node has 1 incoming edge
    m.addConstrs(vars.sum('*', c) == 1 for c in nodes)  # each node has 1 outgoing edge\

    # PART 3: optimise the model
    m._vars = vars
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)  # optimise while adding lazy subtour-eliminating constraints

    # PART 4: reconstruct tour from m.vars (matrix X in literature)
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)  # get the edges selected as tour
    tour = subtour(selected)  # get tour from selected edges

    assert len(tour) == len(nodes)  # sanity check that the tour is actually complete

    """
    TESTING ONLY (needs tspy package):
    tsp = TSP()
    tsp.read_mat(D)
    two_opt = TwoOpt_solver(initial_tour='NN', iter_num=100)
    tour1 = tsp.get_approx_solution(two_opt)
    print(tour1)
    print(get_cost(tour1, tsp))
    """

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


# find the shortest subtour among edges selected by the solution (for a final solution, this is our desired tour)
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


#-----------------GUROBI WITH CAPACITY------------
#solver
def lp_cap(D,cap,maxtime):
    # PART 1: initial preparations
    # set default maximum optimisation time if none nad-set
    if(maxtime==None):
        maxtime=30.0

    # the array of node indices, is auxiliary
    global nodes
    nodes = []
    for i in range(1, len(D)):
        nodes.append(i)

    # copy capacity into global variable to let other functions use it
    global gurcap
    gurcap = cap

    # convert distance matrix into the dictionary format gurobi uses
    dist = {(i, j): D[i][j] for i, j in product(nodes, nodes) if i != j}

    #set output falg to 0 to stop gurobi from printing
    env=gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()

    # create gurobi model
    m = gp.Model(env=env)

    # create matrix indicating the trip (commonly known as X in literature)
    # obj set to -1 as we need to MAXIMISE the number of Xij=`1
    vars = m.addVars(dist.keys(), vtype=GRB.BINARY, name='e')

    # dummy variables used for constraints
    u = {}
    for i in nodes:
        u[i] = m.addVar(lb=0, ub=gurcap, vtype="C", name="u(%s)" % i)

    # PART 2: add initial constraints
    m.addConstrs(vars.sum(c, '*') <= 1 for c in nodes)  # each node has at most 1 incoming edge
    m.addConstrs(vars.sum('*', c) <= 1 for c in nodes)  # each node has at most 1 outgoing edge
    # so that each chain is a single pipette's journey...
    m.addConstr(
        gp.quicksum(vars[i, j] * dist[(i, j)] + vars[j, i] * dist[(j, i)] for i, j in combinations(nodes, 2)) == 0)

    # cycle and too-long chain elimination
    if (gurcap > 1.5):
        m.addConstrs(u[j] - u[i] >= 1 - gurcap * (1 - vars[i, j]) for i, j in combinations(nodes, 2))
        m.addConstrs(u[i] - u[j] >= 1 - gurcap * (1 - vars[j, i]) for i, j in combinations(nodes, 2))

    # PART 3: optimise the model
    m.Params.TIME_LIMIT = maxtime
    m.setObjective(gp.quicksum(vars[i, j] + vars[j, i] for i, j in combinations(nodes, 2)), GRB.MAXIMIZE)
    m.optimize()  # optimise while adding lazy subtour-eliminating constraints

    # PART 4: reconstruct tour from m.vars (matrix X in literature)
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)  # get the edges selected as tour
    chains, iscycle = recover(selected)  # get tour from selected edges

    for i in range(0, len(chains)):
        if (iscycle[i] == True):
            print('Error! Cycle detected!')
            break
    # assert len(tour) == len(nodes)  # sanity check that the tour is actually complete

    """
    TESTING ONLY (needs tspy package):
    tsp = TSP()
    tsp.read_mat(D)
    two_opt = TwoOpt_solver(initial_tour='NN', iter_num=100)
    tour1 = tsp.get_approx_solution(two_opt)
    print(tour1)
    print(get_cost(tour1, tsp))
    """

    return chains


# find the shortest subtour among edges selected by the solution (for a final solution, this is our desired tour)
def recover(edges):
    iscycle = []  # which of tours are cycles
    chains = []  # all the chains covering the cycle
    unvisited = np.ones((len(nodes)+1),dtype=bool)  #tells if a node was visited, position[0] not needed
    unvisited[0]=False  # first False just ensures compatibility with early versions of program)
    while True:
        # see if no more unvisited nodes are left and if so, quit
        noneleft=True
        for i in range(1,len(unvisited)):
            if (unvisited[i]):
                noneleft=False
                unvisited[i]=False
                break
        if(noneleft):
            break

        # start getting a chain using the first found unvisited node
        thischain=[nodes[i-1]]

        while True:
            # get next node in chain
            next = edges.select(thischain[-1], '*')
            if(next==[]):  # if we got to the end, this isn't a cycle
                iscycle.append(False)
                break
            elif(next[0][1]==thischain[0]):  # if we got back to the beginning, it's a cycle and we've walked it all
                iscycle.append(True)
                break
            else:  # if neither, just add the regular next node in list
                thischain.append(next[0][1])
                unvisited[next[0][1]]=False

        if not (iscycle[-1]): # if it's not a cycle, recover its beginning
            next = edges.select('*',thischain[0])
            while (next!=[]):
                thischain.insert(0,next[0][0])
                unvisited[next[0][0]] = False
                next = edges.select('*', thischain[0])

        chains.append(thischain.copy())

    return chains, iscycle


#-----------------MAIN (TESTING ONLY)------------
def main():
    """
    A = [[0, 60, 79, 37, 10, 61],
         [50, 0, 22, 48, 63, 54],
         [79, 42, 0, 49, 70, 38],
         [37, 48, 49, 0, 38, 45],
         [10, 63, 70, 38, 0, 53],
         [61, 54, 38, 45, 53, 0]]
    """
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
    print(gurobi_cap_nolazy(A,2))
    # print(time.time()-t1)


if __name__ == "__main__":
        main()