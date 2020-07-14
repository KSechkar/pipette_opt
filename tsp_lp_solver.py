# TSP SOLVER USING LINEAR PROGRAMMING
# cvxpy parts based on gist.github.com/AshNguyen/3bbaef5203f574cc5a2bad255ca59069
#pure gurobi based on gurobi.github.io/modeling-examples/traveling_salesman/tsp.html
# By Kirill Sechkar
# v0.0.1, 7.7.20

import cvxpy as cvx
import numpy as np
import time

import gurobipy as gp
from itertools import combinations,product


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


#-----------------PURE GUROBI------------
"""def tsp_lp(D):
    nodes=[]
    for i in range(0,len(D)):
        nodes.append(i)

    dist = {(i, j): D[i][j] for i, j in product(nodes, nodes) if i != j}

    m = gp.Model()

    # variables telling if node 'i' is adjacent to node 'j' on the tour
    vars = m.addVars(dist.keys(), obj=dist, vtype=gp.GRB.BINARY, name='e')
    for i, j in vars.keys(): #REMOVE ASAP
        m.addConstr(vars[j, i] == vars[i, j])  # edge in opposite direction

    # Constraints: two edges incident to each city
    m.addConstrs(vars.sum(c, '*') == 2 for c in nodes)

    # Optimize the model
    m._vars = vars
    m.Params.lazyConstraints = 1
    m.optimize(subtourelim)

    #get tour
    vals = m.getAttr('x', vars)
    selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)

    tour = subtour(selected,nodes)
    assert len(tour) == len(nodes)
    print(tour)
    return 0

# callback - lazy subtour elimination
def subtourelim(model, where,nodes):
    if where == gp.GRB.Callback.MIPSOL:
        # make a list of edges selected in the solution
        vals = model.cbGetSolution(model._vars)
        selected = gp.tuplelist((i, j) for i, j in model._vars.keys()
                             if vals[i, j] > 0.5)
        # find the shortest cycle in the selected edge list
        tour = subtour(selected)
        if len(tour) < len(nodes):
            # add subtour elimination constr. for every pair of cities in subtour
            model.cbLazy(gp.quicksum(model._vars[i, j] for i, j in combinations(tour, 2))
                         <= len(tour)-1)

# Given a tuplelist of edges, find the shortest subtour
def subtour(edges,nodes):
    unvisited = nodes[:]
    cycle = nodes[:] # Dummy - guaranteed to be replaced
    while unvisited:  # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        while neighbors:
            current = neighbors[0]
            thiscycle.append(current)
            unvisited.remove(current)
            neighbors = [j for i, j in edges.select(current, '*')
                         if j in unvisited]
        if len(thiscycle) <= len(cycle):
            cycle = thiscycle # New shortest subtour
    return cycle"""


#-----------------MAIN (TESTING ONLY)------------
def main():
    A = [[0, 60, 79, 37, 10, 61],
         [50, 0, 22, 48, 63, 54],
         [79, 42, 0, 49, 70, 38],
         [37, 48, 49, 0, 38, 45],
         [10, 63, 70, 38, 0, 53],
         [61, 54, 38, 45, 53, 0]]
    surelymore = 0
    for i in range(0, len(A)):
        for j in A[i]:
            surelymore += j
    for i in range(0, len(A)):
        A[i][i] = surelymore
    #t1=time.time()
    print(tsp_lp_mtz(A))
    #print(time.time()-t1)


if __name__ == "__main__":
        main()