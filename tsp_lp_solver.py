# TSP SOLVER USING LINEAR PROGRAMMING
# based on gist.github.com/AshNguyen/3bbaef5203f574cc5a2bad255ca59069
# By Kirill Sechkar
# v0.0.1, 7.7.20

import cvxpy as cvx
import numpy as np
import time
import gurobipy


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
def tsp_lp(D, solving):
    # get number of nodes
    numnodes = len(D)
    # matrix indicating the trip
    if(solving): # if we are solving the TSP, we need X to be boolean
        X = cvx.Variable((numnodes, numnodes), boolean=True)
    # otherwise, we are trying to estimate the lower bound using a relaxed problem (non-boolean X)
    else:
        X = cvx.Variable((numnodes, numnodes), nonneg=True)

    # objective
    obj = cvx.Minimize(sum(sum(cvx.multiply(D, X))))

    # basic condition
    Ones=np.ones(numnodes)
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
    if(solving):
        return recover(X.value)
    else:
        return opt

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
    print(tsp_lp(A,True))
    #print(time.time()-t1)


if __name__ == "__main__":
        main()