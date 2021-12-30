# from numpy import load
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from multiprocessing import Pool


DLT_RAW = 'DLT triangulation un-normalized'
DLT_nN = 'DLT triangulation normalized with N'
DLT_nK = 'DLT triangulation normalized with inv(K)'
GT = 'ground truth 3D points'
SOC_RAW = 'SOCP triangulation un-normalized'
SOC_nN = 'SOCP triangulation normalized with N'
SOC_nK = 'SOCP triangulation normalized with inv(K)'

MRE = {DLT_RAW: 10.76624334,
       DLT_nN: 49038,
       DLT_nK: 0.32154451,
       GT: 0.3204574579,
       SOC_RAW: 2.45645747,
       SOC_nN: None,
       SOC_nK: 2.138747459}

# #
# #
# # def load_data(file='Drinking Fountain Somewhere In Zurich.npz'):
# #     data = load(file)
# #     # add ones to x
# #     x = data['x']
# #     ones = np.ones((x.shape[0], 1, x.shape[2]))
# #
# #     ret = {item: data[item] for item in data.files}
# #     ret['x'] = np.concatenate((x, ones), axis=1)
# #     return ret
# #
# #
# # data = load_data()
# #
# import matplotlib.pyplot as plt
# import numpy as np
#
# data =  [
#             [                'Freeze', 'Wind', 'Flood', 'Quake', 'Hail'],
#             [ 'explanation',  1, 2,   3,  4,  5],
#             # ['10 year',  58230, 381139,   78045,   99308, 160454],
#             # ['20 year',  89135,  80552,  152558,  497981, 603535],
#             # ['30 year',  78415,  81858,  150656,  193263,  69638],
#             # ['40 year', 139361, 331509,  343164,  781380,  52269],
#         ]
# # Pop the headers from the data array
# column_headers = data.pop(0)
# row_headers = [x.pop(0) for x in data]
#
# fig, ax = plt.subplots()
# ax.set_axis_off()
#
# the_table = ax.table(cellText=data,
#                       rowLabels=row_headers,
#                       rowColours= ["palegreen"] * len(row_headers),
#                       colColours= ["palegreen"] * len(column_headers),
#                       colLabels=column_headers,
#                       loc='center')
# ax.set_title('TITLE',
#              fontweight="bold")
#
# # plt.show()
#
#
# x = np.array([1, 2, 3])
# y = np.array([2, 3, 4])
#
# res = np.einsum('i,j->ij', x, y)
#
# c = 0



# # Generate a random feasible SOCP.
# m = 3  # number of constraints
# n = 10  # size of the result vector
# p = 5
# n_i = 5  # A's rows dimension
# np.random.seed(2)
# f = np.random.randn(n)
# A = []
# b = []
# c = []
# d = []
# x0 = np.random.randn(n)
# for i in range(m):
#     A.append(np.random.randn(n_i, n))
#     b.append(np.random.randn(n_i))
#     c.append(np.random.randn(n))
#     d.append(np.linalg.norm(A[i] @ x0 + b, 2) - c[i].T @ x0)
# F = np.random.randn(p, n)
# g = F @ x0
#
# # Define and solve the CVXPY problem.
# x = cp.Variable((n,))
# # We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
# soc_constraints = [
#     cp.SOC(c[i].T @ x + d[i], A[i] @ x + b[i]) for i in range(m)
# ]
# prob = cp.Problem(cp.Minimize(f.T @ x),
#                   soc_constraints + [F @ x == g])
# prob.solve()
# print(x.value)

# v = x.value
# for i in range(m):
#     Ai = A[i]
#     ci = c[i]
#
#     left = np.linalg.norm(A[i] @ v + b[i])
#     right = c[i].T @ v + d[i]
#     print(f"||A{i} @ v + b{i}||<= c{i}.T @ v + d{i}", f"{left} <= {right}", "right - left=", right - left)
# # cp.SOC(c[i].T @ x + d[i], A[i] @ x + b[i]) for i in range(m)


# f = np.zeros(5)
# f[-1] = 1
#
# n = 3
# A = np.array([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]])
# B = np.zeros((n,n+2))
# B[:,:-1] = A
# # print(A)
# # print(B)
#
# gamma = 2
# n = 4
# c = np.random.randn(n)
# ai = np.ones(n+1)  # 5
# ai[0:-1] = gamma * c
# # print(c)
# # print(ai)
#
#
# x = np.random.randn(4, 3)
# x = np.expand_dims(x, axis=2)
# print(x.shape)


# # 2e - plot bar chart
#
# mre1, mre2, mre3, mre4, mre5 = 0.32, 10, 0.32, 0.35, 1.35
# # creating the dataset
# data = {'GT': mre1, 'DLT Unnormalized': mre2, 'DLT Normalized': mre3,
#         'SOCP Unnormalized': mre4, 'SOCP Normalized': mre5}
# methods = list(data.keys())
# mres = list(data.values())
#
# fig, ax = plt.subplots(figsize=(10, 5))
#
# # creating the bar plot
# bar = ax.bar(methods, mres, color='maroon', width=0.4)
# plt.ylabel("MRE")
# plt.title("Mean Reprojection Error")
# # fig.autofmt_xdate()
#
# for index, data in enumerate(mres):
#     ax.text(x=index , y=data , s=f"{data}" ,ha='center', va='bottom')
#
# plt.show()


# n_points = 5230
# groups = [(i, min(i+540, n_points)) for i in range(0, n_points, 540)]
# print(groups)



def show_bar(title):
    # creating the dataset
    data = {'GT': MRE[GT], 'DLT': MRE[DLT_RAW], 'DLT Normalized': MRE[DLT_nK],
            'SOCP': MRE[SOC_RAW], 'SOCP Normalized': MRE[SOC_nK]}
    methods = list(data.keys())
    mres = list(data.values())

    fig, ax = plt.subplots(figsize=(10, 5))

    # creating the bar plot
    bar = ax.bar(methods, mres, color='maroon', width=0.4)
    plt.ylabel("MRE")
    plt.xlabel("triangulation scheme")
    plt.title(title)
    # fig.autofmt_xdate()

    for index, data in enumerate(mres):
        ax.text(x=index, y=data, s=f'{data:.4f}', ha='center', va='bottom')
    plt.show()

show_bar("2.e Mean Reprojection Error after Triangulation")



######## DEBUG
# max_iter = 30
def job(args):
    s, e = args
    # points from s to e
    print(f"process {s} to {e}")
    Xs = [SOCP_triangulate(x[visible_points[:, j], :, j], P[visible_points[:, j], :, :], max_iter=30)
          for j in range(s, e)]
    # for j in range(s, e):
    #     X = SOCP_triangulate(x[visible_points[:, j], :, j], P[visible_points[:, j], :, :], max_iter)
    #     Xs.append(X)
    return {'start': s,
            'end': e,
            'result': Xs}


if __name__ == '__main__':
    pass
    # # A, c = get_triangulation_parameters(P, x)  # A: [n_cam, n_points, 3, 4]
    # # n_cam, n_points, _, _ = A.shape
    # #
    # # # triangulate for each point
    # # gamma = 1
    # # for j in range(n_points):
    # #     A_see = A[visible_points[:, j], j, :, :]  # shape [n_see_cam, 3, 4]
    # #     c_see = c[visible_points[:, j], :]  # shape [n_see_cam, 4]
    # #     Xj = solve_SOCP(A_see,c_see, gamma)


    # normalize
    # K_inv = np.array([np.linalg.inv(k) for k in K])
    # P, x = normalize_cam_points(P,x,K_inv)

    # triangulate for each point

    # for j in range(n_points):
    #     X = SOCP_triangulate(x[visible_points[:, j], :, j], P[visible_points[:, j], :, :], max_iter)
        # print("finish ", j)

    # ##### Parallel Computation
    # g = int(n_points / 10)  # group size
    # groups = [(i, min(i+g, n_points)) for i in range(0, n_points, g)]
    # print("Parallel Computation: \n\tnumber of points", n_points, "\n\tgroup size", g, "\n\tnumber of groups", len(groups))
    #
    # with Pool(3) as p:
    #     data = p.map(job, groups)
    #     data = {d['start']: d['result'] for d in data}
    #     # order the points
    #     points3d = []
    #     for g in groups:
    #         points3d += data[g[0]]


    # exit(0)
######## DEBUG - END


def dana_check(A,c, gamma):
    try:
        n_cam = A.shape[0]  # number of constraints = n_cam
        x = cp.Variable((4,))  # [X s]
        s = cp.Variable((1,))
        # We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
        soc_constraints = [cp.SOC(gamma * c[i, :].T @ x, A[i, :, :] @ x) for i in range(n_cam)]
        obj = cp.Minimize(s)
        prob = cp.Problem(obj, soc_constraints)
        prob.solve()

        stat = prob.status
        val = prob.value
        X = x.value
        S = s.value

        b = 0

    except Exception as e:
        print("Exception:", str(e))
        return None



def solve_SOCP_old(A,c, gamma):
    """
    Solve SOC feasibility problem FOR A SINGLE POINT todo.
    Use cvxpy to solve a second order cone program of the form:
    minimize f^Tx s.t ||B[i]x + b[i]|| <= a[i]^Tx+d[i] for i=1,...
    :param A: a list of numpy arrays,    for a single point, shape [n_see_cam, 3, 4] todo
    :param c: a list of numpy arrays                         shape [n_see_cam, 4] todo
    :return: If there is a solution return it. else, return None.
     A solution is an array of size [4,1]. Make sure it's least coordinate is 1.

    Hint: use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
    """
    try:
        # SOCP setup
        n_cam = A.shape[0]  # number of constraints = n_cam
        n = c.shape[1] + 1  # size of the result vector = 5
        r = A.shape[1]  # A's rows dimension = 3
        col = A.shape[2]  # 4
        f = np.zeros((n, ))
        f[-1] = 1  # f = [0,0,0,0,1]
        B, a = [], []
        # loop over the cameras
        for i in range(n_cam):
            Bi = np.zeros((r, col + 1))  # shape [3, 5]
            Bi[:, :-1] = A[i, :, :]  # Bi = [Ai 0]
            B.append(Bi)

            ai = np.ones((n,))  # shape [5, ]
            ai[0:-1] = gamma * c[i, :]  # ai = [gamma*ci, 1]
            a.append(ai)

        # Define and solve the CVXPY problem
        x = cp.Variable((n,))  # [X s]
        # We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
        soc_constraints = [cp.SOC(a[i].T @ x, B[i] @ x) for i in range(n_cam)]
        obj = cp.Minimize(f.T @ x)
        prob = cp.Problem(obj, soc_constraints)
        prob.solve()

        ########################
        def find_X():
            # objective: find X s.t. for all i: ||Ai @ X|| <= gamma ci.T @ X
            x = cp.Variable((4,))
            soc_constraints = [cp.SOC(gamma * c[i, :].T @ x, A[i, :, :] @ x) for i in range(n_cam)]
            obj = cp.Minimize(cp.Constant((1,)))   # cp.expressions.constants.Constant
            prob = cp.Problem(obj, soc_constraints)
            prob.solve()

            ##### VERIFY THE SOLUTION
            # print("status:", prob.status)
            X = x.value
            # eps = 1e-6
            # print("v", v)
            # if v[0] < eps and v[1] < eps and v[2] < eps and v[3] < eps:
            #     print("solution is zero")
            #     v[0] = 0
            #     v[1] = 0
            #     v[2] = 0
            #     v[3] = 0

            # for i in range(n_cam):
            #     Ai = A[i, :, :]
            #     ci = c[i, :]
            #     left = np.linalg.norm(Ai @ X)
            #     right = gamma * ci.T @ X
            #     print(f"||A{i} @ X||<=gamma * c{i}.T @ X", f"{left} <= {right}", "right - left=", right - left)

            return utils.pflat(X), X

        s = prob.value
        # print("status:", prob.status)
        # print("optimal value", prob.value)

        if s > 0 or prob.status == cp.INFEASIBLE:
            # infeasible
            print("infeasible. prob.status=", prob.status,", s=", s)
            # find_X()
            return None, None, None

        # if prob.status == cp.OPTIMAL:
        #     return utils.pflat(x.value[0:-1])

        # problem is feasible


        return (*find_X(), prob.status)
        #########################


        # if x.value is None or x.value[-1] > 0:
        #     # there is no solution for this gamma
        #     return None

        # if x.value is None: # TODO: WHY??
        #     # print("x.value = None")
        #     return None, None
        #
        # if x.value[-1] > 0:
        #     # print("s > 0")
        #     return None, None
        #
        # solution = utils.pflat(x.value[0:-1])

#         ##### VERIFY THE SOLUTION
#         print("status:", prob.status)
#         v = x.value
#         eps = 1e-7
#         if v[0] < eps and v[1] < eps and v[2] < eps and v[3] < eps and v[4] < eps:
#             v[0] = 0
#             v[1] = 0
#             v[2] = 0
#             v[3] = 0
#             v[4] = 0
#
#         # X = v[0:-1]
#         # s = v[-1]  # prob.value
#         for i in range(n_cam):
#             # Ai = A[i, :, :]
#             # ci = c[i, :]
#             # print(f"||A{i} @ X||<=gamma c{i}.T @ X + s", f"{np.linalg.norm(Ai @ X)} <= {gamma * ci @ X + s}" )
#
#             Bi = B[i]
#             ai = a[i]
#             left = np.linalg.norm(Bi @ v)
#             right = ai.T @ v
#             print(f"||B{i} @ v||<= a{i}.T @ v", f"{left} <= {right}", "right - left=", right - left)
# # soc_constraints = [cp.SOC(a[i].T @ x, B[i] @ x) for i in range(n_cam)]
#         # We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
#
#         # print(v[0] <= 1e-7)
#         ##### VERIFY THE SOLUTION -- END

        # return solution, x.value[-1]

    except Exception as e:
        print("Exception:", str(e))
        return None, None, None
