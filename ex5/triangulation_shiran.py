import numpy as np
import cvxpy as cp
import utils_shiran as utils
from cvxpy import error

def DLT_triangulation(P, x, visible_points):
    """
    Use direct linear transformation in order to triangulate the points.
    :param P: ndarray of shape [n_cam, 3, 4], the cameras
    :param x: ndarray of shape [n_cam, 3, n_points], the projected image points
    :param visible_points: boolean matrix of shape [n_cam, n_points], what cameras see what points
    :return: X ndarray of shape [4, n_points], the predicted 3D points
    """

    num_points_3D_points = x.shape[2]
    X_Sol = np.zeros((4, num_points_3D_points))
    for i in range(num_points_3D_points):
        points = x[:, :, i][visible_points[:, i]]
        Cameras = P[visible_points[:, i]]
        num_points = points.shape[0]
        M = np.zeros((3 * num_points, 4 + num_points))
        k = 0
        for j in range(num_points):
            M[k: k + 3, 0: 4] = Cameras[j]
            M[k, 4 + j] = -points[j, 0]
            M[k + 1, 4 + j] = -points[j, 1]
            M[k + 2, 4 + j] = -1
            k = k + 3
        U, S, V = np.linalg.svd(M)
        X = utils.pflat(V[-1, :4])
        X_Sol[:, i] = X
    return X_Sol


def get_triangulation_parameters(P, x):
    """
    Given an array of cameras and projected points, find the parameters needed for SOCP.
    :param P: cameras array of size [n_cam, 3, 4].
    :param x: projected points of size [n_cam, 2]
    :return: A,c - the variables for SOCP.
    Each variable should be a list in the length of the number of constraints (number of cameras that see the point).
    The i'th constraint will correspond to: ||A[i] * x|| <= gamma* c[i]^T * x
    """
    num_cameras = P.shape[0]
    A = np.zeros((2, 4, num_cameras))
    c = np.zeros((4, num_cameras))
    for j in range(num_cameras):
        A[0, :, j] = x[j, 0]*P[j, 2, :] - P[j, 0, :]
        A[1, :, j] = x[j, 1] * P[j, 2, :] - P[j, 1, :]
        c[:, j] = P[j, 2, :]
    return A, c


def solve_SOCP(A,c, gamma):
    """
    Use cvxpy to solve a second order cone program of the form:
    minimize f^TX s.t ||A[i]x + b[i]|| <= c[i]^Tx+d[i] for i=1,...
    :param A: a list of numpy arrays
    :param c: a list of numpy arrays
    :return: If there is a solution return it. else, return None.
     A solution is an array of size [4,1]. Make sure it's least coordinate is 1.

    Hint: use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
    """
    num_constraints = A.shape[2]
    X = cp.Variable((4,))
    soc_constraints = [
        cp.SOC(gamma*c[:, i].T @ X, A[:, :, i] @ X) for i in range(num_constraints)
    ]
    problem = cp.Problem(cp.Minimize(cp.expressions.constants.Constant(2)), soc_constraints)
    try:
        problem.solve()
    except error.SolverError:
        return None

    # ##### VERIFY THE SOLUTION
    # x = utils.pflat(X.value) #X.value
    # for i in range(num_constraints):
    #     Ai = A[:, :, i]
    #     ci = c[:, i]
    #     left = np.linalg.norm(Ai @ x)
    #     right = gamma * ci.T @ x
    #     print(f"||A{i} @ X||<=gamma * c{i}.T @ X", f"{left} <= {right}", "right - left=", right - left)

    return utils.pflat(X.value)


def SOCP_triangulate_Single_Point(x, P, max_iter, low=0, high=1, tol=1e-2):
    """
    Get a single 3D point from a list of cameras that see it and corresponding image points.
    Use an iterative algorithm similar to bisection.
    :param x: ndarray of shape [n_i, 3] the matching image points for the 3D points we are looking for
    :param P: ndarray of shape [n_i, 3, 4] cameras that see the 3D point
    :param max_iter: maximal number of iterations before stopping.
    :param low: minimal current value for gamma.
    :param high: maximal current value for gamma.
    :param tol: the boundaries around gamma which are close enough for a good solution.
    :return: Xi ndarray of shape [4,1] the point we got in the last feasible iteration
    """

    gamma = high
    best_U = None
    # x = np.concatenate([x, np.ones((x.shape[0], 1))], axis=1)
    A, c = get_triangulation_parameters(P, x)
    curr_iter = 0
    while high - low > tol and curr_iter < max_iter:
        U = solve_SOCP(A, c, gamma)
        if U is not None:
            max_error_U = utils.max_single_point_errors(P, x, Xi=U)
        if U is not None and max_error_U < gamma:
            high = max_error_U
            best_U = U
        else:
            low = gamma
            high = 2*high
        gamma = (high+low)/2
        curr_iter+=1
    if best_U is None:
        raise Exception("No feasible solution found")
    return best_U


def SOCP_triangulate(x,P, visible_points, max_iter, low=0, high=1, tol=1e-2):
    num_points = x.shape[2]
    X_Sol = np.zeros((4, num_points))
    for i in range(num_points):
        Pi = P[visible_points[:, i]]
        xi = x[visible_points[:, i], :, i]
        Ui = SOCP_triangulate_Single_Point(xi, Pi, max_iter=max_iter, low=low, high=high, tol=tol)
        X_Sol[:, i] = Ui

    return X_Sol