import numpy as np
import cvxpy as cp
import utils
import warnings
warnings.filterwarnings("ignore")


def DLT_triangulation(P, x, visible_points):
    """
    Use direct linear transformation in order to triangulate the points.
    :param P: ndarray of shape [n_cam, 3, 4], the cameras
    :param x: ndarray of shape [n_cam, 3, n_points], the projected image points
    :param visible_points: boolean matrix of shape [n_cam, n_points], what cameras see what points
    :return: X ndarray of shape [4, n_points], the predicted 3D points
    """
    def triangulation(P, x):
        """
        Find a single 3D point.

        :param P: ndarray of shape [n_see_cam, 3, 4], the cameras that see the point
        :param x: ndarray of shape [n_see_cam, 3], the projected image point on all the cameras that see.
        :return: X ndarray of shape [4, ], the predicted 3D point
        """
        # build the matrix M
        x = np.transpose(x)
        n = P.shape[0]
        M = np.zeros((3 * n, 4 + n))
        for i in range(n):
            row = 3 * i
            M[row: row + 3, 0:4] = P[i, :, :]
            col = 4 + i
            M[row: row + 3, col] = -x[:, i]

        # solve using SVD & set up the camera matrix
        U, S, Vt = np.linalg.svd(M)
        V = np.transpose(Vt)
        v = V[:, -1]
        X = v[0:4]
        return utils.pflat(X)

    n_cam, _, n_points = x.shape
    X = np.zeros((4, n_points))

    # loop over the number of 3D points
    for i in range(n_points):
        # take the cameras that see the current point
        P_see = P[visible_points[:, i], :, :]
        x_see = x[visible_points[:, i], :, i]
        X[:, i] = triangulation(P_see, x_see)
    return X


def get_triangulation_parameters(P, x):
    """
    Given an array of cameras and projected points, find the parameters needed for SOCP.
    :param P: cameras array of size [n_cam, 3, 4].
    :param x: projected points of size [n_cam, 3, n_points]
    :return: A,c - the variables for SOCP.
    A ndarray of shape [n_cam, n_points, 3, 4]
    c ndarray of shape [n_cam, 4]

    Each variable should be a list in the length of the number of constraints (number of cameras that see the point).
    The i'th constraint will correspond to: ||A[i] * x|| <= gamma* c[i]^T * x
    """
    # A_i = x_i outer product P_i^3 - P_i
    P3 = P[:, -1, :]  # shape [n_cam, 4]
    A = np.einsum('ipj,ic->ijpc', x, P3) - np.expand_dims(P, axis=1)  # shape [n_cam, n_points 3, 4]

    # c_i = P_i^3 transpose
    c = P3
    return A, c


def solve_SOCP(A,c, gamma):
    """
    Solve SOC feasibility problem FOR A SINGLE POINT todo.
    Use cvxpy to solve a second order cone program of the form:
    minimize f^Tx s.t ||B[i]x + b[i]|| <= a[i]^Tx+d[i] for i=1,...
    :param A: a list of numpy arrays,    for a single point, shape [n_see_cam, d, 4] todo
    :param c: a list of numpy arrays                         shape [n_see_cam, 4] todo
    :return: If there is a solution return it. else, return None.
     A solution is an array of size [4,1]. Make sure it's least coordinate is 1.

    Hint: use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
    """
    try:
        # objective: find X s.t. for all i: ||Ai @ X|| <= gamma ci.T @ X
        x = cp.Variable((4,))
        n_cam = A.shape[0]
        # We use cp.SOC(t, x) to create the SOC constraint ||x||_2 <= t.
        soc_constraints = [cp.SOC(gamma * c[i, :].T @ x, A[i, :, :] @ x) for i in range(n_cam)]
        obj = cp.Minimize(cp.Constant((1,)))  # don't minimize anything
        prob = cp.Problem(obj, soc_constraints)
        prob.solve()
        X = x.value
        return utils.pflat(X)
    except:
        return None


def SOCP_triangulate(x, P, max_iter, low=0, high=1, tol=1e-2):
    """
    Get a single 3D point from a list of cameras that see it and corresponding image points.
    Use an iterative algorithm similar to bisection.
    :param x: ndarray of shape [n_i, 3] the matching image points for the 3D points we are looking for
    :param P: ndarray of shape [n_i, 3, 4] cameras that see the 3D point
    :param max_iter: maximal number of iterations before stopping.
    :param low: minimal current value for gamma.
    :param high: maximal current value for gamma.
    :param tol: the boundaries around gamma which are close enough for a good solution.
    :return: Xj ndarray of shape [4,1] the point we got in the last feasible iteration
    """
    x = np.expand_dims(x, axis=2)  # shape [n_i, 3, 1]
    A, c = get_triangulation_parameters(P, x)
    A = np.squeeze(A)  # reduce dimension

    # algorithm in the exercise
    X = None
    curr_iter = 0
    gamma = high
    while high - low > tol and curr_iter < max_iter:
        curr_iter += 1
        Xj = solve_SOCP(A, c, gamma)
        if Xj is not None:
            max_error = utils.max_single_point_errors(P, x, Xj)
        if Xj is not None and max_error <= gamma:
            high = max_error
            X = Xj
        else:
            low = gamma
            high = 2 * high
        gamma = (high + low) / 2
    return X

