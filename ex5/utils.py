import numpy as np
import plotly
import plotly.graph_objects as go


def max_single_point_errors(P, x, Xi):
    """
    Get the maximum reprojection error of the 3D predicted point U with the cameras P and real image points u.
    this is calculated on all the given cameras (and not checking which cameras see which point).
    :param P: the cameras, with shape [n_cam, 3, 4]
    :param x: image points of size [n_cam, 3]
    :param Xi: 3D point of size [4, 1]
    :return: The maximum reprojection error over all the cameras. (a scalar)
    The reprojection error on camera i is:
    ei = sqrt(sum((u_i - P_i*U)**2))
    """
    if np.ndim(Xi) == 1:
        Xi = np.expand_dims(Xi, axis=1)
    if np.ndim(x) == 2:
        x = np.expand_dims(x, axis=2)
    re = reprojection_errors(P, Xi, x)
    return np.max(re)


def normalize_cam_points(P,x,N):
    """
    Normalize the camera matrices and the image points with normalization matrices N.
    :param P: ndarray of shape [n_cam, 3, 4], the cameras
    :param x: ndarray of shape [n_cam, 3, n_points], the projected image points
    :param N: ndarray of shape [n_cam, 3, 3], the normalization matrices
    :return: norm_P: ndarray of shape [n_cam, 3, 4], the normalized cameras
            norm_x: ndarray of shape [n_cam, 3, n_points], the normalized image points
    """
    norm_P = np.einsum('irc,ick->irk', N, P)  # N * P
    norm_x = np.einsum('ikj,ijp->ikp', N, x)  # N * x
    return norm_P, norm_x


def reprojection_errors(P, X, x, visible_points=None):
    """
    Projects the 3D points in X to the cameras P and computes the distance to the real image points x.
    :param P: ndarray of shape [n_cam, 3, 4], the cameras
    :param X: ndarray of shape [4, n_points], the predicted 3D points
    :param x: ndarray of shape [n_cam, 3, n_points], the projected image points
    :param visible_points: boolean matrix of shape [n_cam, n_points], what cameras see what points
    :return: errors: ndarray of shape [n_cam, n_points], in the ij entry has ||x_ij - pflat(P_i*X_j)||.
    The errors in the non-visible entries should be np.nan
    """
    projections = np.einsum('ikj,jp->ikp', P, X)  # P * X
    projections = pflat(projections)  # shape [n_cam, 3, n_points]

    # calculate the norm error
    errors = np.sqrt(np.sum((x - projections) ** 2, axis=1))

    if visible_points is not None:
        errors[~visible_points] = np.nan
    return errors


def decompose_camera_matrix(P, K):
    """
    Decompose camera matrices to R and t s.t P[i] = K*R^T[I -t]
    :param P: ndarray of shape [n_cam, 3, 4], the cameras
    :param K: ndarray of shape [n_cam, 3, 3], the calibration matrices
    :return: R: ndarray of shape [n_cam, 3, 3]
            t: ndarray of shape [n_cam, 3]
    """
    Rt = np.linalg.inv(K) @ P
    Rs = np.transpose(Rt[:, :, :3],(0,2,1))
    ts = np.squeeze(-Rs @ np.expand_dims(Rt[:, 0:3, 3], axis=-1))
    return Rs, ts


def pflat(x):
    """
    :param x: ndarray of shape [n_cam, dim, n_points] OR [dim, n_points] OR [dim, ]
    :return: x in homogeneous coordinates.
    """
    if np.ndim(x) == 1:
        return x/x[-1]
    elif np.ndim(x) == 2:
        return x / np.expand_dims(x[-1, :], axis=0)
    return x / np.expand_dims(x[:, -1, :], axis=1)


def plot_cameras(P, K, X, title='reconstruction'):
    """
    Plot a 3D image of the points and cameras
    :param P: ndarray of shape [n_cam, 3, 4], the cameras
    :param K: ndarray of shape [n_cam, 3, 3], the calibration matrices
    :param X: ndarray of shape [4, n_points], the predicted 3D points
    :param title: the name of the plot
    """
    R,t = decompose_camera_matrix(P, K)
    data = []
    data.append(get_3D_quiver_trace(t, R[:, :3, 2], color='#86CE00', name='cam_learn'))
    data.append(get_3D_scater_trace(t.T, color='#86CE00', name='cam_learn', size=1))
    data.append(get_3D_scater_trace(X[:3,:], '#3366CC', '3D points', size=0.5))

    fig = go.Figure(data=data)
    path = title+'.html'
    plotly.offline.plot(fig, filename=path, auto_open=False)


def get_3D_quiver_trace(points, directions, color='#bd1540', name=''):
    assert points.shape[1] == 3, "3d cone plot input points are not correctely shaped "
    assert len(points.shape) == 2, "3d cone plot input points are not correctely shaped "
    assert directions.shape[1] == 3, "3d cone plot input directions are not correctely shaped "
    assert len(directions.shape) == 2, "3d cone plot input directions are not correctely shaped "

    trace = go.Cone(
        name=name,
        x=points[:, 0],
        y=points[:, 1],
        z=points[:, 2],
        u=directions[:, 0],
        v=directions[:, 1],
        w=directions[:, 2],
        sizemode='absolute',
        sizeref=0.5,
        showscale=False,
        colorscale=[[0, color], [1, color]],
        anchor="tail"
    )

    return trace


def get_3D_scater_trace(points, color, name,size=0.5):
    assert points.shape[0] == 3, "3d plot input points are not correctely shaped "
    assert len(points.shape) == 2, "3d plot input points are not correctely shaped "

    trace = go.Scatter3d(
        name=name,
        x=points[0, :],
        y=points[1, :],
        z=points[2, :],
        mode='markers',
        marker=dict(
            size=size,
            color=color,
        )
    )

    return trace
