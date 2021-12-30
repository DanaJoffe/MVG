from utils import *
from triangulation import *
import matplotlib.pyplot as plt


data = np.load('Drinking Fountain Somewhere In Zurich.npz')
P = data['P']  # shape [n_cam, 3, 4]
X = data['X']  # shape [4, n_points]
x = data['x']  # shape [n_cam, 3, n_points]
visible_points = data['visible_points']  # shape [n_cam, n_points]
K = data['K']

# add ones to x
ones = np.ones((x.shape[0], 1, x.shape[2]))
x = np.concatenate((x, ones), axis=1)

plot = {'histogram': True,  # question 1b
        'MREs': True,  # question 1e
        'bar': True  # question 2e
        }

# triangulation methods
DLT_RAW = 'DLT triangulation un-normalized'
DLT_nN = 'DLT triangulation normalized with N'
DLT_nK = 'DLT triangulation normalized with inv(K)'
GT = 'ground truth 3D points'
SOC_RAW = 'SOCP triangulation un-normalized'
SOC_nN = 'SOCP triangulation normalized with N'
SOC_nK = 'SOCP triangulation normalized with inv(K)'

# save all the mean reprojected errors, for various methods
MRE = {DLT_RAW: None,
       DLT_nN: None,
       DLT_nK: None,
       GT: None,
       SOC_RAW: None,
       SOC_nN: None,
       SOC_nK: None}

# save all the triangulated 3D points, for various methods
points_3D = MRE.copy()

n_points = x.shape[2]


# ========================> 1b - histogram <========================
re_gt = reprojection_errors(P, X, x, visible_points)  # ground truth
MRE[GT] = np.mean(re_gt[visible_points])
if plot['histogram']:
    fig, ax = plt.subplots(tight_layout=True)
    ax.hist(re_gt[visible_points].flatten(), 100)
    plt.title("1.b Reprojection Errors Histogram")
    plt.xlabel("reprojection error")
    plt.ylabel("# of points")
    plt.show()


# ========================> 1e - triangulate <========================
# normalization matrix
def get_normalization_mat(x):
    m = np.mean(x, 2)
    s = np.std(x, 2)
    return np.array([[[1/s[i,0], 0, -m[i,0]/s[i,0]],
                      [0, 1/s[i,1], -m[i,1]/s[i,1]],
                      [0, 0, 1]]
                     for i in range(x.shape[0])])

# Triangulate: Unnormalized cameras and points
points_3D[DLT_RAW] = DLT_triangulation(P, x, visible_points)
MRE[DLT_RAW] = np.mean(reprojection_errors(P, points_3D[DLT_RAW], x, visible_points)[visible_points])

# Triangulate: Cameras and points are normalized with K^−1
norm_Pk, norm_xk = normalize_cam_points(P, x, np.linalg.inv(K))
points_3D[DLT_nK] = DLT_triangulation(norm_Pk, norm_xk, visible_points)
MRE[DLT_nK] = np.mean(reprojection_errors(P, points_3D[DLT_nK], x, visible_points)[visible_points])

# Calc mean reprojection error after triangulation, with normalization by N matrix (zero mean, std 1)
norm_Pn, norm_xn = normalize_cam_points(P, x, get_normalization_mat(x))
points_3D[DLT_nN] = DLT_triangulation(norm_Pn, norm_xn, visible_points)
MRE[DLT_nN] = np.mean(reprojection_errors(P, points_3D[DLT_nN], x, visible_points)[visible_points])

# Plot a comparison table
def show_table(title):
    column_headers = ['Ground Truth', 'Unnormalized', 'Normalized with N', 'Normalized with K']
    row_headers = ['MRE']
    mres = [MRE[GT], MRE[DLT_RAW], MRE[DLT_nN], MRE[DLT_nK]]
    data = [[f'{mre:.4f}' for mre in mres]]

    fig, ax = plt.subplots()
    ax.set_axis_off()
    the_table = ax.table(cellText=data,
                         rowLabels=row_headers,
                         rowColours=["palegreen"] * len(row_headers),
                         colColours=["palegreen"] * len(column_headers),
                         colLabels=column_headers,
                         cellLoc='center',
                         loc='center')
    the_table.set_fontsize(40)
    ax.set_title(title,
                 fontweight="bold")
    plt.show()

if plot['MREs']:
    show_table('1.e Mean Reprojection Errors (MREs) with DLT Triangulation')


# ========================> 2e - Triangulate all the points using ’SOCP_triangulate’ <========================
def SOC_triangulation(x, P, **kwargs):
    X_socp = np.zeros(X.shape)
    for j in range(n_points):
        X_socp[:, j] = SOCP_triangulate(x[visible_points[:, j], :, j], P[visible_points[:, j], :, :], **kwargs)
        # report progress every 100 points
        if j % 100 == 0:
            print(f"triangulate point {j}/{n_points}...")
    return X_socp

# Triangulate all the points using ’SOCP_triangulate’ without normalization
points_3D[SOC_RAW] = SOC_triangulation(x, P, max_iter=30)
MRE[SOC_RAW] = np.mean(reprojection_errors(P, points_3D[SOC_RAW], x, visible_points)[visible_points])

# Triangulate all the points using ’SOCP_triangulate’ with normalization
points_3D[SOC_nK] = SOC_triangulation(norm_xk, norm_Pk, max_iter=30)  # , tol=1e-6, high=1e-1)
MRE[SOC_nK] = np.mean(reprojection_errors(P, points_3D[SOC_nK], x, visible_points)[visible_points])

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
    plt.title(title)
    # fig.autofmt_xdate()

    for index, data in enumerate(mres):
        ax.text(x=index, y=data, s=f'{data:.4f}', ha='center', va='bottom')
    plt.show()

if plot['bar']:
    show_bar("2.e Mean Reprojection Error (MRE) after Triangulation")


# ========================> 2f - plot <========================
for name, points3D in points_3D.items():
    if points3D is not None:
        utils.plot_cameras(P, K, points3D, name + ' - reconstruction')
