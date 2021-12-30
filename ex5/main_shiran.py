import numpy as np
import utils_shiran as utils
import matplotlib.pyplot as plt
import triangulation_shiran as triangulation

data = np.load('Drinking Fountain Somewhere In Zurich.npz')
P = data['P']
X = data['X']
x = data['x']
visible_points = data['visible_points']
K = data['K']

# Computer Exercise 1

proj_error_gt = utils.reprojection_errors(P, X, x, visible_points)
mean_err_gt = proj_error_gt[visible_points].mean()
plt.hist(proj_error_gt[visible_points], bins=30)
plt.title("Reprojection Error - Ground Truth 3D Points")
plt.show()
print(mean_err_gt)


X_Sol_DLT = triangulation.DLT_triangulation(P, x, visible_points)
proj_error_dlt = utils.reprojection_errors(P, X_Sol_DLT, x, visible_points)
mean_err_dlt = proj_error_dlt[visible_points].mean()
print(mean_err_dlt)

norm_P, norm_x = utils.normalize_cam_points(P, x, np.linalg.inv(K))
X_Sol__DLT_norm = triangulation.DLT_triangulation(norm_P, norm_x, visible_points)
proj_error_dlt_norm = utils.reprojection_errors(P, X_Sol__DLT_norm, x, visible_points)
mean_err_dlt_norm = proj_error_dlt_norm[visible_points].mean()
print(mean_err_dlt_norm)


# Computer Exercise 2

# X_Sol_SOCP = triangulation.SOCP_triangulate(x, P, visible_points, max_iter=10)
# proj_error_SOCP = utils.reprojection_errors(P, X_Sol_SOCP, x, visible_points)
# mean_err_SOCP = proj_error_SOCP[visible_points].mean()
# print(mean_err_SOCP)



X_Sol_SOCP_norm = triangulation.SOCP_triangulate(norm_x, norm_P, visible_points, max_iter=20)#, high=0.0001, tol=1e-6)
proj_error_SOCP_norm = utils.reprojection_errors(P, X_Sol_SOCP_norm, x, visible_points)
mean_err_SOCP_norm = proj_error_SOCP_norm[visible_points].mean()
print(mean_err_SOCP_norm)

exit(0)

plt.xticks(range(5), ('DLT', 'DLT Normalized', 'SOCP', 'SOCP-Normalized', 'GT'), rotation=20)
bar = plt.bar(np.arange(5), [mean_err_dlt, mean_err_dlt_norm, mean_err_SOCP, mean_err_SOCP_norm, mean_err_gt])
for rect in bar:
    height = rect.get_height()
    plt.text(rect.get_x() + rect.get_width() / 2.0, height, f'{height:.4f}', ha='center', va='bottom')
plt.tight_layout()
plt.show()

utils.plot_cameras(P, K, X_Sol_DLT, "DLT Without Normalization - 3D Reconstruction")
utils.plot_cameras(P, K, X_Sol__DLT_norm, "DLT With Normalization - 3D Reconstruction")
utils.plot_cameras(P, K, X_Sol_SOCP, "SOCP Without Normalization - 3D Reconstruction")
utils.plot_cameras(P, K, X_Sol_SOCP_norm, "SOCP With Normalization - 3D Reconstruction")