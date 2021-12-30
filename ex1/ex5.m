%% Dana Joffe 312129240
%% Computer Exercise 5a
im = imread('compEx5.jpg');
load compEx5.mat

figure
imagesc(im)
colormap gray
axis equal
hold on
plot(corners(1 ,[1: end 1]) , corners(2 ,[1: end 1]), 'y.-');
title("5.1 corners on top of image1")

P1 = K * [1 0 0 0; 0 1 0 0; 0 0 1 0];

% normalize the corner points
n_corners = pflat(inv(K) * corners);
figure
plot(n_corners(1 ,[1: end 1]) , n_corners(2 ,[1: end 1]), '.-');
title("5.2 normalized corner points")
axis equal
axis ij

kalib_P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];

% compute the camera center and principal axis of the camera.
c_center = pflat(null(kalib_P1)); c_center  = c_center(1:3, :);
p_axis = det(kalib_P1(:, 1:3)) * kalib_P1(3, 1:3)'; p_axis = pflat(p_axis);

set(0,'DefaultFigureVisible','on');
figure
plot3(c_center(1), c_center(2), c_center(3), 'r+')
title("5.3 3D world with 2 cameras")
hold on
quiver3(c_center(1) ,c_center(2) ,c_center(3) ,p_axis(1) ,p_axis(2) ,p_axis(3), 1, 'r')
hold on

% compute the 3D points in the plane v that project onto the corner points
plane_pi = pflat(v);
pi_ = plane_pi(1:3); 
s = -pi_' * n_corners;
U = [n_corners; s]; U = pflat(U);

plot3( U(1 ,[1: end 1]) , U(2 ,[1: end 1]),U(3 ,[1: end 1]) , 'b.-');
hold on
plot3( n_corners(1 ,[1: end 1]) , n_corners(2 ,[1: end 1]),n_corners(3 ,[1: end 1]) , 'r.-');
hold on

%% Computer Exercise 5b
% plot camera2
c_center2 = [2 0 0];
R = [sqrt(3)/2 0 1/2 ; 0 1 0; -1/2 0 sqrt(3)/2];
t = c_center2';
P2 = R * [eye(3) -t];
p_axis2 = det(P2(:, 1:3)) * P2(3, 1:3)'; p_axis2 = pflat(p_axis2);

plot3(c_center2(1), c_center2(2), c_center2(3), 'g+')
hold on
quiver3(c_center2(1) ,c_center2(2) ,c_center2(3) ,p_axis2(1) ,p_axis2(2) ,p_axis2(3), 1, 'g')
hold on
axis equal
axis ij

% project the 3D corner points into the image using the camera matrix
p_corners = P2 * U; p_corners = pflat(p_corners);
pCorWorld = R' * p_corners + t;
plot3(pCorWorld(1 ,[1: end 1]) , pCorWorld(2 ,[1: end 1]),pCorWorld(3 ,[1: end 1]) , 'g.-');

legend('camera1 center', 'camera1 principal axis', 'points in plane v', 'image1 plane', 'camera2 center', 'camera2 principal axis', 'image2 plane')
xlabel('X') 
ylabel('Y')
zlabel('Z')

% transform the normalized corner points to the new image
H = R - t * pi_';
t_corners = H * n_corners; t_corners = pflat(t_corners);

figure
plot(t_corners(1 ,[1: end 1]) , t_corners(2 ,[1: end 1]), '.-');
title("5.4 transformed corner points using homography")
axis equal
axis ij

figure
plot(p_corners(1 ,[1: end 1]) , p_corners(2 ,[1: end 1]), '.-');
title("5.5 projected corner points using camera matrix")
axis equal
axis ij

% Transform the original image and corner points using normalized homography
H_tot = K * H * inv(K);
tform = maketform ('projective', H_tot');
[new_im ,xdata , ydata] = imtransform(im ,tform, 'size', size(im)); 
new_corners = pflat(H_tot * corners);

figure
imagesc(xdata ,ydata , new_im); 
colormap gray
hold on
plot(new_corners(1 ,[1: end 1]) , new_corners(2 ,[1: end 1]), 'y.-');
title("5.6 transformed image & corners using normalized homography")
axis equal
axis ij
