%% Dana Joffe 312129240
%% Computer Exercise 4
im1 = imread('compEx4im1.jpg');
im2 = imread('compEx4im2.jpg');
load compEx4.mat

u = pflat(U);

figure
imagesc(im1)
title("4.1 Points in U projected into camera P1")
colormap gray
hold on

% project the points in U into the cameras P1
u_proj_1 = pflat(P1 * u);
plot(u_proj_1(1 ,:) ,u_proj_1(2 ,:), '.', 'Markersize', 2)

figure
imagesc(im2)
colormap gray
title("4.2 Points in U projected into camera P2")
hold on

% project the points in U into the cameras P2
u_proj_2 = pflat(P2 * u);
plot(u_proj_2(1 ,:) ,u_proj_2(2 ,:), '.', 'Markersize', 2)

% Compute the camera centers and principal axes of the cameras.
c_center1 = pflat(null(P1)); c_center1  = c_center1(1:3, :);
c_center2 = pflat(null(P2)); c_center2 = c_center2(1:3, :);

p_axis1 = det(P1(:, 1:3)) * P1(3, 1:3)'; p_axis1 = pflat(p_axis1);
p_axis2 = det(P2(:, 1:3)) * P2(3, 1:3)'; p_axis2 = pflat(p_axis2);

figure
h = zeros(1,5);
h(1) = plot3(u(1 ,:) ,u(2 ,:) ,u(3 ,:) , '.', 'Markersize', 2);
title("4.3 3D points, camera centers & principal axes")
axis equal
hold on

% plot camera centers
h(2) = plot3(c_center1(1), c_center1(2), c_center1(3), 'r+');
hold on
h(3) = plot3(c_center2(1), c_center2(2), c_center2(3), 'g+');
hold on

% plot principal axes
s = 1;
h(4) = quiver3(c_center1(1) ,c_center1(2) ,c_center1(3) ,p_axis1(1) ,p_axis1(2) ,p_axis1(3), s, 'r');
hold on
h(5) = quiver3(c_center2(1) ,c_center2(2) ,c_center2(3) ,p_axis2(1) ,p_axis2(2) ,p_axis2(3), 8*s, 'g');
hold on
legend(h(2:end),'camera1 center','camera2 center','camera1 principal axis', 'camera2 principal axis')
