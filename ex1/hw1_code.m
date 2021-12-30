%% Dana Joffe 312129240
%% Computer Exercise 1
load compEx1.mat

a = pflat(x2D);
figure
plot(a(1 ,:) ,a(2 ,:) , '.')
title("1.1 pflat on 2D points")
axis equal

b = pflat(x3D);
figure
plot3(b(1 ,:) ,b(2 ,:) ,b(3 ,:) , '.')
title("1.2 pflat on 3D points")
axis equal

%% Computer Exercise 2
im = imread('compEx2.JPG ');
figure
imagesc(im)
colormap gray
title("2. intersection of parallel lines")
hold on

load compEx2.mat
plot(p1(1 ,:) ,p1(2 ,:) , 'r.', 'MarkerSize', 10, 'LineWidth', 2)
plot(p2(1 ,:) ,p2(2 ,:) , 'y.', 'MarkerSize', 10, 'LineWidth', 2)
plot(p3(1 ,:) ,p3(2 ,:) , 'g.', 'MarkerSize', 10, 'LineWidth', 2)

% compute the lines going through the points
l1 = cross(p1(:,1), p1(:,2));
rital(l1, 'r')
l2 = cross(p2(:,1), p2(:,2));
rital(l2, 'y')
l3 = cross(p3(:,1), p3(:,2));
rital(l3, 'g')

% Q: Do these lines appear to be parallel (in 3D)?
% A: Yes, but after being projected on the image plane - they intersect.

inter23 = pflat(cross(l2, l3));
plot(inter23(1 ,:) ,inter23(2 ,:) , 'b.', 'MarkerSize', 10, 'LineWidth', 2)

% Compute the distance between the first line and the the intersection point
d = abs(l1'*inter23) / sqrt(l1(1)^2 + l1(2)^2)

%% Computer Exercise 3
load compEx3.mat

figure
plot([startpoints(1 ,:); endpoints(1 ,:)], [startpoints(2 ,:); endpoints(2 ,:)], 'b-');
title("3.1 original grid")
axis equal
 
H1 = [sqrt(3) -1 1 ; 1 sqrt(3) 1 ; 0 0 2];
H2 = [1 -1 1 ; 1 1 0 ; 0 0 1];
H3 = [1 1 0 ; 0 2 0 ; 0 0 1];
H4 = [sqrt(3) -1 1 ; 1 sqrt(3) 1 ; 1/4 1/2 2];


homo_endpoints = [endpoints; ones(1, length(endpoints))];
homo_startpoints = [startpoints; ones(1, length(startpoints))];

% H1
trans_end = pflat(H1*homo_endpoints);
trans_start = pflat(H1*homo_startpoints);

figure
subplot(2, 2, 1)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H1 - Euclidean")
axis equal

% H2
trans_end = pflat(H2*homo_endpoints);
trans_start = pflat(H2*homo_startpoints);

% figure
subplot(2, 2, 2)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H2 - Similarity")
axis equal

% H3
trans_end = pflat(H3*homo_endpoints);
trans_start = pflat(H3*homo_startpoints);

% figure
subplot(2, 2, 3)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H3 - Affine")
axis equal

% H4
trans_end = pflat(H4*homo_endpoints);
trans_start = pflat(H4*homo_startpoints);

% figure
subplot(2, 2, 4)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H4 - Projective")
axis equal
suptitle('3.2 Transformations') 

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

