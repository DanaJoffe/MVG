%% Dana Joffe 312129240
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