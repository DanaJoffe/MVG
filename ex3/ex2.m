% Dana Joffe 312129240

load compEx1data.mat
load Ex1results.mat

% compute the camera matrices
e2 = null(F');
P = {[eye(3) [0 0 0]'], ...
    [[0, -e2(3), e2(2) ; e2(3), 0, -e2(1); -e2(2), e2(1), 0] * F e2]};

% normalize
x_norm = {N{1} * x{1}, N{2} * x{2}};
P_norm = {N{1} * P{1}, N{2} * P{2}};

% triangulate
n = length(x_norm{1});
X = zeros(4,n);
for i=1:n
    X(:,i) = triangulation(P_norm, [x_norm{1}(:,i) x_norm{2}(:,i)]);
end

% plot images, image points, and projected 3D points
im = {imread('kronan1.JPG'), imread('kronan2.JPG')};
proj = {pflat(P{1} * X), pflat(P{2} * X)};

figure
imagesc(im{1})
colormap gray
hold on
axis equal
plot(proj{1}(1, :), proj{1}(2, :),'yo');
hold on
plot(x{1}(1, :), x{1}(2 , :),'c.');
legend('projected points', 'measured points')
title("2.1 Triangulation - camera 1")

figure
imagesc(im{2})
colormap gray
hold on
axis equal
plot(proj{2}(1, :), proj{2}(2, :),'yo');
hold on
plot(x{2}(1, :), x{2}(2 , :),'c.');
legend('projected points', 'measured points')
title("2.2 Triangulation - camera 2")

figure
plot3(X(1 ,:) ,X(2 ,:) ,X(3 ,:) , '.')
grid on
title("2.3 Triangulation - 3D points")


function X=triangulation(P, x)
% the same function from homework 2 computer exercise 4.
% x (matrix) are different measures of the same 3D point and P (cells) are the camera
% matrices. Return X - the 3D point.
    n = length(P); % number of cameras
    M = zeros(3 * n, 4 + n);
    for i=1:n
        row = 3 * i - 2;
        M(row:row+2, 1:4) = P{i};
        col = 4 + i;
        M(row:row+2, col) = -x(:, i);
    end
    % solve using SVD & set up the camera matrix
    [~,~,V] = svd(M);
    v = V(:, end);
    
    X = v(1:4, 1);
    X = pflat(X);
end