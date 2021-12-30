% Dana Joffe 312129240
load cube_matches.mat
load Ex3results.mat

%% a
im = {imread('cube1.JPG'), imread('cube2.JPG')};
n = length(x1); % number of points to triangulate
x = {[x1; ones(1,n)], [x2; ones(1,n)]}; % x in homogenous coordinates

% triangulate
X = zeros(4,n);
for i=1:n
    X(:,i) = triangulation(P, [x{1}(:,i) x{2}(:,i)]);
end

proj = {};
proj{1} = pflat(P{1} * X);
proj{2} = pflat(P{2} * X);

figure
imagesc(im{1})
colormap gray
hold on
axis equal
plot(proj{1}(1, :), proj{1}(2, :),'ro');
hold on
plot(x{1}(1, :), x{1}(2 , :),'*');
legend('projected points', 'measured points')
title("4.a Triangulation - camera 1")

figure
i=2;
imagesc(im{i})
colormap gray
hold on
axis equal
plot(proj{i}(1, :), proj{i}(2, :),'ro');
hold on
plot(x{i}(1, :), x{i}(2 , :),'*');
legend('projected points', 'measured points')
title("4.a Triangulation - camera 2")


% Normalize points and cameras
x_norm = {inv(K{1}) * x{1}, inv(K{2}) * x{2}};
P_norm = {inv(K{1}) * P{1}, inv(K{2}) * P{2}};

% triangulate
Xn = zeros(4,n);
for i=1:n
    Xn(:,i) = triangulation(P_norm, [x_norm{1}(:,i) x_norm{2}(:,i)]);
end

projn = {};
projn{1} = pflat(P{1} * Xn);
projn{2} = pflat(P{2} * Xn);

figure
imagesc(im{1})
colormap gray
hold on
axis equal
plot(projn{1}(1, :), projn{1}(2, :),'ro');
hold on
plot(x{1}(1, :), x{1}(2 , :),'*');
legend('projected points', 'measured points')
title("4.a Triangulation with normalization - camera 1")

figure
i=2;
imagesc(im{i})
colormap gray
hold on
axis equal
plot(projn{i}(1, :), projn{i}(2, :),'ro');
hold on
plot(x{i}(1, :), x{i}(2 , :),'*');
legend('projected points', 'measured points')
title("4.a Triangulation with normalization - camera 2")

%% b
goodpoints = sqrt(sum((x1 - projn{1}(1:2 ,:)).^2)) < 3 & sqrt(sum((x2 - projn{2}(1:2 ,:)).^2)) < 3;
Xn_rm = Xn(:, goodpoints);

figure
plot3(Xn_rm(1 ,:) ,Xn_rm(2 ,:) ,Xn_rm(3 ,:) , '.','Markersize' ,2);
hold on
plotcams(P)
hold on
plot3( [ Xmodel(1, startind ); Xmodel(1, endind )],...
    [ Xmodel(2, startind ); Xmodel(2, endind )],...
    [ Xmodel(3, startind ); Xmodel(3, endind )], 'b-');

axis equal
title("4.b Remaining 3D points, model & cameras")


function X=triangulation(P, x)
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
        
    % find the correct solution s.t. the points are in front of the cameras
    proj = P{1}*X;
    X = X * sign(proj(3,1));
    X = pflat(X);
end