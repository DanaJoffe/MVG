% Dana Joffe 312129240
load compEx3data.mat

% construct normaization matrices
N = cell(1,2);
mean1 = mean(x{1}(1:2 ,:) ,2); % x and y of camera 1
std1 = std(x{1}(1:2 ,:) ,0 ,2); % x and y of camera 1
N{1} = [1/std1(1), 0, -mean1(1)/std1(1); 0, 1/std1(2), -mean1(2)/std1(2); 0 0 1];

mean2 = mean(x{2}(1:2 ,:) ,2);
std2 = std(x{2}(1:2 ,:) ,0 ,2);
N{2} = [1/std2(1), 0, -mean2(1)/std2(1); 0, 1/std2(2), -mean2(2)/std2(2); 0 0 1];

x_norm = cell(size(x));
x_norm{1} = N{1} * x{1}; x_norm{2} = N{2} * x{2}; 

figure
subplot(1, 2, 1)
plot(x_norm{1}(1 ,:) ,x_norm{1}(2 ,:) , '.');
title("camera 1")
subplot(1, 2, 2)
plot(x_norm{2}(1 ,:) ,x_norm{2}(2 ,:) , '.');
title("camera 2")
suptitle('3. Normalized measured projections')

%% a
Xmodel_hom = [Xmodel; ones(1, length(Xmodel))]; % to homogeneous coordinates
P = {};  n = {}; s = {};
[P{1}, n{1}, s{1}]=resectioning(x{1}, Xmodel_hom, N{1});
[P{2}, n{2}, s{2}]=resectioning(x{2}, Xmodel_hom, N{2});

%% b
proj_pt = {};
% Project model points into camera 1
im = imread('cube1.JPG');
proj_pt{1} = pflat(P{1} * Xmodel_hom);
visible = isfinite(x{1}(1 ,:));

figure
imagesc(im)
colormap gray
hold on
axis equal

title("3.b camera 1")
plot(proj_pt{1}(1, visible ), proj_pt{1}(2, visible ),'ro'); % Plots a red ’o’ at each visible point in xproj
hold on
plot(x{1}(1, visible ), x{1}(2 , visible ),'*'); % Plots a ’*’ at each point coordinate
legend('projected model points', 'measured image points')

% Project model points into camera 2
im = imread('cube2.JPG');
proj_pt{2} = pflat(P{2} * Xmodel_hom);
visible = isfinite(x{2}(1 ,:));

figure
imagesc(im)
colormap gray
hold on
axis equal

title("3.b camera 2")
plot(proj_pt{2}(1, visible ), proj_pt{2}(2, visible ),'ro'); % Plots a red ’o’ at each visible point in xproj
hold on
plot(x{2}(1, visible ), x{2}(2 , visible ),'*'); % Plots a ’*’ at each point coordinate
legend('projected model points', 'measured image points')

% plot Xmodel in 3D with lines
figure
plot3( [ Xmodel(1, startind ); Xmodel(1, endind )],...
    [ Xmodel(2, startind ); Xmodel(2, endind )],...
    [ Xmodel(3, startind ); Xmodel(3, endind )], 'b-');
title("3.b 3D model points & cameras")
axis equal
hold on
plotcams(P)
grid on

% Compute the inner parameters of the cameras
K = {};
[k, ~] = rq(P{1}(:,1:3)); K{1} = k./k(3,3);
[k, ~] = rq(P{2}(:,1:3)); K{2} = k./k(3,3);

%% c
RMSE = {sqrt(mean( sum( (x{1}(1:2, :) - proj_pt{1}(1:2, :) ).^2, 1) )), ...
        sqrt(mean( sum( (x{2}(1:2, :) - proj_pt{2}(1:2, :) ).^2, 1) ))};

% repeat the resectioning but without normalization
RMSE_ = {};
[P_, ~, ~]=resectioning(x{1}, Xmodel_hom, eye(3));
proj_pt_ = pflat(P_ * Xmodel_hom);
RMSE_{1} = sqrt(mean( sum( (x{1}(1:2, :) - proj_pt_(1:2, :) ).^2, 1) ));
   
[P_, ~, ~]=resectioning(x{2}, Xmodel_hom, eye(3));
proj_pt_ = pflat(P_ * Xmodel_hom);
RMSE_{2} = sqrt(mean( sum( (x{2}(1:2, :) - proj_pt_(1:2, :) ).^2, 1) ));

%% save cameras for Exercise 4
% save('Ex3results.mat','P', 'Xmodel', 'K', 'startind', 'endind')


function [P, n, s]=resectioning(x, X, N)
% x are the measured points and X are the 3D corresponding points, in homogenous coordinates. 
% N is the normalization matrix
% return the camera matrix P, minimal norm n and minimal singular value s
    n = length(X);
    x_norm = N * x;
    M = zeros(3 * n, 12 + n);
    for i=1:n
        row = 3 * i - 2;
        M(row, 1:4) = X(:, i);
        M(row+1, 5:8) = X(:, i);
        M(row+2, 9:12) = X(:, i);
        col = 12 + i;
        M(row:row+2, col) = -x_norm(:, i);
    end

    % solve using SVD & set up the camera matrix
    [~,S,V] = svd(M);
    v = V(:, end);
    P = reshape(v(1:12) ,[4 3])';
    
    % re-normalize
    P = inv(N) * P;
    
    % find the correct solution s.t. the points are in front of the camera
    proj = P*X;
    P = P * sign(proj(3,1));

    % calculate the norm and minimal singular value
    n = norm(M * v);
    s = min(S(S > 0));
end