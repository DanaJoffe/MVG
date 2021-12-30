% Dana Joffe 312129240

%% Computer Exercise 3
load compEx1data.mat
load Ex1results.mat
load compEx3data.mat

% calc essential matrix E
Kin = inv(K);
y = {Kin * x{1}, Kin * x{2}};

[E, sv, n] = get_essential_matrix(y);
estimated_E = E / E(3,3);

% check that the epipolar constraints are roughly fulfilled
figure
plot(diag(y{2}' * E * y{1}));
title("3.1 Epipolar constraints for E")
xlabel("Points index")
ylabel("Epipolar constraint value")

% compute the fundamental matrix from the essential matrix
F = Kin' * E * Kin;
F = F / F(3,3);

% draw epipolar lines on image 2
l = F * x{1};
ind = randi([1 length(x{1})],1,20);
x2_samples = x{2}(:, ind);
epip_lines = l(:, ind);

im = imread('kronan2.JPG');
figure
imagesc(im)
hold on
plot(x2_samples(1 ,:) ,x2_samples(2 ,:) , 'y.', 'MarkerSize', 12);
hold on
rital(epip_lines, 'g')
colormap gray
axis equal
title("3.2 Epipolar lines on image 2")

% distance between all the points and their corresponding epipolar lines
l = l./ sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3 1]));
distances = abs(sum(l.*x{2}));
mean_dist = mean(distances);

figure
hist(distances ,100);
title("3.3 Distance between points and their corresponding epipolar lines")
xlabel("Distance from epipolar line")
ylabel("# of points")


%% Computer Exercise 4

% create a valid essential matrix from an approximate solution
[U,S,V] = svd(E);
if det(U*V') >0
    E = U*diag ([1 1 0])*V';
else
    V = -V;
    E = U* diag([1 1 0])*V';
end

% compute four camera solutions
u3 = U(:, 3);
W = [0, -1, 0; 1, 0, 0; 0, 0, 1];

P2s = {[U*W*V' u3], [U*W'*V' u3], [U*W*V', -u3], [U*W'*V', -u3]};
P1 = [eye(3) [0 0 0]'];

% triangulate for each of the four camera solutions
X = cell(1, 4);
best = 0;
max_points_infront = -1;
for i=1:4
   X{i} = triangulate(y, {P1, P2s{i}});
   
   % pick the best solution
   proj = [P1 * X{i} P2s{i} * X{i}];
   points_in_front = length(find(proj(3, :)>0));
   if points_in_front > max_points_infront
       max_points_infront = points_in_front;
       best = i;
   end
end

P = {K * P1, K * P2s{best}};
X = X{best};

% plot image the points and the projected 3D-points i
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
title("4.1 Triangulation using E: project on camera 1")

figure
imagesc(im{2})
colormap gray
hold on
axis equal
plot(proj{2}(1, :), proj{2}(2, :),'yo');
hold on
plot(x{2}(1, :), x{2}(2 , :),'c.');
legend('projected points', 'measured points')
title("4.2 Triangulation using E: project on camera 2")

figure
plot3(X(1 ,:) ,X(2 ,:) ,X(3 ,:) , '.')
hold on
plotcams(P)
grid on
title("4.3 Triangulation using E: 3D points")


function X=triangulate(x, P)
% ========> EX4 <========
% Param: x - 2 cells of image points correspondences.
%        P - 2 cells of cameras.
% Return all 3D points, corresponding to the image points in x.

    % normalize
    N = {get_normalization_matrix(x{1}), get_normalization_matrix(x{2})};
    x = {N{1} * x{1}, N{2} * x{2}};
    P = {N{1} * P{1}, N{2} * P{2}};

    n = length(x{1});
    X = zeros(4,n);
    for i=1:n
        X(:,i) = triangulation(P, [x{1}(:,i) x{2}(:,i)]);
    end
end

function X=triangulation(P, x)
% ========> EX4 <========
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
function [E, sv, nor]=get_essential_matrix(x)
% ========> EX3 <========
% Param: x - 2 cells of corresponding image points
% Return: E - essential matrix, sv - singular value of the SVD solution, 
%         nor - norm of Mv.
    len = length(x{1});
    M = zeros(len, 9);
    for i=1:len
        xx = x{2}(:,i) * x{1}(:,i)';
        M(i ,:) = xx(:)';
    end

    % solve using SVD & set up the normalized fundamental matrix
    [~,S,V] = svd(M);
    v = V(:, end);
    E = reshape(v ,[3 3]);

    % calculate the norm and minimal singular value
    nor = norm(M * v);
    sv = min(S(S > 0));

    % ensure that det(E)=0 and that E has two nonzero equal singular values
    [U,S,V] = svd(E);
%     s = (S(1,1) + S(2,2)) / 2;
    s = 1;
    S = diag([s, s, 0]);
    E = U * S * V';
end

function N=get_normalization_matrix(x)
% ========> EX4 <========
% Param: x is a matrix with columns being homogenous coordinates of 2D points.
% Return: a normalization matrix that subtract the mean and re-scale using 
% the standard deviation.
    m = mean(x(1:2 ,:) ,2);
    s = std(x(1:2 ,:) ,0 ,2);
    N = [1/s(1), 0, -m(1)/s(1); 0, 1/s(2), -m(2)/s(2); 0 0 1];
end