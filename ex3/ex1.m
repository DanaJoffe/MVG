% Dana Joffe 312129240

load compEx1data.mat

% calc normalization matrices
N = {get_normalization_matrix(x{1}), get_normalization_matrix(x{2})};
x_norm = {N{1} * x{1}, N{2} * x{2}};

% calc the fundamental matrix
[F, Fn, s, n] = get_fundamental_matrix(x, N);
F

figure
plot(diag(x_norm{2}' * Fn * x_norm{1}));
title("1.1 Epipolar constraints for normalized F")
xlabel("Points index")
ylabel("Epipolar constraint value")

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
title("1.2 Epipolar lines on image 2")

% distance between all the points and their corresponding epipolar lines
l = l./ sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3 1]));
distances = abs(sum(l.*x{2}));
mean_dist = mean(distances);

figure
hist(distances ,100);
title("1.3 Distance between points and their corresponding epipolar lines")
xlabel("Distance from epipolar line")
ylabel("# of points")

% Repeat the question without normalization
[F_, ~, ~, ~] = get_fundamental_matrix(x, {eye(3), eye(3)});
l = F_ * x{1};
l = l./ sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3 1]));
distances = abs(sum(l.*x{2}));
mean_dist2 = mean(distances);

%% save data for Computer Exercise 2
% save('Ex1results.mat','F', 'N')

%% functions

function [F, Fn, s, n]=get_fundamental_matrix(x, N)
% Param: x - 2 cells of corresponding image points, N - 2 cells of
% normalization matrices.
% Return: F - fundamental matrix, Fn - normalized fundamental matrix, s -
% singular value of the SVD solution, n - norm of Mv
    x_norm = {N{1} * x{1}, N{2} * x{2}};
    len = length(x{1});
    M = zeros(len, 9);
    for i=1:len
        xx = x_norm{2}(:,i) * x_norm{1}(:,i)';
        M(i ,:) = xx(:)';
    end

    % solve using SVD & set up the normalized fundamental matrix
    [~,S,V] = svd(M);
    v = V(:, end);
    Fn = reshape(v ,[3 3]);

    % calculate the norm and minimal singular value
    n = norm(M * v);
    s = min(S(S > 0));

    % enforce det(Fn)=0 
    [U,S,V]=svd(Fn);
    S(3,3)=0;
    Fn = U*S*V';

    % compute the un-normalized fundamental matrix F
    F = N{2}' * Fn * N{1};
    F = F / F(3,3);
end
function N=get_normalization_matrix(x)
% Param: x is a matrix with columns being homogenous coordinates of 2D points.
% Return: a normalization matrix that subtract the mean and re-scale using 
% the standard deviation.
    m = mean(x(1:2 ,:) ,2);
    s = std(x(1:2 ,:) ,0 ,2);
    N = [1/s(1), 0, -m(1)/s(1); 0, 1/s(2), -m(2)/s(2); 0 0 1];
end