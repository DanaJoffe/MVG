% Dana Joffe 312129240

%% a - Total Least Squares

load compEx1data.mat

X_ = pflat(X);
plane = total_least_squares(X_);
rms_a = RMS(plane,X_)

%% b - RANSAC for plane fitting
k=3;
N=10;
dist_th = 0.1;
con_set = [];
ransac_plane = [];
for i=1:N
   % form a plane using a random minimal set of points
   randind = randi([1 length(X_)],1,k);
   plane = null(X_(:, randind)'); % Computes a plane from a sample set
   plane = plane ./ norm(plane(1:3)); % normalize the plane to have length 1
   
   % evaluate plane
   inliers = abs(plane'*X_) <= dist_th;
   
   if sum(inliers) > length(con_set)
       con_set = X_(:,inliers);
       ransac_plane = plane;
   end
end

% find a least square solution on the consensus set
rms_b = RMS(ransac_plane,con_set)
inliers_amount = length(con_set)

% plot the absolute distances between the plane and the points
figure
hist(abs(ransac_plane' * con_set) ,100);
title("1.b Distance between RANSAC plane and inliers")
xlabel("Distance from the plane")
ylabel("# of points")

%% c

% find a least square solution on the consensus set
plane = total_least_squares(con_set);
rms_c = RMS(plane,con_set)

% plot the absolute distances between the plane and the points
figure
hist(abs(plane' * con_set) ,100);
title("1.c Distance between TLS plane and inliers")
xlabel("Distance from the plane")
ylabel("# of points")

% plot the projection of the inliers into the images
x_proj = pflat(P{1} * con_set);
im1 = imread('house1.jpg');
figure
imagesc(im1)
hold on
plot(x_proj(1, :), x_proj(2, :),'y.');
colormap gray
axis equal
title("1.c Inliers projected into image 1")

x_proj = pflat(P{2} * con_set);
im2 = imread('house2.jpg');
figure
imagesc(im2)
hold on
plot(x_proj(1, :), x_proj(2, :),'y.');
colormap gray
axis equal
title("1.c Inliers projected into image 2")

%% d

% compute a homography from camera 1 to camera 2
P_ = {inv(K) * P{1}, inv(K) * P{2}};

plane_pi = pflat(plane);
pi_ = plane_pi(1:3);
H = P_{2}(1:3,1:3) - P_{2}(:,4) * pi_';
x_trans = K * pflat(H * pflat(inv(K) * x));

figure
imagesc(im1)
hold on
plot(x(1, :), x(2, :),'y.', 'MarkerSize', 12);
colormap gray
axis equal
title("1.d  Points x in image 1")

figure
imagesc(im2)
hold on
plot(x_trans(1, :), x_trans(2, :),'y.', 'MarkerSize', 12);
colormap gray
axis equal
title("1.d  Transformed points in image 2")


function dist=RMS(plane,X)
% Compute the RMS distance between 3D-points and a plane.
% param X: 3D points in homogeneuse coordinates.
% param plane: normalized plane.
%     X = pflat(X);
%     plane = plane ./ norm(plane(1:3));
    dist = sqrt(sum((plane'*X).^2)/ size(X ,2));
end
function plane=total_least_squares(X)
% Solve the total least squares problem.
% param X: 3D points in homogeneuse coordinates
% return: the plane's coefficients, after normalization.
    mean_X = mean(X,2);
    Xtilde = X - mean_X;
    M = Xtilde(1:3,:) * Xtilde(1:3 ,:)';
    [V,D] = eig(M);
    [~,min_ind] = min(diag(D));
    t = V(:,min_ind);
    d = - t' * mean_X(1:3);
    plane = [t; d];
    plane = plane ./ norm(plane(1:3)); % normalize the plane to have length 1
end