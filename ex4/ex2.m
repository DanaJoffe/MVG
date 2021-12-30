% Dana Joffe 312129240

load compEx2data.mat

%% Find the homography from a to b: b ~ Ha

x1 = [xA; ones(1, length(xA))];
x2 = [xB; ones(1, length(xB))];

% use RANSAC to find inliers
k=4;
N=40;
dist_th = 5;
con_set = [];
ransac_h = [];

for i=1:N
   % form an homography using a random minimal set of points
   randind = randi([1 length(x1)],1,k);
   H = DLT_homography(x1(:, randind), x2(:, randind));
   
   % evaluate the homography
   x2_est = pflat(H*x1);
   distances = sqrt(sum((x2 - x2_est).^2));
   inliers = distances <= dist_th;
   
   if sum(inliers) > sum(con_set)
       con_set = inliers;
       ransac_h = H;
   end
end

inliers_amount = sum(con_set)
x1_in = x1(:,con_set);
x2_in = x2(:,con_set);
H = DLT_homography(x1_in, x2_in);

%% Transform the images to a common coordinate system

A = imread('a.jpg');
B = imread('b.jpg');
bestH = H';

tform = maketform('projective',bestH); % Creates a transfomation that matlab can use for images. Notice that H is transposed.
transfbounds = findbounds(tform ,[1, 1; size(A,2), size(A,1)]); % Finds the bounds of the transformed image

% Computes bounds of a new image such that both the old ones will fit.
xdata = [min([transfbounds(: ,1); 1]) max([transfbounds(: ,1); size(B ,2)])];
ydata = [min([transfbounds(: ,2); 1]) max([transfbounds(: ,2); size(B ,1)])];

% Put B inside bigger boundaries for both images
[ newA ] = imtransform(A,tform ,'xdata',xdata ,'ydata',ydata ); % Transform the image using bestH
tform2 = maketform ('projective',eye(3));
[ newB ] = imtransform(B,tform2 ,'xdata',xdata ,'ydata',ydata ,'size',size(newA));

% Writes both images in the new image.
newAB = newB;
newAB( newB < newA ) = newA( newB < newA );

figure
imagesc(newAB)
title("2. Stitched images")


function N=norm_mat(x)
% return: normalization matrix
    m = mean(x(1:2 ,:) ,2);
    s = std(x(1:2 ,:) ,0 ,2);
    N = [1/s(1), 0, -m(1)/s(1); 0, 1/s(2), -m(2)/s(2); 0 0 1];
end

function H=DLT_homography(x1, x2)
% param x1: matrix, columns contain homogenous points in 2D.
% param x2: matrix, columns contain the corresponding points to x1.
% return: homography from x1 to x2: x2 ~ H * x1
    N1 = norm_mat(x1);
    N2 = norm_mat(x2);
    x1n = N1 * x1;
    x2n = N2 * x2;
    
    n = length(x1);
    M = zeros(3 * n, 9 + n);
    for i=1:n
        row = 3 * i - 2;
        M(row, 1:3) = x1n(:, i)';
        M(row+1, 4:6) = x1n(:, i)';
        M(row+2, 7:9) = x1n(:, i)';
        col = 9 + i;
        M(row:row+2, col) = -x2n(:, i);
    end

    % solve using SVD & set up the homography
    [~,S,V] = svd(M);
    v = V(:, end);
    Hn = reshape(v(1:9) ,[3 3])';
    
    % re-normalize
    H = inv(N2) * Hn * N1;
    
%     % calculate the norm and minimal singular value
%     n = norm(M * v);
%     s = min(S(S > 0));
end