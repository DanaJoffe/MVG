% Dana Joffe 312129240
load compEx1data.mat

%% a
figure
plot3(X(1 ,:) ,X(2 ,:) ,X(3 ,:) , '.','Markersize' ,2);
hold on
plotcams(P)
title("1.a Structure & motion uncalibrated reconstruction")
axis equal

%% b
i = 1;
P1 = P{i};
proj_pt = pflat(P1 * X);
im_pt1 = x{i};
im1 = imread(imfiles{i});
visible = isfinite(x{i}(1 ,:)); % Determines which of the points are visible in image i

figure
imagesc(im1)
title("1.b camera 1")
colormap gray
hold on

plot(proj_pt(1, visible ), proj_pt (2, visible ),'ro'); % Plots a red ’o’ at each visible point in xproj
plot(im_pt1(1, visible ), im_pt1(2 , visible ),'*'); % Plots a ’*’ at each point coordinate

axis equal
legend('projected points', 'image points')

%% c

% T1 transformation
T1_inv = inv(T1);
X_new1 = pflat(T1 * X);
P_new1 = cell(size(P));
for i=1:length(P)
   P_new1{i} = P{i} * T1_inv;
end

figure
plot3(X_new1(1 ,:) ,X_new1(2 ,:) ,X_new1(3 ,:) , '.','Markersize' ,2);
hold on
plotcams(P_new1)
title("1.c Structure & motion after T1 transformation")
axis equal

% T2 transformation
T2_inv = inv(T2);
X_new2 = pflat(T2 * X);
P_new2 = cell(size(P));
for i=1:length(P)
   P_new2{i} = P{i} * T2_inv;
end

figure
plot3(X_new2(1 ,:) ,X_new2(2 ,:) ,X_new2(3 ,:) , '.','Markersize' ,2);
hold on
plotcams(P_new2)
title("1.c Structure & motion after T2 transformation")
axis equal

%% d

proj_pt_T1 = pflat(P1 * X_new1);

figure
imagesc(im1)
title("1.d Project T1*X into P1")
colormap gray
hold on

plot(proj_pt_T1(1, visible ), proj_pt_T1 (2, visible ),'ro'); % Plots a red ’o’ at each visible point in xproj
plot(im_pt1(1, visible ), im_pt1(2 , visible ),'*'); % Plots a ’*’ at each point coordinate

axis equal
legend('projected T1 points', 'image points')

% project into transformed camera 1
proj_pt_T1 = pflat(P_new1{1} * X_new1);
figure
imagesc(im1)
title("1.d Project T1*X into P1*inv(T1)")
colormap gray
hold on

plot(proj_pt_T1(1, visible ), proj_pt_T1 (2, visible ),'ro'); % Plots a red ’o’ at each visible point in xproj
plot(im_pt1(1, visible ), im_pt1(2 , visible ),'*'); % Plots a ’*’ at each point coordinate

axis equal
legend('projected T1 points', 'image points')
