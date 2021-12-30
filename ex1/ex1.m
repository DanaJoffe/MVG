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
