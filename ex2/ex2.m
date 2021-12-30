% Dana Joffe 312129240
load compEx1data.mat

% find K, R for camera P1
[K,R] = rq(P{1}(1:3,1:3));
K = K./K(3,3);
t = inv(K)*P{1}(:,4);

% find K, R for projected camera P1 * T1^-1
P1_T1 = P{1} * inv(T1);
[K1,R1] = rq(P1_T1(1:3,1:3));
K1 = K1./K1(3,3);
t1 = inv(K1)*P1_T1(:,4);

% find K, R for projected camera P1 * T2^-1
P1_T2 = P{1} * inv(T2);
[K2,R2] = rq(P1_T2(1:3,1:3));
K2 = K2./K2(3,3);
t2 = inv(K2)*P1_T2(:,4);
