%% Dana Joffe 312129240
%% Computer Exercise 3
load compEx3.mat

figure
plot([startpoints(1 ,:); endpoints(1 ,:)], [startpoints(2 ,:); endpoints(2 ,:)], 'b-');
title("3.1 original grid")
axis equal
 
H1 = [sqrt(3) -1 1 ; 1 sqrt(3) 1 ; 0 0 2];
H2 = [1 -1 1 ; 1 1 0 ; 0 0 1];
H3 = [1 1 0 ; 0 2 0 ; 0 0 1];
H4 = [sqrt(3) -1 1 ; 1 sqrt(3) 1 ; 1/4 1/2 2];


homo_endpoints = [endpoints; ones(1, length(endpoints))];
homo_startpoints = [startpoints; ones(1, length(startpoints))];

% H1
trans_end = pflat(H1*homo_endpoints);
trans_start = pflat(H1*homo_startpoints);

figure
subplot(2, 2, 1)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H1 - Euclidean")
axis equal

% H2
trans_end = pflat(H2*homo_endpoints);
trans_start = pflat(H2*homo_startpoints);

% figure
subplot(2, 2, 2)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H2 - Similarity")
axis equal

% H3
trans_end = pflat(H3*homo_endpoints);
trans_start = pflat(H3*homo_startpoints);

% figure
subplot(2, 2, 3)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H3 - Affine")
axis equal

% H4
trans_end = pflat(H4*homo_endpoints);
trans_start = pflat(H4*homo_startpoints);

% figure
subplot(2, 2, 4)
plot([trans_start(1 ,:); trans_end(1 ,:)], [trans_start(2 ,:); trans_end(2 ,:)], 'b-');
title("H4 - Projective")
axis equal
suptitle('3.2 Transformations')
