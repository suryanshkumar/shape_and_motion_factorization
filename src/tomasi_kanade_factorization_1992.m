% Author: Suryansh Kumar
% ETH Zurich

% add the dataset and important function file path.
dataset_folder = '../dataset/';
utility_folder = '../util/';
dataset_file = 'point_trajectory.mat';

addpath(dataset_folder);
addpath(utility_folder);

% Step 1: get the measurment matrix
W = giveme_measurement_matrix(dataset_folder, dataset_file);

% Step 2: substract the translation component.
Wm = mean(W, 2);
W_tilde = W - Wm*ones(1, size(W, 2));

% Step 3: perform the rank 3 factorization of the measurement matrix.
[O1, S, O2] = svd(W_tilde);

    % First three columns of O1. (%This part is intensionally indented)
    O1_prime = O1(:, 1:3);
    % Submatrix containing the top three singular values of S.
    S1_prime = S(1:3, 1:3);
    % First three rows of O2' matrix.
    O2 = O2';
    O2_prime = O2(1:3, :);
    
% Step 4: create the R_hat and S_hat matrix.
R_hat = O1_prime*sqrtm(S1_prime);
S_hat = sqrtm(S1_prime)*O2_prime;

% Here, comes the complex part of the algorithm.
% The above decomposition is not unique, but the 
% point to note is that, its a linear transformed 
% bases of the same space. Hence, we must solve for such a transformation.

% Step 5: metric constraint. (Compute Q)

% Note that Q*Q^T is symmetric. Thus we have to solve for only six 
% variables under the orthonormal constraint.

% solve for each frame.
L = [];
Ro = [];

for i = 1:size(R_hat, 1)/2
    G_k1 = R_hat(2*i-1, :);
    a = G_k1(1); b = G_k1(2); c = G_k1(3);
   
    G_k2 = R_hat(2*i, :);
    d = G_k2(1); e = G_k2(2); f = G_k2(3);
   
    L = [L;
         a*a a*b+b*a a*c+c*a b*b c*b+b*c c*c;
         d*d d*e+e*d d*f+f*d e*e f*e+e*f f*f;
         a*d a*e+b*d a*f+c*d b*e b*f+e*c c*f];
       
    Ro = [Ro; 1; 1; 0];
end

QQT = L \ Ro;

LT = [QQT(1), QQT(2), QQT(3);
      QQT(2), QQT(4), QQT(5);
      QQT(3), QQT(5), QQT(6)];

% metric correction matrix.
Q = chol(LT);

% Step 6: recover rotation and structure
R = R_hat*Q;
S = Q \ S_hat;


% Misc: plot the 3D points.
figure,
plot3(S(1, :), S(2, :), S(3, :), 'r.');



