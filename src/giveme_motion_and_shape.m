% % Author: Suryansh Kumar
% % ETH Zurich


function [R, S] = giveme_motion_and_shape(W)

% Step 1: mean centralize the meaurement matrix.
Wm = mean(W, 2);
W_tilde = W - Wm*ones(1, size(W, 2));

% Step 2: perform the rank 3 factorization of the measurement matrix.
[O1, S, O2] = svd(W_tilde);

    % First three columns of O1. (%This part is intensionally indented)
    O1_prime = O1(:, 1:3);
    % Submatrix containing the top three singular values of S.
    S1_prime = S(1:3, 1:3);
    % First three rows of O2' matrix.
    O2 = O2';
    O2_prime = O2(1:3, :);
    
% Step 3: create the R_hat and S_hat matrix.
R_hat = O1_prime*sqrtm(S1_prime);
S_hat = sqrtm(S1_prime)*O2_prime;

% Notes: Here, comes the complex part of the algorithm.
% The above decomposition is not unique, but the 
% point to note is that, its a linear transformed 
% bases of the same space. Hence, we must solve 
% for such a transformation.

% Step 4: metric constraint. (Compute Q)
% Note that Q*Q^T is symmetric. Accordingly,
% we first solve for the six variables under the orthonormality 
% constraint and later use cholesky factorization to recover Q.

% solve for each frame.
L = [];
Ro = [];

for i = 1:size(R_hat, 1)/2
    G_k1 = R_hat(2*i-1, :);
    r11 = G_k1(1); r12 = G_k1(2); r13 = G_k1(3);
   
    G_k2 = R_hat(2*i, :);
    r21 = G_k2(1); r22 = G_k2(2); r23 = G_k2(3);
   
    L = [L;
         r11*r11 r11*r12+r12*r11 r11*r13+r13*r11 r12*r12 r13*r12+r12*r13 r13*r13;
         r21*r21 r21*r22+r22*r21 r21*r23+r23*r21 r22*r22 r23*r22+r22*r23 r23*r23;
         r11*r21 r11*r22+r12*r21 r11*r23+r13*r21 r12*r22 r12*r23+r22*r13 r13*r23];
       
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

% Note that one can use other optimization methods to solve
% for Q and may end up getting a better solution.

end



