% % Author: Suryansh Kumar
% % ETH Zurich


function [R, S] = giveme_motion_and_shape(W)

% Step 1: substract the translation component.
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

% Note that one can use other optimization methods to solve
% for Q and may end up getting a better solution.

end



