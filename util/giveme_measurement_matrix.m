% Author: Suryansh Kumar
% ETH Zurich

function [W] = giveme_measurement_matrix(folder_path, file_name)

file_path = fullfile(folder_path, file_name);

if exist(file_path, 'file')
    point_trajectory = load(file_path);
else
    sprintf('file %s does not exits', file_path);
end

% read the x and y corrdinates of the points in all the frames.
U = point_trajectory.Xs;
V = point_trajectory.Ys;

% Create a measurement matrix.

% This step is slightly different from the paper.
% In the paper, it mentions the augumention of U and V to
% create 2F X P matrix. Nonetheless, it does not changes the overall idea
% of the algorithm. Our choice is inspired by many extention of this 
% factorization mathod that uses alternate u, v rows rather than augmenting
% U and V.
[F, P] = size(U);
W = zeros(2*F, P);

W(1:2:end, :) = U;
W(2:2:end, :) = V;

end