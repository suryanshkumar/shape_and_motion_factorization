% Author: Suryansh Kumar
% ETH Zurich

% important path for the code.
dataset_folder = './dataset/';
utility_folder = './util/';
src_folder = './src/';
dataset_file = 'point_trajectory.mat';

% add the path
addpath(dataset_folder);
addpath(utility_folder);
addpath(src_folder);

% get the measurment matrix
W = giveme_measurement_matrix(dataset_folder, dataset_file);

% solve for the shape using tomasi and kannade 1992 algorithm.
[R, S] = giveme_motion_and_shape(W);

% plot the reconstructed 3D points.
figure,
plot3(S(1, :), S(2, :), S(3, :), 'r.');