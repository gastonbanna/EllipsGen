%% main.m - Ellipsoid Pair Data Generation and Export
% 
% This script serves as the main code to run in order to generate and store
% geometric and transformation data for pairs of ellipsoidal particles. 
% The generated dataset is saved to a CSV file for further analysis. 
% The script can also visualizes a selected pair of ellipses.
%
% Author: Gaston Banna
% Date: 28 Septembre 2023
%
% -------------------------------------------------------------------------
% Functionality:
% - Generates N random pairs of ellipsoidal particles.
% - Computes their geometric parameters and transformations.
% - Stores the results in a structured data matrix.
% - Saves the data as 'data.csv'.
% - Visualizes a randomly selected pair of ellipses.
%
% Dependencies:
% - Configurations.m: Computes ellipsoid configurations.
% - Matrices.m: Computes transformation matrices.
% - ellipsegraph.m: Plots ellipses.
% - test_appartenance.m: Tests point membership in ellipses.
%
% Output:
% - 'data.csv' : Table of generated ellipse parameters.
%
% -------------------------------------------------------------------------

%% Initialization
N = 10000; % Number of ellipse pairs

% Preallocate arrays for efficiency
A = zeros(2, 2);
B = zeros(1, 10);
x = zeros(N, 2);
y = zeros(N, 2);
param = zeros(N, 10);
data = zeros(N, 38);

%% Generate Random Ellipse Pairs
for i = 1:N
    % Randomly generate ellipse parameters
    gamma_i = 1 + 19*rand();
    omega_i = 1 + 19*rand();
    theta_i = pi*rand();
    gamma_j = 1 + 19*rand();
    theta_j = pi*rand();
    c_j = rand(2, 1); % Center of the second ellipse
    phi = -pi + 2*pi*rand(); % Random rotation angle
    epsilon_bar = sign(-1 + 2*rand())*10^(-(1 + 1*rand())); % "Distance"
    
    % Compute ellipse configurations
    [E_i, E_j, x_i, x_j, epsilon] = Configurations(gamma_i, omega_i, ...
        theta_i, gamma_j, theta_j, c_j, phi, epsilon_bar);
    
    % Store position and parameters
    A = [x_i, x_j];
    B = [E_i, E_j];
    x(i, :) = A(1, :);
    y(i, :) = A(2, :);
    param(i, :) = B;
    
    % Store ellipse parameters and computed properties in the dataset
    data(i, 1:5) = E_i;
    data(i, 6:7) = [gamma_i, omega_i];
    data(i, 8:9) = [cos(theta_i), sin(theta_i)];
    [Q_i, o_i] = Matrices(E_i); % Compute transformation matrix
    data(i, 10:12) = [Q_i(1,1), Q_i(1,2), Q_i(2,2)];
    data(i, 13:16) = [norm(o_i), atan2(o_i(2), o_i(1)), ...
        cos(atan2(o_i(2), o_i(1))), sin(atan2(o_i(2), o_i(1)))];
    
    data(i, 17:21) = E_j;
    data(i, 22) = gamma_j;
    data(i, 23:25) = [1, cos(theta_j), sin(theta_j)];
    [Q_j, o_j] = Matrices(E_j);
    data(i, 26:28) = [Q_j(1,1), Q_j(1,2), Q_j(2,2)];
    data(i, 29:32) = [norm(o_j), atan2(o_j(2), o_j(1)), ...
        cos(atan2(o_j(2), o_j(1))), cos(atan2(o_j(2), o_j(1)))];
    
    % Store contact point and perturbation information
    data(i, 33:38) = [x_i(1), x_i(2), x_j(1), x_j(2), ...
        sign(epsilon_bar)*norm(x_j - x_i), epsilon_bar];
end

%% Export Data to CSV File
parameters = {'a_i', 'b_i', 'theta_i', 'o_i_x', 'o_i_y', 'gamma_i', ...
    'omega_i', 'cos_theta_i', 'sin_theta_i', 'Q_i_11', 'Q_i_12', ...
    'Q_i_22', 'r_i', 'phi_i', 'cos_phi_i', 'sin_phi_i', 'a_j', 'b_j', ...
    'theta_j', 'o_j_x', 'o_j_y', 'gamma_j', 'omega_j', 'cos_theta_j', ...
    'sin_theta_j', 'Q_j_11', 'Q_j_12', 'Q_j_22', 'r_j', 'phi_j', ...
    'cos_phi_j', 'sin_phi_j', 'x_i_x', 'x_i_y', 'x_j_x', 'x_j_y', ...
    'epsilon','epsilon_bar'};
T = array2table(data);
T.Properties.VariableNames = parameters;
writetable(T,'data.csv');

%% Visualization: Display a Random Pair of Ellipses
paire = randi(N);
E_i = param(paire, 1:5);
E_j = param(paire, 6:10);
x_i = [x(paire, 1); y(paire, 1)];
x_j = [x(paire, 2); y(paire, 2)];
f = ellipsegraph(1, E_i, E_j, 3, x_i, x_j);

% Test point membership for the displayed pair
[p1inE, p1inH, p1ok] = membership_test(x_i, E_i, E_j);
[p2inE, p2inH, p2ok] = membership_test(x_j, E_j, E_i);


%3D2
