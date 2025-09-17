%% Configurations.m is a function that generates a pair of ellipses and 
%% computes their contact points
%
% -------------------------------------------------------------------------
% Inputs:
%   gamma_i, omega_i - Shape parameters of the first ellipse
%   theta_i - Orientation of the first ellipse
%   gamma_j - Shape parameter of the second ellipse
%   theta_j - Orientation of the second ellipse
%   o_j - Center of the second ellipse
%   phi - Contact angle parameter
%   epsilon_bar - Distance parameter, using the contact point on (i)
%
% Outputs:
%   E_i - Parameters of the first ellipse [a, b, theta, o_x, o_y]
%   E_j - Parameters of the second ellipse [a, b, theta, o_x, o_y]
%   x_i - Contact point on the first ellipse
%   x_j - Contact point on the second ellipse
%
% -------------------------------------------------------------------------

function [E_i, E_j, x_i, x_j, epsilon] = Configurations(gamma_i, ...
    omega_i, theta_i, gamma_j, theta_j, o_j, phi, epsilon_bar)

    % Rotation matrix function
    R = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
    % Diagonal matrix function for ellipse shape
    D = @(a, b) [1 / a^2, 0; 0, 1 / b^2];

    % Compute semi-axes of the ellipses
    a_i = sqrt(omega_i * gamma_i);
    b_i = sqrt(omega_i / gamma_i);
    a_j = sqrt(gamma_j);
    b_j = 1 / a_j;

    % Compute the transformed shape matrix of the first ellipse (i) in the 
    % second ellipse (j)'s frame
    Q_i_hat = D(1/sqrt(a_j), 1/sqrt(b_j)) * R(-theta_j) * R(theta_i) * ...
        D(a_i, b_i) * R(-theta_i) * R(theta_j) * D(1/sqrt(a_j), ...
        1/sqrt(b_j));
    [a_i_hat, b_i_hat, theta_i_hat] = parametrize(Q_i_hat);

    % Compute the initial contact point in the local frame of (i)
    x_i_bar = [a_i_hat * cos(phi); b_i_hat * sin(phi)];
    g_x_i_b = [b_i_hat * cos(phi); a_i_hat * sin(phi)];
    n_i_bar = g_x_i_b / norm(g_x_i_b); % Normal vector at the contact point
    o_j_bar = x_i_bar + (1 + epsilon_bar) * n_i_bar; % Second ellipse pos

    % Transformations to global coordinates
    T1 = @(x_bar) R(theta_i_hat) * (x_bar - o_j_bar);
    T2 = @(x_hat) R(theta_j) * D(1/sqrt(a_j), 1/sqrt(b_j)) * x_hat + o_j;

    % Compute center of the first ellipse in the transformed frame
    o_i_bar = [0; 0];
    o_i_hat = T1(o_i_bar);
    o_i = T2(o_i_hat);

    % Define ellipse parameters in global coordinates
    E_i = [a_i, b_i, theta_i, o_i(1), o_i(2)];
    E_j = [a_j, b_j, theta_j, o_j(1), o_j(2)];

    % Solve for contact points on the ellipses
    [x_i, x_j] = Contact_Solve(E_i, E_j);

    % Compute the separation distance
    epsilon = norm(x_j - x_i);
end


function [E_i, E_j, varargout] = Configurations_univ(dimension, ratios_i, angles_i, positions_i, ratios_j, angles_j, positions_j)

% Warnings
if size(ratios_i, 2) ~= size(angles_i, 2) || size(ratios_i, 2) ~= size(positions_i, 2) || size(ratios_i, 2) ~= size(ratios_j, 2) || size(ratios_i, 2) ~= size(angles_j, 2) || size(ratios_i, 2) ~= size(positions_j, 2)
    error('Not enough input arguments.')
else
    N = size(ratios_i, 2);
end

switch dimension
    case '2D'
        % Rotation matrix function
        R = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
        % Diagonal matrix function for ellipse shape
        D = @(a, b) [1 / a^2, 0; 0, 1 / b^2];

        % Read the input parameters
        gamma_i = ratios_i(2, :);
        omega_i = ratios_i(1, :);
        theta_i = angles_i(1, :);
        gamma_j = ratios_j(2, :);
        omega_j = ratios_j(1, :);
        theta_j = angles_j(1, :);
        o_j = positions_j(1:2, :);
        phi = positions_i(1, :);
        epsilon_bar = positions_i(2, :);

        % Compute semi-axes of the ellipses
        a_i = sqrt(omega_i .* gamma_i);
        b_i = sqrt(omega_i ./ gamma_i);
        a_j = sqrt(omega_j .* gamma_j);
        b_j = 1 / a_j;

        % Compute the transformed shape matrix of the first ellipse (i) in the second ellipse (j)'s frame
        Q_i_hat = cell(1, N);
        a_i_hat = zeros(1, N);
        b_i_hat = zeros(1, N);
        theta_i_hat = zeros(1, N);
        for k = 1:N
            Q_i_hat(1, k) = D(1/sqrt(a_j(k)), 1/sqrt(b_j(k))) * R(-theta_j(k)) * R(theta_i(k)) * D(a_i(k), b_i(k)) * R(-theta_i(k)) * R(theta_j(k)) * D(1/sqrt(a_j(k)), 1/sqrt(b_j(k)));
            [a_i_hat(k), b_i_hat(k), theta_i_hat(k)] = parametrize(Q_i_hat(1, k));
        end

        % Compute the initial contact point in the local frame of (i)
        x_i_bar = [a_i_hat .* cos(phi); b_i_hat .* sin(phi)];
        g_x_i_b = [b_i_hat .* cos(phi); a_i_hat .* sin(phi)];
        norms = sqrt(sum(g_x_i_b.^2, 1));
        n_i_bar = g_x_i_b ./ norms; % Normal vector at the contact point
        o_j_bar = x_i_bar + (1 + epsilon_bar) * n_i_bar; % Second ellipse pos

        % Transformations to global coordinates
    T1 = @(x_bar) R(theta_i_hat) * (x_bar - o_j_bar);
    T2 = @(x_hat) R(theta_j) * D(1/sqrt(a_j), 1/sqrt(b_j)) * x_hat + o_j;

    % Compute center of the first ellipse in the transformed frame
    o_i_bar = [0; 0];
    o_i_hat = T1(o_i_bar);
    o_i = T2(o_i_hat);

    % Define ellipse parameters in global coordinates
    E_i = [a_i, b_i, theta_i, o_i(1), o_i(2)];
    E_j = [a_j, b_j, theta_j, o_j(1), o_j(2)];

    % Solve for contact points on the ellipses
    [x_i, x_j] = Contact_Solve(E_i, E_j);

    % Compute the separation distance
    epsilon = norm(x_j - x_i);



    case '3D'
        % Rotation matrix function
        R = @(theta_x, theta_y, theta_z) eul2rotm([theta_z; theta_y; theta_x]);
        % Diagonal matrix function for ellipse shape
        D = @(a, b, c) [1 / a^2, 0, 0; 0, 1 / b^2, 0; 0, 0, 1 / c^2];
    otherwise
end





%% Algorithme de generation d'une paire d'ellipse et du point de contact construit

    R = @(theta_x, theta_y, theta_z) eul2rotm([theta_z; theta_y; theta_x]);
    D = @(a, b, c) [1 / a^2, 0, 0; 0, 1 / b^2, 0; 0, 0, 1 / c^2];

    a_i = nthroot(omega_i * gamma_1_i^2 * gamma_2_i, 3);
    b_i = nthroot(omega_i * gamma_2_i / gamma_1_i, 3);
    c_i = nthroot(omega_i / gamma_2_i^2 / gamma_1_i, 3);
    D_i = D(a_i, b_i, c_i);
    R_i = R(theta_x_i, theta_y_i, theta_z_i);
    Q_i = R_i * D_i * R_i';

    a_j = nthroot(gamma_1_j^2 * gamma_2_j, 3);
    b_j = nthroot(gamma_2_j / gamma_1_j, 3);
    c_j = nthroot(1 / gamma_2_j^2 / gamma_1_j, 3);
    D_j = D(a_j, b_j, c_j);
    D_J = D_j^(-1/2);
    R_j = R(theta_x_j, theta_y_j, theta_z_j);

    Q_i_hat = D_J * R_j' * Q_i * R_j * D_J;
    [a_i_hat, b_i_hat, c_i_hat, theta_x_i_hat, theta_y_i_hat, theta_z_i_hat, D_i_hat, R_i_hat] = parametrize(Q_i_hat);
    
    a_i_bar = a_i_hat;
    b_i_bar = b_i_hat;
    c_i_bar = c_i_hat;
    theta_x_i_bar = 0;
    theta_y_i_bar = 0;
    theta_z_i_bar = 0;
    o_i_bar = [0; 0; 0];

    E_i_bar = [a_i_bar, b_i_bar, c_i_bar, theta_x_i_bar, theta_y_i_bar, theta_z_i_bar, o_i_bar(1), o_i_bar(2), o_i_bar(3)];

    a_j_bar = 1;
    b_j_bar = 1;
    c_j_bar = 1;
    theta_x_j_bar = 0;
    theta_y_j_bar = 0;
    theta_z_j_bar = 0;

    x_i_bar = [a_i_bar*cos(phi)*cos(psi); b_i_bar*sin(phi)*cos(psi); c_i_bar*sin(psi)];
%     n_i_bar = normalvec(E_i_bar, x_i_bar);
    n_i_bar = [b_i_bar*c_i_bar*cos(phi)*cos(psi)*sin(psi); c_i_bar*a_i_bar*sin(phi)*cos(psi)*sin(psi); a_i_bar*b_i_bar*sin(psi)^2];
    n_i_bar = n_i_bar/norm(n_i_bar);
    o_j_bar = x_i_bar + (1 + epsilon_bar)*n_i_bar;

    T1 = @(x_bar) R_i_hat * (x_bar - o_j_bar);
    T2 = @(x_hat) R_j * D_J * x_hat + o_j;
%     T0 = @(x) T2(T1(x));
    T  = @(x) R_j * D_J * R_i_hat * (x - o_j_bar) + o_j;

% %     o_i_hat = T1(o_i_bar);
% %     o_i = T2(o_i_hat);
%     o_i = T(o_i_bar);
    o_i = o_j - R_j * D_J * R_i_hat * o_j_bar;

%     x_i_hat = T1(x_i_bar);
%     x_i = T2(x_i_hat);
    x_i = T(x_i_bar);

    E_j_bar = [a_j_bar, b_j_bar, c_j_bar, theta_x_j_bar, theta_y_j_bar, theta_z_j_bar, o_j_bar(1), o_j_bar(2), o_j_bar(3)];
    E_i = [a_i, b_i, c_i, theta_x_i, theta_y_i, theta_z_i, o_i(1), o_i(2), o_i(3)];
    E_j = [a_j, b_j, c_j, theta_x_j, theta_y_j, theta_z_j, o_j(1), o_j(2), o_j(3)];

    [x_i_C, x_j_C, fval1, fval2, exitflag1, exitflag2, output1, output2] = Adapted_Choi_Solve(E_i, E_j);
    [x_i_C_bar, x_j_C_bar, fval1_hat, fval2_hat, exitflag1_hat, exitflag2_hat, output1_hat, output2_hat] = Adapted_Choi_Solve(E_i_bar, E_j_bar);

    x_i_C_hat = T1(x_i_C_bar);
    x_i_C_T = T2(x_i_C_hat);

    x_j_C_hat = T1(x_j_C_bar);
    x_j_C_T = T2(x_j_C_hat);



    varargout(1) = {fval1};
    varargout(2) = {fval2};
    varargout(3) = {exitflag1};
    varargout(4) = {exitflag2};
    varargout(5) = {output1};
    varargout(6) = {output2};
    varargout(7) = {fval1_hat};
    varargout(8) = {fval2_hat};
    varargout(9) = {exitflag1_hat};
    varargout(10) = {exitflag2_hat};
    varargout(11) = {output1_hat};
    varargout(12) = {output2_hat};

%% Calculs des separations

epsilon = sign(epsilon_bar)*norm(x_i_C - x_j_C);

end
