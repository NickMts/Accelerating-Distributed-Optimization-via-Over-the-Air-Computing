% -------------------------------------------------------------------------
% Distributed OTA Optimization for Power Grid Management with PEVGs
%
% This script performs a distributed optimization for power allocation 
% among Plug-in Electric Vehicle Groups (PEVGs) in a power grid. The goal 
% is to optimize the load distribution based on battery capacities, satisfaction 
% parameters, and power prices, while satisfying constraints such as grid 
% capacity.
%
% Variables:
%   N        - Number of PEVGs (users)
%   fj       - Objective function value
%   con      - Constraint violations
%   pp       - Power price array
%   zeta, theta - Constants for projection update
%   alpha, beta - Step-size constants for the optimization process
%   precoding - Precoding factor for transmissions
%   lambda   - Lagrange multipliers
%   x, y     - Primal variables (x: power allocation, y: load satisfaction)
%   b        - Battery capacity
%   s        - Satisfaction parameter
%   participate - Tracks PEVG participation
%
% The script also uses a projection algorithm to ensure that the power 
% allocation and satisfaction values remain within the specified constraints.
% -------------------------------------------------------------------------

clear;

% System parameters
N = 20;  % Number of PEVGs (users)
fj = 0;  % Objective function value accumulator
conj = 0;  % Constraint violations accumulator
targ = 0;  % Target objective value accumulator
N0 = -90;  % Noise parameter
counter = 5000;  % Number of iterations for optimization
alpha = 1;  % Step-size constant (alpha)
beta = 2;   % Step-size constant (beta)
participate = zeros(1, counter);  % Tracks participation
pp = zeros(1, 20);  % Price per unit accumulator
zeta = 2;  % Projection constant (zeta)
theta = 2;  % Projection constant (theta)
precoding = 1 * 10^6;  % Precoding factor
C = 99;  % Grid capacity (maximum total load)

% Main optimization loop for 20 iterations
for j = 1:20
    j  % Display the current iteration number

    % Initialize system parameters for each iteration
    y = 100 * ones(1, N);  % Initial load satisfaction (constant for each PEVG)
    b = unifrnd(35, 65, 1, N);  % Random battery capacity for each PEVG
    s = 1 + rand(1, N);  % Random satisfaction parameter [1, 2]
    p = 17;  % Initial price per unit of power
    ymax = 0.5 * (b - p).^2 ./ s;  % Maximum load satisfaction
    
    % Initial values for primal-dual variables
    x = rand(1, N);  % Initial power allocation
    lambda = 5 * rand(1, N);  % Initial Lagrange multipliers (random)

    f_0 = zeros(1, counter);  % Store objective function value for each iteration
    X = zeros(N, counter);  % Store power allocation over iterations
    Y = zeros(N, counter);  % Store load satisfaction over iterations
    Lambda = zeros(N, counter);  % Store Lagrange multipliers over iterations
    a = zeros(1, counter);  % Store step sizes

    for i = 1:counter
        % Compute the step size for this iteration
        a(i) = alpha / (beta + i);
        Dmax = zeta + theta * sqrt(sum(a(1:end-1)));  % Projection parameter
        D(i) = Dmax;
        xprev = x;  % Store previous power allocation

        % Update primal variable x (power allocation)
        vector = [lambda .* (0 - (b - s .* x - p)), lambda .* (1 * ones(1, N) - 0)];

        % Path loss model for communication
        dist = 10 + 10 * rand(1, N);  % Random distances for each PEVG
        T0_dB = -25;  % Path loss in dB
        T0 = 10^(T0_dB / 10);  % Convert to linear scale
        aa = 2;  % Path loss exponent
        L1 = sqrt(T0 * dist.^(-aa));  % Path loss model
        
        % Rician factor and channel generation
        e = 10;  % Rician factor
        e1 = e / (1 + e);
        e2 = 1 / (1 + e);
        users_angle = rand(N, 1);  % Random angle for each PEVG
        h_LOS = exp(1j * pi * sin(users_angle) .* (0:K-1));  % LOS channel component
        h_NLOS = sqrt(1 / 2) * (randn(N, K) + 1j * randn(N, K));  % NLOS component
        H = (e1 * h_LOS + e2 * h_NLOS);  % Total channel gain

        % Update primal variable x with Lagrangian multipliers
        mltx = lambda .* (0 - (b - s .* x - p));
        mlty = lambda .* (1 * ones(1, N) - 0);

        % Communication participation check using channel conditions
        cntr = 0;  % Counter for participants
        for jj = 1:N
            mult = norm([vector(jj), vector(jj + N)])^2 / precoding * snr / Pmax;
            if abs(H(jj))^2 < mult
                mltx(jj) = 0;  % No participation for this PEVG
                mlty(jj) = 0;
                cntr = cntr + 1;
            end
        end
        participate(i) = participate(i) + cntr;  % Track participation
        
        % Update primal variables x and y
        mltx = mltx + sqrt(precoding) / snr * wgn(1, N, N0);
        mlty = mlty + sqrt(precoding) / snr * wgn(1, N, N0);
        x = x - a(i) * (mltx);

        % Projection of x (ensure power allocation constraint)
        flag = 1;
        while flag ~= 0
            x = max(x, 0);  % Ensure non-negative values
            C_temp = sum(x);  % Check if the total power exceeds capacity
            if C_temp > C
                len = length(x(x > 0));
                x(x > 0) = x(x > 0) - (C_temp - C) / len;
            end
            if length(x(x >= 0)) == N && sum(x) <= C
                flag = 0;  % Stop projection when condition is satisfied
            end
        end

        % Update y (load satisfaction) with projection
        yprev = y;
        y = max(0, y - a(i) * (-1 * ones(1, N) + mlty));  % Update load satisfaction
        y = min(y, ymax);  % Ensure satisfaction doesn't exceed ymax

        % Update Lagrange multipliers lambda
        lambda = lambda + a(i) * (yprev - (b .* x - 0.5 * s .* x.^2 - p * x));
        lambda = max(0, lambda);  % Ensure non-negative lambda
        lambda = min(lambda, Dmax);  % Apply projection to lambda

        % Store results for plotting
        X(:, i) = x;
        Y(:, i) = y;
        Lambda(:, i) = lambda;
        f_0(i) = sum(y);
        con(i) = norm(max(-(b .* x - 0.5 * s .* x.^2 - p * x - y), 0)) / N;
        target(i) = sum(b .* x - 0.5 * s .* x.^2 - p * x);
    end

    % Accumulate results over all iterations
    fj = f_0 + fj;
    targ = targ + target;
    conj = conj + con;
end

% Normalize results by number of iterations
fj = fj / j;
targ = targ / j;
conj = conj / j;
participate = (50 - participate / counter) / 50;

% Plot the results
figure;
plot(1:counter, fj, 1:counter, targ);
xlabel('Iteration');
ylabel('Objective Value');
title('Objective Function and Target over Iterations');

figure;
semilogy(1:counter, conj);
xlabel('Iteration');
ylabel('Constraint Violation');
title('Constraint Violations over Iterations');
