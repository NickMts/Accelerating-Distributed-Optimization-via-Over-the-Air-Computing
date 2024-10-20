% -------------------------------------------------------------------------
% Distributed Optimization for Power Allocation Among PEVGs
%
% This script performs a distributed optimization for power allocation among
% Plug-in Electric Vehicle Groups (PEVGs) in a power grid. The optimization
% focuses on balancing power allocation and satisfaction, while considering
% constraints such as grid capacity and the battery capacity of each PEVG.
%
% Variables:
%   N        - Number of PEVGs (users)
%   fj       - Objective function value accumulator
%   conj     - Constraint violations accumulator
%   targ     - Target objective value accumulator
%   b        - Battery capacity of each PEVG
%   s        - Satisfaction parameter for each PEVG
%   p        - Price per unit of power
%   C        - Grid capacity (maximum total load)
%   ymax     - Maximum load satisfaction for each PEVG
%   x        - Power allocation for each PEVG
%   y        - Load satisfaction for each PEVG
%   lambda   - Lagrange multipliers for optimization
%   a        - Step size for optimization
%
% The script uses a projection algorithm to ensure power allocation (x) and 
% satisfaction (y) remain within the constraints, ensuring efficient power
% distribution.
% -------------------------------------------------------------------------

clear;

% System parameters
N = 20;  % Number of PEVGs (users)
fj = 0;  % Objective function value accumulator
conj = 0;  % Constraint violations accumulator
targ = 0;  % Target objective value accumulator

% Main optimization loop for 20 iterations
for j = 1:20
    j  % Display the current iteration number

    % Initialize system parameters for each iteration
    b = unifrnd(35, 65, 1, N);  % Random battery capacity for each PEVG
    s = 1 + rand(1, N);  % Random satisfaction parameter [1, 2]
    p = 17;  % Price per unit of power
    C = 99;  % Grid capacity (maximum total load)
    ymax = 0.5 * (b - p).^2 ./ s;  % Maximum load satisfaction for each PEVG

    % Initialize primal-dual variables
    x = rand(1, N);  % Initial power allocation
    y = 100 * ones(1, N);  % Initial load satisfaction (constant for each PEVG)
    lambda = rand(1, N);  % Initial Lagrange multipliers

    % Parameters for the optimization
    counter = 5000;  % Number of iterations for optimization
    alpha = 1;  % Step size constant
    beta = 2;  % Step size constant
    zeta = 2;  % Projection constant
    theta = 2;  % Projection constant

    % Storage for results
    f_0 = zeros(1, counter);  % Objective function value for each iteration
    X = zeros(N, counter);  % Store power allocation over iterations
    Y = zeros(N, counter);  % Store load satisfaction over iterations
    Lambda = zeros(N, counter);  % Store Lagrange multipliers over iterations
    a = zeros(1, counter);  % Store step sizes

    % Iterative optimization process
    for i = 1:counter
        i;  % Display the current iteration number
        
        % Step size calculation
        a(i) = alpha / (beta + i);
        Dmax = zeta + theta * sqrt(sum(a(1:end-1)));  % Projection parameter
        D(i) = Dmax;
        
        % Update primal variable x (power allocation)
        xprev = x;  % Store previous power allocation
        x = x - a(i) * (zeros(1, N) + lambda .* (0 - (b - s .* x - p)));  % Gradient step
        
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
        
        % Update primal variable y (load satisfaction) with projection
        yprev = y;  % Store previous load satisfaction
        y = max(0, y - a(i) * (-1 * ones(1, N) + lambda .* (1 * ones(1, N) - 0)));  % Gradient step
        y = min(y, ymax);  % Ensure satisfaction doesn't exceed ymax

        % Store results for plotting
        X(:, i) = x;
        Y(:, i) = y;
        
        % Update Lagrange multipliers lambda
        lambda = lambda + a(i) * (yprev - (b .* xprev - 0.5 * s .* xprev.^2 - p * xprev));
        
        % Projection of lambda
        lambda = max(0, lambda);  % Ensure non-negative lambda
        lambda = min(lambda, Dmax);  % Apply projection to lambda
        Lambda(:, i) = lambda;

        % Compute objective function and constraints
        f_0(i) = sum(y);  % Objective function: sum of load satisfaction
        con(i) = norm(max(-(b .* x - 0.5 * s .* x.^2 - p * x - y), 0)) / N;  % Constraint violation
        target(i) = sum(b .* x - 0.5 * s .* x.^2 - p * x);  % Target objective value
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
