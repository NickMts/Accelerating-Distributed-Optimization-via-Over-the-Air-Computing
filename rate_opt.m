% -------------------------------------------------------------------------
% Distributed Resource Allocation with Projection Algorithm
% 
% This script performs distributed optimization for resource allocation 
% in a wireless communication system with N users and K subcarriers. The 
% optimization includes power and bandwidth allocation across subcarriers 
% while ensuring that the constraints are satisfied.
%
% The projection algorithm ensures that the power (p) and bandwidth (w) 
% constraints are met after each iteration.
%
% Variables:
%   N        - Number of users
%   K        - Number of subcarriers
%   Pofdm    - Total available power for OFDM
%   N0       - Noise power spectral density
%   B        - Total bandwidth (Hz)
%   alpha, beta - Constants for step-size adaptation
%   Rth      - Minimum rate threshold for each user
%   zeta, theta - Constants for projection update
%   counter  - Number of iterations for the optimization
%   obj      - Array to store objective values for each iteration
%   con      - Array to store constraint violations for each iteration
%   iter     - Number of outer iterations
%
% Projection Steps:
% - Projection for `w` ensures that the bandwidth allocation across each 
%   subcarrier sums to 1.
% - Projection for `p` ensures that the total power allocation does not 
%   exceed Pofdm.
% -------------------------------------------------------------------------

clear;

% System parameters
N = 10;          % Number of users
K = 64;          % Number of subcarriers
Pofdm = 1;       % Total OFDM power
N0 = 3.9811 * 10^(-21); % Noise power spectral density
B = 10^6;        % Total bandwidth (Hz)
alpha = 1;       % Constant for step size
beta = 1000;     % Constant for step size
Rth = 2.5 * ones(1, N);  % Minimum rate threshold for each user
zeta = 2;        % Projection constant
theta = 1;       % Projection constant
counter = 100;   % Number of iterations for the optimization
obj = zeros(1, counter);  % Array to store objective values
con = zeros(1, counter);  % Array to store constraint violations
iter = 100;      % Number of outer iterations

% Main optimization loop
for j = 1:iter
    % Channel and path loss models
    H = reshape(Hh(j, :, :), [N, K]); % Channel matrix reshaped from input
    for i = 1:N
        H(i, :) = abs(H(i, :) * L1(i)).^2; % Apply path loss
    end
    
    % Initialize variables for each iteration
    w = reshape(W(j, :, :), [N, K]);  % Weighting factors (bandwidth allocation)
    p = reshape(P(j, :, :), [N, K]);  % Power allocation
    lambda = reshape(LL(j, :, :), [1, N]); % Lagrange multipliers
    a = 0; % Step size
    
    % Iterative optimization process
    for i = 1:counter
        [i, j]  % Display current iteration
        a(i) = alpha / (beta + i);  % Step-size adaptation
        D(i) = zeta + theta * sqrt(sum(a(1:end-1))); % Update projection parameter
        wprev = w;  % Store previous bandwidth allocation

        L = zeros(N, K);  % Initialize Lagrange multiplier matrix
        for jj = 1:N
            L(jj, :) = (-1 - lambda(jj)) / (K * 100);  % Update Lagrange multipliers
        end
        
        % Update bandwidth allocation w and power allocation p
        w = w - a(i) * L .* (log(1 + H .* p ./ (w * N0 * B / K)) / log(2) - ...
            H .* p ./ (N0 * B / K * log(2) * w .* (1 + H .* p ./ (w * N0 * B / K))));
        p = p - a(i) * L .* (H ./ (N0 * B / K * log(2) * (1 + H .* p ./ (wprev * N0 * B / K))));
        
        % --- Projection for variable w ---
        for jj = 1:K
            ww = w(:, jj);  % Select column for projection
            flag = 1;
            while flag ~= 0
                ww = max(ww, 10^(-15));  % Ensure non-negative values
                C_temp = sum(ww);  % Sum of weights for projection
                if abs(C_temp - 1) > 10^(-14)
                    len = length(ww(ww >= 0));
                    ww(ww >= 0) = ww(ww >= 0) - (C_temp - 1) / len;
                end
                if length(ww(ww >= 0)) == N && abs(sum(ww) - 1) <= 10^(-14)
                    flag = 0;  % Stop projection when condition is satisfied
                end
            end
            w(:, jj) = ww;  % Update bandwidth allocation after projection
        end
        
        % --- Projection for variable p ---
        pp = [];  % Initialize flat power allocation array
        for jj = 1:N
            pp = [pp, p(jj, :)];  % Flatten power allocation for projection
        end
        flag = 1;
        while flag ~= 0
            pp = max(pp, 10^(-15));  % Ensure non-negative values
            C_temp = sum(pp);  % Sum of power allocation
            if C_temp > Pofdm  % If total power exceeds limit
                len = length(pp(pp > 0));
                pp(pp > 0) = pp(pp > 0) - (C_temp - Pofdm) / len;
            end
            if length(pp(pp >= 0)) == N * K && sum(pp) <= Pofdm
                flag = 0;  % Stop projection when condition is satisfied
            end
        end
        for jj = 1:N
            p(jj, :) = pp((jj - 1) * K + 1:jj * K);  % Reshape power allocation
        end
        
        % Update Lagrange multipliers lambda based on rate constraints
        bb = 0;
        for jj = 1:N
            rate_diff = Rth(jj) - sum(1 / K * w(jj, :) .* log2(1 + p(jj, :) .* H(jj, :) ./ (w(jj, :) * N0 * B / K)));
            lambda(jj) = lambda(jj) + a(i) * rate_diff;  % Update lambda
            bb = max(bb, max(0, rate_diff));  % Track maximum constraint violation
        end
        lambda = max(0, lambda);  % Ensure non-negative lambda
        lambda = min(lambda, D(i));  % Apply projection to lambda
        
        % Store objective and constraint values
        obj(i) = obj(i) + sum(sum(1 / K * w .* log2(1 + p .* H ./ (w * N0 * B / K))));
        con(i) = con(i) + bb;
    end
end

% Average results over iterations
obj = obj / iter;
con = con / iter;

% Plot the objective function and constraint violation over iterations
figure;
plot(1:counter, obj);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function Over Iterations');

figure;
plot(1:counter, con);
xlabel('Iteration');
ylabel('Constraint Violation');
title('Constraint Violation Over Iterations');
