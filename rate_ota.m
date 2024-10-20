% -------------------------------------------------------------------------
% Distributed Resource Allocation using Over-the-Air Computation
% 
% This script simulates a distributed optimization framework for resource
% allocation in a wireless communication system with N users and K subcarriers.
% The goal is to maximize the sum rate of users while satisfying power constraints
% using iterative optimization with over-the-air (OTA) computation.
% The optimization includes power and bandwidth allocation, and uses projections
% to ensure the constraints are met for each iteration.
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
%   snr      - Signal-to-noise ratio
%   precoding - Precoding factor for transmissions
%   Pmax     - Maximum allowed power per user
%   obj      - Array to store objective values for each iteration
%   con      - Array to store constraint violations for each iteration
%   participate - Tracks user participation in each iteration
%
% The algorithm includes a projection step to ensure the power (p) and bandwidth (w)
% constraints are met after each iteration.
% -------------------------------------------------------------------------

clear;

% System parameters
N = 10;        % Number of users
K = 64;        % Number of subcarriers
Pofdm = 1;     % Total OFDM power
N0 = 3.9811 * 10^(-21) / Pofdm; % Noise power spectral density
B = 10^6;      % Total bandwidth (Hz)
alpha = 1; beta = 80000; % Step-size constants
Rth = 2.5 * ones(1, N);  % Minimum rate threshold for each user
zeta = 1; theta = 1;     % Projection constants
counter = 100;           % Number of iterations for the optimization
obj = zeros(1, counter); % Objective value for each iteration
con = zeros(1, counter); % Constraint violation for each iteration
snr = 10^0;              % Signal-to-noise ratio
precoding = 1 * 10^4;    % Precoding factor
Pmax = 1;                % Maximum power per user
iter = 100;              % Number of outer iterations
participate = zeros(1, counter); % Tracks participation

% Iterative optimization loop
for j = 1:iter
    % Generate initial random channel states and power/bandwidth allocations
    H = reshape(Hh(j, :, :), [N, K]); % Reshape channel matrix H
    L11 = reshape(L1(j, :, :), [1, N]); % Reshape path loss values
    
    % Adjust channel matrix with path loss
    for i = 1:N
        H(i, :) = abs(H(i, :) * L11(i)).^2;
    end
    
    % Initialize variables for each user
    w = reshape(W(j, :, :), [N, K]);  % Weighting factors (bandwidth allocation)
    p = reshape(P(j, :, :), [N, K]);  % Power allocation
    lambda = reshape(LL(j, :, :), [1, N]); % Lagrange multipliers
    
    % Iterative optimization
    for i = 1:counter
        % Step-size calculation
        a = alpha / (beta + i);
        D = zeta + theta * sqrt(sum(a(1:end-1))); % Projection parameter

        % Store previous values of w
        wprev = w;
        L = zeros(N, K); % Initialize variable for Lagrangian update
        
        % Compute Lagrangian multiplier updates over the air
        for jj = 1:N
            L(jj, :) = (-1 - lambda(jj)) / (K * 1);
        end
        
        % Projection and participation check
        vector = (1/(K*1) + L) .* (log(1 + H .* p ./ (w * N0 * B / K)) / log(2) - ...
            H .* p ./ (N0 * B / K * log(2) * w .* (1 + H .* p ./ (w * N0 * B / K)))); 
        cntr = 0;
        
        % Participation and projection check for each user
        for jj = 1:N
            mult = norm([vector(jj, :)])^2 / (precoding * Pmax);
            if abs(H(jj, :)) < mult
                cntr = cntr + 1;
                L(jj, :) = -1 / (K * 1); % If not participating, reset Lagrange multiplier
            end
        end
        
        participate(i) = cntr + participate(i); % Track participation
        
        % Update w and p based on the Lagrangian multiplier
        w = w - a * (L .* (log(1 + H .* p ./ (w * N0 * B / K)) / log(2) - ...
            H .* p ./ (N0 * B / K * log(2) * w .* (1 + H .* p ./ (w * N0 * B / K)))) + ...
            sqrt(precoding) / snr * wgn(N, K, -90));
        p = p - a * (L .* (H ./ (N0 * B / K * log(2) * (1 + H .* p ./ (wprev * N0 * B / K)))) + ...
            sqrt(precoding) / snr * wgn(N, K, -90));
        
        % --- Projection for variable w ---
        for jj = 1:K
            ww = w(:, jj); % Select the column for projection
            flag = 1;
            while flag ~= 0
                ww = max(ww, 10^(-15)); % Ensure non-negative values
                C_temp = sum(ww); % Calculate sum of weights
                if abs(C_temp - 1) > 10^(-14) % Projection step
                    len = length(ww(ww >= 0)); 
                    ww(ww >= 0) = ww(ww >= 0) - (C_temp - 1) / len;
                end
                if length(ww(ww >= 0)) == N && abs(sum(ww) - 1) <= 10^(-14)
                    flag = 0; % Stop when the projection condition is satisfied
                end
            end
            w(:, jj) = ww; % Update w after projection
        end
        
        % --- Projection for variable p ---
        pp = reshape(p', 1, []); % Flatten power matrix for easier processing
        flag = 1;
        while flag ~= 0
            pp = max(pp, 10^(-15)); % Ensure non-negative values
            C_temp = sum(pp); % Calculate total power
            if C_temp > Pofdm % If power exceeds the maximum, project
                len = length(pp(pp > 0));
                pp(pp > 0) = pp(pp > 0) - (C_temp - Pofdm) / len;
            end
            if length(pp(pp >= 0)) == N * K && sum(pp) <= Pofdm
                flag = 0; % Stop when projection condition is satisfied
            end
        end
        p = reshape(pp, K, N)'; % Reshape back into the matrix form
        
        % Update lambda values based on rate constraints
        bb = 0;
        for jj = 1:N
            rate_diff = Rth(jj) - sum(1 / K * w(jj, :) .* log2(1 + p(jj, :) .* H(jj, :) ./ (w(jj, :) * N0 * B / K)));
            lambda(jj) = lambda(jj) + a * rate_diff;
            bb = max(bb, max(0, rate_diff)); % Track the maximum constraint violation
        end
        lambda = max(0, lambda); % Ensure non-negative lambda
        lambda = min(lambda, D); % Projection for lambda
        
        % Store objective function and constraint violation
        obj(i) = obj(i) + sum(sum(1 / K * w .* log2(1 + p .* H ./ (w * N0 * B / K))));
        con(i) = con(i) + bb;
    end
end

% Average the results over iterations
obj = obj / iter;
con = con / iter;
participate = (N - participate / iter) / N;

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
