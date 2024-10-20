% -------------------------------------------------------------------------
% Distributed Wireless System Simulation
% This script simulates a distributed wireless communication system with N
% users and K subcarriers. The script generates random user distances,
% models path loss using a Rician fading channel model, and initializes
% several matrices for signal processing tasks. The outputs include path
% loss values (L1), channel coefficients (Hh), power allocation matrix (P),
% and weighting matrix (W).
%
% Variables:
%   N     - Number of users
%   K     - Number of subcarriers
%   T0_dB - Reference path loss in dB
%   T0    - Reference path loss in linear scale
%   aa    - Path loss exponent
%   e     - Rician factor (set to 0 for Rayleigh fading)
%   L1    - Path loss values for each user
%   Hh    - Channel coefficients (LOS and NLOS components)
%   W     - Weighting matrix for each user on each subcarrier
%   P     - Power allocation matrix for each user on each subcarrier
%   LL    - Placeholder matrix for future calculations (currently initialized to zero)
% -------------------------------------------------------------------------

clear;
 
% Define system parameters
N = 10;       % Number of users
K = 64;       % Number of subcarriers
numIterations = 100;  % Number of iterations for the simulation

% Preallocate matrices for efficiency
L1 = zeros(numIterations, N, K);  % Path loss values
Hh = zeros(numIterations, N, K);  % Channel coefficients
W = zeros(numIterations, N, K);   % Weighting matrix
P = zeros(numIterations, N, K);   % Power allocation matrix
LL = zeros(numIterations, N, K);  % Placeholder for future use

% Constants for path loss and Rician fading model
T0_dB = -25;    % Reference path loss in dB
T0 = 10^(T0_dB / 10);  % Convert dB to linear scale
aa = 2.2;       % Path loss exponent
e = 0;          % Rician factor (0 for Rayleigh fading)
e1 = e / (1 + e);   % Scaling factor for LOS component
e2 = 1 / (1 + e);   % Scaling factor for NLOS component

% Simulation loop for 'numIterations' iterations
for counter = 1:numIterations
    % Generate random distances for each user (uniformly between 10 and 20 units)
    dist = 10 + 10 * rand(1, N);
    
    % Calculate path loss for each user using the path loss model
    L1(counter, :, :) = sqrt(T0 * dist.^(-aa));

    % Generate channel coefficients for Rician fading model
    % LOS component
    users_angle = rand(N, 1);  % Random angle for each user
    h_LOS = exp(1j * pi * sin(users_angle) .* (0:K-1));  % LOS channel
    
    % NLOS component (Rayleigh fading, complex Gaussian noise)
    h_NLOS = sqrt(1/2) * (randn(N, K) + 1j * randn(N, K));  
    
    % Combined channel matrix (Hh) with LOS and NLOS components
    Hh(counter, :, :) = (e1 * h_LOS + e2 * h_NLOS);
    
    % Initialize weighting matrix W (randomly distributed)
    W(counter, :, :) = 0.1 * rand(N, K) / N;
    
    % Initialize power allocation matrix P (randomly distributed)
    P(counter, :, :) = 0.1 * rand(N, K) / (N * K);
    
    % Placeholder matrix LL (currently set to zero)
    LL(counter, :, :) = 0 * rand(1, N);
end
