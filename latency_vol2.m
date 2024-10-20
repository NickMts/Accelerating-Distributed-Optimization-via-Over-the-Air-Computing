% -------------------------------------------------------------------------
% Simulation of OMA and OTA Performance in a Wireless Communication System
% 
% This script simulates and compares the performance of Orthogonal Multiple
% Access (OMA) and Over-the-Air (OTA) computing in a wireless communication
% system with N users. The simulation calculates the latency (time delay)
% for data transmission using both OMA and OTA approaches, using a path loss
% model and a Rician fading channel model.
%
% Variables:
%   D     - Total data size to transmit (bits)
%   b     - Number of bits per transmission block
%   B     - Total bandwidth (Hz)
%   N_0   - Noise power spectral density (W/Hz)
%   N     - Number of users
%   K     - Number of subcarriers (set to 1 for single-carrier systems)
%   P     - Power per subcarrier
%   iter  - Number of iterations for the latency calculation
%   oma_l - Latency values for OMA
%   ota_l - Latency values for OTA
%   fdma_l - Latency for FDMA (currently unused)
%   H     - Channel coefficients (LOS and NLOS components)
% -------------------------------------------------------------------------

clear;

% Define system parameters
D = 5 * 2 * 640;  % Data size (bits)
b = 17;           % Bits per transmission block
B = 10^6;         % Total bandwidth (Hz)
N_0 = 4 * 10^(-21); % Noise power spectral density (W/Hz)
N = 50;           % Number of users
K = 1;            % Number of subcarriers
B = B / K;        % Bandwidth per subcarrier
P = 1 / K;        % Power per subcarrier
iter = 100;       % Number of iterations for the simulation

% Initialize latency arrays
oma_l = zeros(1, 1000);  % Latency for OMA
ota_l = zeros(1, 1000);  % Latency for OTA
fdma_l = 0;              % Placeholder for FDMA latency (unused)

% Simulation loop for 'iter' iterations
for i = 1:100
    for j = 1:iter
        % Path loss model
        dist = 10 + 10 * rand(1, N);  % Random distance for each user (10 to 20 units)
        T0_dB = -25;                  % Reference path loss in dB
        T0 = 10^(T0_dB / 10);         % Convert to linear scale
        aa = 2.2;                     % Path loss exponent
        L1 = sqrt(T0 * dist.^(-aa))'; % Path loss for each user
        
        % Rician fading model
        e = 10;             % Rician factor
        e1 = e / (1 + e);   % Scaling factor for LOS component
        e2 = 1 / (1 + e);   % Scaling factor for NLOS component
        
        % Generate channel coefficients
        users_angle = rand(N, 1);                          % Random angle for each user
        h_LOS = exp(1j * pi * sin(users_angle) .* (0:K-1)); % LOS component
        h_NLOS = sqrt(1/2) * (randn(N, K) + 1j * randn(N, K));  % NLOS component
        
        % Combined channel matrix (H) with path loss and Rician fading
        H = (abs(e1 * h_LOS + e2 * h_NLOS) .* L1).^2;  % Channel gain
        
        % Calculate latency for OMA
        oma_l(j) = oma_l(j) + (b * D + 64) * sum(1 ./ (B * log2(1 + P .* H ./ (B * N_0))));
    end
end

% Average OMA latency over all iterations
oma_l = oma_l / 100;

% OTA latency (constant)
B = 10^6;  % Reset bandwidth to initial value
ota_l = D / B * ones(1, 100);  % OTA latency is constant across iterations
ota = ota_l(1) * (1:100);  % Total OTA latency over iterations

% Calculate cumulative OMA latency
oma = zeros(1, 100);
for i = 1:100
    oma(i) = sum(oma_l(1:i));  % Cumulative OMA latency
end

% Plot results
plot(1:100, oma, 1:100, ota);
xlabel('Iteration');
ylabel('Cumulative Latency');
legend('OMA Latency', 'OTA Latency');
title('Cumulative Latency: OMA vs OTA');
grid on;
