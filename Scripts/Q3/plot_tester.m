%% identification.m ------------------------------------------------------
% This script contains all the necessary steps for the identification of
% the ARMAX model, such as the processing of the plant data recording in
% .mat files, followed by an iteration of the polynomial orders of the
% model and finally the identification of the system itself. Lastly, the
% state-space matrices are calculated.

clear; close all; clc;

%% Collection of identification data --------------------------------------
% Datasets folder (one level up from this script, in ../Datasets)
folder = fileparts(mfilename('fullpath'));
addpath(fullfile(folder, '..', 'Datasets'));
% Datasets to verify the model
dataset = {'L04_1v_good_v.mat', 'L04_1.25v.mat', 'PRBS_0.8v_0.4f.mat', 'L04_0.5v.mat', 'SW_0.2v_0.4f.mat', 'KP_LAB3_4.mat'};  
Ts = 0.02;        
col_u = 1;        
col_theta = 2;    
col_alpha = 3;

% Gains defined in L2/L3
Ke = 4.57;
Kp = -36;

% Filter parameter
af = 0.8;

%%  Verification ---------------------------------------------------------------------
% Cycle for the graphs
for k = 1:numel(dataset)
    S = load(dataset{k});
    X = S.input;                      
    t = X.time(:);
    sigs = X.signals.values;

    % Graphs
    figure('Name', dataset{k});
    subplot(3,1,1); plot(t, sigs(:, col_u)); grid on;
    title(dataset{k}, 'Interpreter', 'none');
    xlabel('Time [s]'); ylabel('u [V]');
    subplot(3,1,2); plot(t, sigs(:, col_theta)); grid on;
    xlabel('Time [s]'); ylabel('\theta_e [V]');
    subplot(3,1,3); plot(t, sigs(:, col_alpha)); grid on;
    xlabel('Time [s]'); ylabel('\alpha_e [V]');

end
