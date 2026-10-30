%% identification.m ------------------------------------------------------
% This script contains all the necessary steps for the identification of
% the system. The work displayed below conveys the scripts necessary for
% Labs L4 - L7


% Data processing, ARMAX model identification
% from plant data, and conversion to state model is performed in one file (of type .m).
% This file starts by reading plant data stored in a .mat file and produces as output
% the matrices that define the state model that become available in the MATLAB
% Workspace.

clear; close all; clc;

%% Constants -------------------------------------------------------------
dataset = {'L04_1v_good_v.mat', 'L04_1.25v.mat'}; % Datasets to verify the model   
Ts = 0.02;        
col_u = 1;        
col_theta = 2;    
col_alpha = 3;

% Gains defined in L2/L3
Ke = 4.57;
Kp = -36;

%% L4 ---------------------------------------------------------------------


% Cicle for the graphs
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

%% L5 ---------------------------------------------------------------------
t_dis = 10; % The first 10s are to be avoided (transitory)
af = 0.8; % Pole of the lower-pass filter
