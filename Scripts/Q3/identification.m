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
dataset = {'L04_1v_good_v.mat', 'L04_1.25v.mat'};  
Ts = 0.02;        
col_u = 1;        
col_theta = 2;    
col_alpha = 3;

% Gains defined in L2/L3
Ke = 4.57;
Kp = -36;

%%  ---------------------------------------------------------------------
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

%% Processing of the ARMAX model ------------------------------------------
t_dis = 10; % The first 10s are to be avoided (transitory)
file = 'L04_1v_good_v.mat';

% Loading Data
S = load(file);
X = S.input;
t = X.time(:);
sigs = X.signals.values;

utrend = sigs(:, col_u);        
thetae = sigs(:, col_theta);    
alphae = sigs(:, col_alpha); 

% Calculation of ytrend
jumps = [0; diff(thetae)]; % Notices the 10V jumps at (-5V -> 5V)
thetae = thetae - 10 * cumsum(round(jumps / 10)); % Makes thetae continuous
ytrend = thetae * Kp + alphae * Ke;

% Filter Parameters
af = 0.8; % Pole of the lower-pass filter
Afilt = [1 -af];
Bfilt = (1 - af) * [1 -1];
yf = filter(Bfilt, Afilt, ytrend);       

% Choice of the model orders (iterated based on the printed table)
na = 2; % AR part
nb = 4; % X part
nc = 2 ; % MA part
nk = 1; % pure delay

% Signal Processing
idx = t >= t(1) + t_dis; % Mask for discarding the first t_dis [s]
tt = t(idx);
yf = yf(idx);
u = utrend(idx) - mean(utrend(idx));

% Figure
figure('Name', 'Processed Data');
subplot(2,1,1); plot(tt, u); grid on;
xlabel('Time [s]'); ylabel('u (detrended) [V]');
subplot(2,1,2); plot(tt, yf); grid on;
xlabel('Time [s]'); ylabel('differentiated, filtered y [deg]');

% Model fitting function
fit = @(y, ys) 100 * (1 - norm(y - ys) / norm(y - mean(y))); % Definition of the fitting function
z = iddata(yf, u, Ts); % Necessary to save Ts         
opt = armaxOptions('Display', 'off');

% Sweep cycle
fprintf('\n  na  nb  nc | fit [%%]\n');
fprintf('--------------------\n');
res = [];
for na_i = 1:4
    for nb_i = 1:4
        th_i = armax(z, [na_i nb_i na_i nk], opt); % nc = na
        [den_i, num_i] = polydata(th_i);
        ys_i = filter(num_i, den_i, u); % Getting the value 
        fitting = fit(yf, ys_i);
        if fitting > 0 && fitting < 100
            fprintf(' %3d %3d %3d | %6.1f\n', na_i, nb_i, na_i, fitting);
        else
            fprintf(' %3d %3d %3d | %6.1f\n', na_i, nb_i, na_i, NaN);
        end
    end
end

% Creation of the definitive ARMAX model
th = armax(z, [na nb nc nk], opt);
fprintf("\n");
disp(th);
[den_armax , num_armax] = polydata(th);        

ys = filter(num_armax, den_armax, u);
fprintf('\nFinal model fitting (na=%d nb=%d nc=%d): %.1f %%\n', na, nb, nc, fit(yf, ys));

% Database plot
figure('Name', 'Model vs Data');
plot(tt, yf, tt, ys); grid on;
xlabel('Time [s]'); ylabel('differentiated angle [deg]');
legend('data', 'model');

% Zoomed Data plot
figure('Name', 'Model vs Data');
plot(tt, yf, tt, ys); grid on;
xlabel('Time [s]'); ylabel('differentiated angle [deg]');
xlim([15 17]);
legend('data', 'model');

%% State-space model -----------------------------------------------------
% Adding the integrator
den_int = conv(den_armax, [1 -1]);
[num,den] = eqtflength (num_armax, den_int);

% Calculating the state-space matrices
[A, B, C, D] = tf2ss(num, den);   % x(k+1)=Ax+Bu ; y=Cx+Du
sys = ss(A, B, C, D, Ts);

% Printing the state-space matrices
fprintf('\nState-space matrices (n = %d):\n', size(A, 1))
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);
disp('D ='); disp(D);

% Confirming the asymptotic stability of the system with and without the 
% integrator
fprintf("ARMAX poles:\n");
disp(abs(roots(den_armax)));
fprintf("State-space poles:\n");
disp(abs(roots(den)));

%% Validation of the state-space model ------------------------------------
% Processing of the validation data-set
vds = 'L04_1.25v.mat';
V = load(vds);
VI = V.input;
t = VI.time(:);
sigs_val = VI.signals.values;
u_val = sigs_val(:, col_u);
thetae_val  = sigs_val(:, col_theta);
alphae_val  = sigs_val(:, col_alpha);

jumps_val = [0; diff(thetae)];
thetae = thetae_val - 10 * cumsum(round(jumps_val / 10));   % tira saltos do pot.
ytrend_val = thetae_val * Kp + alphae * Ke;                  % angulo total [graus]

yf_val = filter(Bfilt, Afilt, ytrend);                   % diferenciar + filtra
yf_v = yf_val(idx);
u_v  = u_val(idx) - mean(u_val(idx));

%% 4. VALIDACAO (a): modelo ARMAX vs dados NAO usados na identificacao
ys_v = filter(num_armax, den_armax, u_v);
fprintf('\nFIT DE VALIDACAO: %.1f %%\n', fit(yf_v, ys_v));

figure('Name', 'Validacao');
plot(tt, yf_v, tt, ys_v); grid on;
xlabel('Time [s]'); ylabel('differentiated angle [deg]');
legend('data (validation)', 'model');
title('Validation dataset (not used for identification)');

%% 5. VALIDACAO (b): resposta temporal simulada (plausivel?)
% Degrau no motor: como o motor e um integrador, o angulo deve subir em rampa.
figure('Name', 'Degrau');
tsim = 0:Ts:5;
step(sys, tsim); grid on;
xlabel('Time'); ylabel('Angle [deg]');
title('Step response of the identified model (u = 1 V)');

%% 6. VALIDACAO (c): resposta em frequencia (plausivel?)
% Em baixas frequencias: -20 dB/dec e fase -90 graus (integrador).
figure('Name', 'Bode');
bode(sys); grid on;

figure('Name', 'Polos e zeros');
pzmap(sys); grid on;
fprintf('\nZeros do modelo:\n'); disp(tzero(sys));