%% controller_design.m  -  L8: LQ design on the identified model
%
% This script requires the script modelo_estado.mat generated using
% identification.m, which contains the state space model A, B, C and D
% matrices and the sampling time to reduce sampling mistakes 
%
clear; close all; clc;

Kp = 36;  Ke = 4.46;

%% LQ design on the identified model

% Loading Parameters ------------------------------------------------------
pasta = fileparts(mfilename('fullpath'));
addpath(fullfile(pasta, '..', 'Datasets'));

model_file = 'modelo_estado.mat'; % Imported from identification.m (PATH)
R = 100;                            
R_locus = logspace(-4, 6, 80); % Necessary log curve for the root-locus    
R_cmp = [1 10 100 1000]; % Used for comparison between actuator weights         
r_amp = 10;                       
t_end = 2;                        

% Model loading and debug -------------------------------------------------
M = load(model_file);
A = M.A;  B = M.B;  C = M.C;
if isfield(M, 'Ts_f')
    Ts = M.Ts_f;
elseif isfield(M, 'Ts')
    Ts = M.Ts;
else
    error('Sampling period not found in %s (expected Ts_f or Ts).', model_file);
end
n = size(A, 1);
t = 0:Ts:t_end;

fprintf('Model: n = %d states, Ts = %.3f s\n', n, Ts);
fprintf('Open-loop poles (one should be 1, the integrator):\n');
disp(eig(A));
fprintf('rank(ctrb) = %d, rank(obsv) = %d (both should be %d)\n', ...
        rank(ctrb(A, B)), rank(obsv(A, C)), n);

% Computing the LQ gain ----- J = sum( x'Qx + u'Ru ) ----------------------
Q = C' * C;
[K, S, e_cl] = dlqr(A, B, Q, R);    % u(k) = -K x(k)

fprintf('\nLQ GAIN for R = %g:\n', R);
disp('K ='); disp(K);

% Closed-loop eigenvalues ---------------------------------------------
fprintf('Closed-loop eigenvalues (eig of A-BK):\n');
disp(e_cl);
fprintf('Moduli (all must be < 1 for stability):\n');
disp(abs(e_cl));

% Root-square locus: closed-loop poles as R varies --------------------
E = zeros(numel(R_locus), n);
for i = 1:numel(R_locus)
    Ki = dlqr(A, B, Q, R_locus(i));
    E(i, :) = eig(A - B * Ki).';
end

figure('Name', 'Root-square locus (R)');
hold on;
th = linspace(0, 2 * pi, 200);
plot(cos(th), sin(th), 'k--');                       % unit circle
for j = 1:n
    scatter(real(E(:, j)), imag(E(:, j)), 15, log10(R_locus), 'filled');
end
plot(real(eig(A)), imag(eig(A)), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);   % open loop
plot(real(e_cl), imag(e_cl), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);        % chosen R
axis equal; grid on;
xlabel('Real part'); ylabel('Imaginary part');
cb = colorbar;  ylabel(cb, 'log_{10}(R)');
legend('unit circle', 'locus', 'open-loop poles', sprintf('closed loop, R = %g', R));
title('Closed-loop poles as R varies');

% Simulation of the state feedback for several R ----------------------
% To follow a reference r we use the pre-gain Nbar (so that the gain r -> y is 1):
%   u(k) = -K x(k) + Nbar r,   Nbar = Nu + K Nx
fprintf('\n   R      rise [s]  settling [s]  overshoot [%%]  max|u| [V]\n');
fprintf('------------------------------------------------------------\n');

figure('Name', 'Step response of the LQ state feedback');
for i = 1:numel(R_cmp)
    % Defined in the LAB guide ----
    Ki = dlqr(A, B, Q, R_cmp(i));
    N = [A - eye(n), B; C, 0] \ [zeros(n, 1); 1];
    Nbar = N(n + 1) + Ki * N(1:n);

    % outputs: 1 = y (angle), 2 = u (control voltage)
    sys_cl = ss(A - B * Ki, B * Nbar, [C; -Ki], [0; Nbar], Ts);
    resp = step(r_amp * sys_cl, t);
    y = resp(:, 1);  u = resp(:, 2);

    info = stepinfo(sys_cl(1));
    fprintf(' %6g   %7.3f   %10.3f   %12.1f   %9.2f\n', R_cmp(i), ...
            info.RiseTime, info.SettlingTime, info.Overshoot, max(abs(u)));

    subplot(2, 1, 1); hold on; plot(t, y);
    subplot(2, 1, 2); hold on; plot(t, u);
end
subplot(2, 1, 1); grid on; plot(t, r_amp * ones(size(t)), 'k--');
xlabel('Time [s]'); ylabel('y [deg]');
legend([arrayfun(@(x) sprintf('R = %g', x), R_cmp, 'UniformOutput', false), {'reference'}]);
title(sprintf('Step reference of %g deg', r_amp));
subplot(2, 1, 2); grid on;
xlabel('Time [s]'); ylabel('u [V]');


%% Kalman filter, LTR and closed-loop simulation --------------------------

% Parameter Declaration  -------------------------------------------------------------
R = 100; % The R concluded as best from the last section                     
Qw = 100*eye(n);                     
Rv = 1; 
Rv_locus = logspace(-4, 6, 80); % range of Rv for the root-locus
Rv_cmp = Rv * [0.01 1 100]; % Rv values to compare in simulation
r_amp = 10; % [deg] reference step
t_end = 3; % [s] simulation length
sigma_v = 0.05;  % [deg] measurement noise used ONLY to visualize the effect of Rv

% Pre-gain so that the gain r -> y is 1 (Section 5.4 of the guide)
N = [A - eye(n), B; C, 0] \ [zeros(n, 1); 1];
Nbar = N(n + 1) + K * N(1:n);

fprintf('Model: n = %d states, Ts = %.3f s\n', n, Ts);
fprintf('LQ gain (R = %g):\n', R);  disp(K);
fprintf('Nbar = %.4f\n', Nbar);

% Kalman gain (M) ---------------------------------------------------------
% Plant with noise:  x(k+1) = A x + B u + w ,  y = C x + v
Mk = dlqe(A, eye(n), C, Qw, Rv);                % current-estimator gain (n x 1)

fprintf('\nKalman gain M (Qw = %g*I, Rv = %g):\n', Qw(1,1), Rv);  disp(Mk);
e_obs = eig((eye(n) - Mk * C) * A);             % error-equation poles of the current estimator
fprintf('Observer poles:\n');                    disp(e_obs);
fprintf('Controller poles (eig(A-BK)):\n');      disp(eig(A - B * K));

% Observer root-square locus in Rv ---------------------------------------
E = NaN(numel(Rv_locus), n);
for i = 1:numel(Rv_locus)
    try
        Mi = dlqe(A, eye(n), C, Qw, Rv_locus(i));
        E(i, :) = eig((eye(n) - Mi * C) * A).';
    catch
        % leave NaN if the Riccati equation fails for this Rv
    end
end

figure('Name', 'Observer root-square locus (Rv)');
hold on;
th = linspace(0, 2 * pi, 200);
plot(cos(th), sin(th), 'k--');                                  % unit circle
for j = 1:n
    scatter(real(E(:, j)), imag(E(:, j)), 15, log10(Rv_locus), 'filled');
end
plot(real(eig(A)), imag(eig(A)), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);   % open loop
plot(real(e_obs), imag(e_obs), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);      % chosen Rv
axis equal; grid on;
xlabel('Real part'); ylabel('Imaginary part');
cb = colorbar;  ylabel(cb, 'log_{10}(R_v)');
legend('unit circle', 'locus', 'open-loop poles', sprintf('observer, R_v = %g', Rv));
title('Observer poles as R_v varies');

% Simulation of the Full LQG loop -----------------------------------------
% Current estimator (Fig. 13) + feedback (Fig. 14), at each sample k:
%   xhat(k|k-1) = A xhat(k-1|k-1) + B u(k-1)
%   xhat(k|k)   = xhat(k|k-1) + M ( y(k) - C xhat(k|k-1) )
%   u(k)        = -K xhat(k|k) + Nbar r(k)
Nk = round(t_end / Ts) + 1;
t  = (0:Nk-1)' * Ts;
r  = r_amp * ones(Nk, 1);
rng(1);                                          % same noise for every Rv
v  = sigma_v * randn(Nk, 1);
ss_idx = t > t_end / 2;                          % second half: "steady state"

% Baseline: ideal LQ with the TRUE state (no estimator)
x = zeros(n, 1);  y_lq = zeros(Nk, 1);  u_lq = zeros(Nk, 1);
for k = 1:Nk
    y_lq(k) = C * x;
    u_lq(k) = -K * x + Nbar * r(k);
    x = A * x + B * u_lq(k);
end

fprintf('\n                 overshoot [%%]  settling [s]  max|u| [V]  std(u) steady [V]\n');
fprintf('----------------------------------------------------------------------------\n');
st = find(abs(y_lq - r_amp) > 0.02 * r_amp, 1, 'last');
if isempty(st), ts = 0; else ts = t(min(st + 1, Nk)); end
fprintf(' LQ (true state)  %10.1f  %12.3f  %10.2f  %16.3f\n', ...
    (max(y_lq) - r_amp) / r_amp * 100, ts, max(abs(u_lq)), std(u_lq(ss_idx)));

% LQG for each Rv (one column per Rv)
Y = zeros(Nk, numel(Rv_cmp));                    % true output
U = zeros(Nk, numel(Rv_cmp));                    % control voltage
for i = 1:numel(Rv_cmp)
    Mi = dlqe(A, eye(n), C, Qw, Rv_cmp(i));
    x  = zeros(n, 1);  xh = zeros(n, 1);  u_prev = 0;
    for k = 1:Nk
        Y(k, i) = C * x; % true output
        ym = Y(k, i) + v(k); % measured output (with noise)
        xp = A * xh + B * u_prev; % prediction  xhat(k|k-1)
        xh = xp + Mi * (ym - C * xp); % correction  xhat(k|k) com Mi
        U(k, i) = -K * xh + Nbar * r(k);
        x = A * x + B * U(k, i); % plant
        u_prev = U(k, i);
    end
    st = find(abs(Y(:, i) - r_amp) > 0.02 * r_amp, 1, 'last');
    if isempty(st), ts = 0; else ts = t(min(st + 1, Nk)); end
    fprintf(' Rv = %-10g %10.1f  %12.3f  %10.2f  %16.3f\n', Rv_cmp(i), ...
        (max(Y(:, i)) - r_amp) / r_amp * 100, ts, ...
        max(abs(U(:, i))), std(U(ss_idx, i)));
end

% System Response (LQR and LQG in chain) ---------------------------------------------------------------
figure('Name', 'LQG simulation');
subplot(2, 1, 1);
plot(t, y_lq, 'k', t, Y, t, r, 'k--'); grid on;
xlabel('Time [s]'); ylabel('y [deg]');
legend([{'LQ (true state)'}, ...
    arrayfun(@(x) sprintf('LQG, R_v = %g', x), Rv_cmp, 'UniformOutput', false), ...
    {'reference'}]);
title(sprintf('Step reference of %g deg (R = %g, Q_w = %g I)', r_amp, R, Qw(1, 1)));
subplot(2, 1, 2);
plot(t, u_lq, 'k', t, U); grid on;
xlabel('Time [s]'); ylabel('u [V]');