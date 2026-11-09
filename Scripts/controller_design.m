%% controller_design.m ----------------------------------------------------
% This script aims to obtain the the matrices A, B, C and D of the state
% space model
% Another .m script file uses these matrices to design
% the gains of a state feedback controller.

%% controller_design.m  -  L8: LQ design on the identified model
%
% Does NOT need the plant. Needs modelo_estado.mat (from L7), with the
% state-space model A, B, C, D and the sampling period.
%
% What it does (what L8 asks):
%   1) loads the model and checks it (open-loop poles, controllability)
%   2) computes the LQ gain K with dlqr, using Q = C'*C and R = 100
%   3) shows the closed-loop eigenvalues of A-BK
%   4) root-square locus: closed-loop poles as R varies
%   5) simulates the state feedback (step response) for several R, to decide
%      what time-response features are "satisfactory"
%
% Needs the Control System Toolbox (dlqr, ss, step, stepinfo).
% NOTE: no local functions (R2015a in the lab does not accept them in scripts).

clear; close all; clc;

%% Constants -------------------------------------------------------------
pasta = fileparts(mfilename('fullpath'));
addpath(fullfile(pasta, '..', 'Datasets'));     % same folder logic as identification.m

model_file = 'modelo_estado.mat';   % saved by L7
R = 100;                            % LQ weight on the control effort (start here)
R_locus = logspace(-4, 6, 80);      % range of R for the root-square locus
R_cmp   = [1 10 100 1000];          % values of R to compare in simulation
r_amp   = 10;                       % [deg] size of the reference step in simulation
t_end   = 2;                        % [s] simulation length

%% 1. Load the model and check it ---------------------------------------
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

%% 2. LQ gain -------------------------------------------------------------
% Cost: J = sum( x'Qx + u'Ru ).  Q = C'C penalizes the output y = Cx.
Q = C' * C;
[K, S, e_cl] = dlqr(A, B, Q, R);    % u(k) = -K x(k)

fprintf('\nLQ GAIN for R = %g:\n', R);
disp('K ='); disp(K);

%% 3. Closed-loop eigenvalues ---------------------------------------------
fprintf('Closed-loop eigenvalues (eig of A-BK):\n');
disp(e_cl);
fprintf('Moduli (all must be < 1 for stability):\n');
disp(abs(e_cl));

%% 4. Root-square locus: closed-loop poles as R varies --------------------
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

%% 5. Simulation of the state feedback for several R ----------------------
% To follow a reference r we use the pre-gain Nbar (so that the gain r -> y is 1):
%   u(k) = -K x(k) + Nbar r,   Nbar = Nu + K Nx
fprintf('\n   R      rise [s]  settling [s]  overshoot [%%]  max|u| [V]\n');
fprintf('------------------------------------------------------------\n');

figure('Name', 'Step response of the LQ state feedback');
for i = 1:numel(R_cmp)
    Ki = dlqr(A, B, Q, R_cmp(i));
    N  = [A - eye(n), B; C, 0] \ [zeros(n, 1); 1];
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

fprintf('\nDecide what you call "satisfactory": overshoot, settling time and\n');
fprintf('control effort (max|u|). Smaller R = faster, but larger u and it excites\n');
fprintf('unmodelled dynamics.\n');