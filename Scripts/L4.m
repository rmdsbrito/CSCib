%% L4.m  -  Verificar os dados recolhidos no bench
%
% Para que serve: depois de cada experiencia, ver se os dados ficaram bons
% ANTES de sair do bench (dados maus nao se podem recolher na L5).
%
% Como usar:
%   1) Depois de correr o Simulink, guarda os dados no MATLAB com:
%         save('sq_1V_0p4Hz.mat', 'input')
%   2) Poe o nome do(s) ficheiro(s) na lista 'ficheiros' aqui em baixo.
%   3) Corre este script (F5).
%   4) Compara os graficos com a Figura 11 do guia.

clear; close all; clc;

%% CONFIGURACAO
ficheiros = {'L04_1v_good_v.mat'};   % lista de ficheiros a verificar
% ex: {'sq_1V_0p4Hz.mat', 'sq_amp2.mat', 'valid.mat'}
Ts = 0.02;        % [s] periodo de amostragem esperado
col_u = 1;        % coluna do comando do motor
col_theta = 2;    % coluna do potenciometro
col_alpha = 3;    % coluna do extensometro

%% VERIFICACAO
for k = 1:numel(ficheiros)
    S = load(ficheiros{k});
    X = S.input;                       % variavel guardada pelo Scope
    t = X.time(:);
    sigs = X.signals.values;

    % Warnings para debug -------------------------------------------------------------------------------
    % Informacao basica
    Ts_real = mean(diff(t));
    fprintf('\n%s\n', ficheiros{k});
    fprintf('  amostras: %d | duracao: %.1f s | Ts medido: %.4f s\n', ...
        numel(t), t(end) - t(1), Ts_real);
    if abs(Ts_real - Ts) > 1e-4
        fprintf('  AVISO: Ts devia ser %.3f s. Confirma o Simulink!\n', Ts);
    end
    if size(sigs, 2) ~= 3
        fprintf('  AVISO: esperava 3 colunas (u, pot, strain) e ha %d.\n', size(sigs, 2));
    end

    % O potenciometro volta a -5 V quando chega a +5 V (volta completa)
    saltos = sum(abs(diff(sigs(:, col_theta))) > 5);
    if saltos > 0
        fprintf('  AVISO: o potenciometro deu %d volta(s) completa(s)! Repete a experiencia.\n', saltos);
    else
        fprintf('  potenciometro: sem voltas completas (OK)\n');
    end
    % --------------------------------------------------------------------------------------------------
    

    % Graficos ------------------------------------------------------------
    figure('Name', ficheiros{k});
    subplot(3,1,1); plot(t, sigs(:, col_u)); grid on;
    title(ficheiros{k}, 'Interpreter', 'none');
    xlabel('Time [s]'); ylabel('u [V]');
    subplot(3,1,2); plot(t, sigs(:, col_theta)); grid on;
    xlabel('Time [s]'); ylabel('\theta_e [V]');
    subplot(3,1,3); plot(t, sigs(:, col_alpha)); grid on;
    xlabel('Time [s]'); ylabel('\alpha_e [V]');
    %----------------------------------------------------------------------
end

fprintf('\nO que deves ver (Figura 11):\n');
fprintf('  u:     onda quadrada\n');
fprintf('  theta: rampa a subir e a descer, sem saltos\n');
fprintf('  alpha: picos nas transicoes do comando\n');