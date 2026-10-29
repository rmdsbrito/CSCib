%% L7.m  -  Revisao da identificacao e modelo final
%
% Nao precisa da planta. Precisa de:
%   - o dataset de identificacao (ex: sq_1V_0p4Hz.mat)
%   - o dataset de validacao     (ex: valid.mat) - NAO usado para identificar
%
% O que faz (o que o guia pede na L7):
%   PARTE A: varre o polo do filtro (af), o periodo de amostragem (Ts) e as
%            ordens (na, nb), e compara os modelos no dataset de validacao.
%   PARTE B: com a configuracao que escolheres, fixa o modelo A, B, C, D e
%            desenha o que o Q3 pede (polos/zeros, frequencia, tempo).
%
% Como usar:
%   1) Poe os nomes dos ficheiros e corre (F5). Le a tabela da PARTE A.
%   2) Poe a configuracao escolhida em 'af_f', 'dec_f', 'na_f', 'nb_f', 'nc_f'.
%   3) Corre de novo. A PARTE B guarda modelo_estado.mat e as figuras.
%
% NOTA: sem funcoes locais (o MATLAB R2015a do lab nao as aceita em scripts).

clear; close all; clc;

%% CONFIGURACAO
ficheiro_ident = 'L04_1v.mat';
ficheiro_valid = 'L04_0.75v.mat';

Kp = 36;  Ke = 4.46;  Ts = 0.02;            % constantes e Ts da recolha
col_u = 1;  col_theta = 2;  col_alpha = 3;
t_descarta = 10;                            % [s]
nk = 1;

% Valores a comparar na PARTE A
lista_af  = [0.6 0.7 0.8 0.9];              % polo do filtro
lista_dec = [1 2];                          % 1 -> Ts=0.02 s ; 2 -> Ts=0.04 s
lista_na  = 2:4;                            % nc = na
lista_nb  = 1:3;

% Configuracao FINAL (escolhes DEPOIS de veres a tabela)
af_f  = 0.8;
dec_f = 1;
na_f = 4;  nb_f = 3;  nc_f = 4;

fit = @(y, ys) 100 * (1 - norm(y - ys) / norm(y - mean(y)));  % ajuste em %
opt = armaxOptions('Display', 'off');

%% CARREGAR OS DOIS DATASETS
ficheiros = {ficheiro_ident, ficheiro_valid};
for k = 1:2
    S = load(ficheiros{k});
    X = S.input;
    t = X.time(:);
    sigs = X.signals.values;
    thetae = sigs(:, col_theta);
    alphae = sigs(:, col_alpha);
    saltos = [0; diff(thetae)];
    thetae = thetae - 10 * cumsum(round(saltos / 10));    % tira saltos do pot.
    dat(k).t = t;
    dat(k).u = sigs(:, col_u);
    dat(k).y = thetae * Kp + alphae * Ke;                 % angulo total [graus]
end

%% PARTE A: COMPARAR CONFIGURACOES
res = [];  n_instaveis = 0;
for af_i = lista_af
    for dec_i = lista_dec
        Ts_i  = Ts * dec_i;
        Afilt = [1 -af_i];
        Bfilt = (1 - af_i) * [1 -1];
        for k = 1:2
            tk = dat(k).t(1:dec_i:end);                   % sub-amostragem
            yk = filter(Bfilt, Afilt, dat(k).y(1:dec_i:end));
            uk = dat(k).u(1:dec_i:end);
            idx = tk >= tk(1) + t_descarta;
            Y{k} = yk(idx);
            U{k} = uk(idx) - mean(uk(idx));
        end
        z_i = iddata(Y{1}, U{1}, Ts_i);
        for na_i = lista_na
            for nb_i = lista_nb
                try
                    th_i = armax(z_i, [na_i nb_i na_i nk], opt);
                    if any(abs(roots(th_i.A)) >= 1)       % modelo instavel
                        n_instaveis = n_instaveis + 1;
                        continue;
                    end
                    f_id  = fit(Y{1}, filter(th_i.B, th_i.A, U{1}));
                    f_val = fit(Y{2}, filter(th_i.B, th_i.A, U{2}));
                    res = [res; af_i Ts_i na_i nb_i na_i f_id f_val]; %#ok<AGROW>
                catch
                end
            end
        end
    end
end

[~, ix] = sort(res(:, 7), 'descend');                    % melhor validacao primeiro
res = res(ix, :);
n_show = min(20, size(res, 1));
fprintf('\nTOP %d configuracoes (ordenadas pelo fit de VALIDACAO)\n', n_show);
fprintf('  af    Ts[s]  na  nb  nc | fit ident [%%] | fit valid [%%]\n');
fprintf('------------------------------------------------------------\n');
fprintf(' %.2f   %.3f  %2d  %2d  %2d | %13.1f | %12.1f\n', res(1:n_show, :)');
fprintf('(%d modelos instaveis descartados)\n', n_instaveis);
save('fits_L7.mat', 'res');

%% PARTE B: MODELO FINAL (configuracao escolhida)
Ts_f = Ts * dec_f;
Afilt = [1 -af_f];
Bfilt = (1 - af_f) * [1 -1];
for k = 1:2
    tk = dat(k).t(1:dec_f:end);
    yk = filter(Bfilt, Afilt, dat(k).y(1:dec_f:end));
    uk = dat(k).u(1:dec_f:end);
    idx = tk >= tk(1) + t_descarta;
    tt{k} = tk(idx);
    Y{k}  = yk(idx);
    U{k}  = uk(idx) - mean(uk(idx));
end

th = armax(iddata(Y{1}, U{1}, Ts_f), [na_f nb_f nc_f nk], opt);
disp(th);
A1 = th.A;  B1 = th.B;                                    % polinomios em q^-1
fprintf('Fit identificacao: %.1f %% | Fit validacao: %.1f %%\n', ...
        fit(Y{1}, filter(B1, A1, U{1})), fit(Y{2}, filter(B1, A1, U{2})));

% Voltar a pôr o integrador (1 - q^-1) e converter para estado
den_int = conv(A1, [1 -1]);
num_int = [B1, zeros(1, length(den_int) - length(B1))];
[A, B, C, D] = tf2ss(num_int, den_int);
sys = ss(A, B, C, D, Ts_f);

fprintf('\nMODELO DE ESTADO FINAL (n = %d, Ts = %.3f s):\n', size(A, 1), Ts_f);
disp('A ='); disp(A);  disp('B ='); disp(B);
disp('C ='); disp(C);  disp('D ='); disp(D);
fprintf('Polos em malha aberta:\n');  disp(eig(A));
fprintf('Zeros:\n');                  disp(tzero(sys));
save('modelo_estado.mat', 'A', 'B', 'C', 'D', 'Ts_f', 'af_f', 'dec_f', ...
     'na_f', 'nb_f', 'nc_f', 'nk', 'Kp', 'Ke');

% Graficos para o Q3
f1 = figure('Name', 'Validacao');
plot(tt{2}, Y{2}, tt{2}, filter(B1, A1, U{2})); grid on;
xlabel('Time [s]'); ylabel('differentiated angle [deg]');
legend('data (validation)', 'model');
title('Validation dataset (not used for identification)');

f2 = figure('Name', 'Polos e zeros');   pzmap(sys);  grid on;
f3 = figure('Name', 'Bode');            bode(sys);   grid on;
f4 = figure('Name', 'Degrau');
step(sys, 0:Ts_f:5);  grid on;
xlabel('Time'); ylabel('Angle [deg]');
title('Open-loop step response of the model (u = 1 V)');

saveas(f1, 'L7_validacao.png');   saveas(f2, 'L7_polos_zeros.png');
saveas(f3, 'L7_bode.png');        saveas(f4, 'L7_degrau.png');
fprintf('\nFiguras guardadas: L7_validacao, L7_polos_zeros, L7_bode, L7_degrau (.png)\n');