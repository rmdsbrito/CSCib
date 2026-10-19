%% L5.m  -  Processar dados e identificar o modelo ARMAX
%
% Nao precisa da planta. Parte de um .mat recolhido na L4.
%
% Como usar:
%   1) Poe o nome do ficheiro em 'ficheiro'.
%   2) Corre (F5). Aparece a TABELA DE FITS por ordem.
%   3) Escolhe a MENOR ordem a partir da qual o fit quase nao melhora.
%   4) Poe essas ordens em na, nb, nc e corre outra vez.
%   5) Copia a tabela e o grafico "modelo vs dados" para o relatorio.
%
% NOTA: sem funcoes locais (o MATLAB R2015a do lab nao as aceita em scripts).

clear; close all; clc;

%% CONFIGURACAO
ficheiro = 'L04_1v.mat';   % dataset a identificar

Kp = 36;          % [graus/V] potenciometro
Ke = 4.46;        % [graus/V] extensometro
Ts = 0.02;        % [s] periodo de amostragem

col_u = 1;        % coluna do comando do motor
col_theta = 2;    % coluna do potenciometro
col_alpha = 3;    % coluna do extensometro

t_descarta = 10;  % [s] transitorio inicial a deitar fora
af = 0.8;         % polo do filtro passa-baixo

% Ordens do modelo final (muda DEPOIS de veres a tabela)
na = 1;  nb = 2;  nc = 1;
nk = 1;           % atraso de 1 amostra

%% 1. CARREGAR OS DADOS
S = load(ficheiro);
X = S.input;
t = X.time(:);
sigs = X.signals.values;

u_bruto = sigs(:, col_u);        % comando do motor [V]
thetae  = sigs(:, col_theta);    % potenciometro [V]
alphae  = sigs(:, col_alpha);    % extensometro [V]

%% 2. ANGULO TOTAL DA PONTA DA BARRA
% Se o potenciometro deu voltas, a tensao salta ~10 V. Isto corrige os saltos.
saltos = [0; diff(thetae)];
thetae = thetae - 10 * cumsum(round(saltos / 10));

% y = theta + alpha  [graus]
ytrend = thetae * Kp + alphae * Ke;

%% 3. DIFERENCIAR E FILTRAR
% Diferenciar remove o integrador do motor (o angulo cresce sem parar).
% O filtro passa-baixo atenua o ruido que a derivada amplifica.
Afilt = [1 -af];
Bfilt = (1 - af) * [1 -1];
yf = filter(Bfilt, Afilt, ytrend);

%% 4. DESCARTAR O INICIO E TIRAR A MEDIA DA ENTRADA
idx = t >= t(1) + t_descarta;
tt  = t(idx);
yf  = yf(idx);
u   = u_bruto(idx) - mean(u_bruto(idx));

figure('Name', 'Dados processados');
subplot(2,1,1); plot(tt, u); grid on;
xlabel('Time [s]'); ylabel('u (detrended) [V]');
subplot(2,1,2); plot(tt, yf); grid on;
xlabel('Time [s]'); ylabel('differentiated, filtered angle [deg]');

%% 5. TABELA DE FITS (testa varias ordens)
fit = @(y, ys) 100 * (1 - norm(y - ys) / norm(y - mean(y)));  % ajuste em %
z   = iddata(yf, u, Ts);                       % formato da toolbox
opt = armaxOptions('Display', 'off');

res = [];
for na_i = 1:4
    for nb_i = 1:4
        try
            th_i = armax(z, [na_i nb_i na_i nk], opt);   % nc = na
            ys_i = filter(th_i.B, th_i.A, u);            % simula o modelo
            res  = [res; na_i nb_i na_i fit(yf, ys_i)]; %#ok<AGROW>
        catch
            res  = [res; na_i nb_i na_i NaN]; %#ok<AGROW>
        end
    end
end

fprintf('\n  na  nb  nc | fit [%%]\n');
fprintf('--------------------\n');
fprintf(' %3d %3d %3d | %6.1f\n', res');
fprintf('\nEscolhe a MENOR ordem a partir da qual o fit quase nao melhora.\n');
save('Datasets/fits_ordens.mat', 'res');

%% 6. MODELO FINAL (com as ordens na, nb, nc definidas em cima)
th = armax(z, [na nb nc nk], opt);
disp(th);
A1 = th.A;        % denominador do modelo
B1 = th.B;        % numerador do modelo

ys = filter(B1, A1, u);
fprintf('Fit do modelo final (na=%d nb=%d nc=%d): %.1f %%\n', na, nb, nc, fit(yf, ys));

figure('Name', 'Modelo vs dados');
plot(tt, yf, tt, ys); grid on;
xlabel('Time [s]'); ylabel('differentiated angle [deg]');
legend('data', 'model');

save('modelo_armax.mat', 'A1', 'B1', 'na', 'nb', 'nc', 'nk', 'af', 'Kp', 'Ke', 'Ts');