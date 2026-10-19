%% L6.m  -  Integrador, espaco de estados e validacao
%
% Nao precisa da planta. Precisa de:
%   - modelo_armax.mat  (guardado pelo L5, com as ordens que escolheste)
%   - valid.mat         (dataset recolhido na L4, NUNCA usado na identificacao)
%
% Como usar:
%   1) Confirma que modelo_armax.mat e o do modelo final que queres (o L5
%      sobrescreve este ficheiro de cada vez que corre).
%   2) Corre (F5).
%   3) Le o fit de validacao e os graficos. Guarda-os para o relatorio.
%
% Saida: A, B, C, D (modelo de estado discreto) em modelo_estado.mat
% NOTA: sem funcoes locais (o MATLAB R2015a do lab nao as aceita em scripts).

clear; close all; clc;

%% CONFIGURACAO
ficheiro_modelo = 'modelo_armax.mat';   % saida do L5
ficheiro_valid  = 'L04_0.75v.mat';          % dados so para validacao
t_descarta = 10;                        % [s] igual ao usado no L5
col_u = 1;  col_theta = 2;  col_alpha = 3;

%% 1. CARREGAR O MODELO ARMAX DO L5
M = load(ficheiro_modelo);
A1 = M.A1;   B1 = M.B1;                 % polinomios em q^-1
af = M.af;   Kp = M.Kp;   Ke = M.Ke;   Ts = M.Ts;
fprintf('Modelo: na=%d nb=%d nc=%d nk=%d\n', M.na, M.nb, M.nc, M.nk);
fprintf('Polos do ARMAX (devem estar dentro do circulo unitario):\n');
disp(abs(roots(A1)));

%% 2. ADICIONAR O INTEGRADOR E CONVERTER PARA ESTADO
% No L5 derivamos a saida (tiramos um polo em z=1). Agora voltamos a poe-lo:
% multiplica-se o denominador por (1 - q^-1).
den_int = conv(A1, [1 -1]);
% O ARMAX usa o operador de atraso (q^-1); o tf2ss usa o operador avanco (q).
% Com numerador e denominador do MESMO comprimento, os coeficientes servem
% nos dois (e o que o eqtflength do guia faz).
num_int = [B1, zeros(1, length(den_int) - length(B1))];

[A, B, C, D] = tf2ss(num_int, den_int);   % x(k+1)=Ax+Bu ; y=Cx+Du
sys = ss(A, B, C, D, Ts);

fprintf('\nModelo de estado (n = %d):\n', size(A, 1));
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('C ='); disp(C);
disp('D ='); disp(D);
fprintf('Polos em malha aberta (um deles deve ser 1, o integrador):\n');
disp(eig(A));
save('modelo_estado.mat', 'A', 'B', 'C', 'D', 'Ts');

%% 3. PROCESSAR O DATASET DE VALIDACAO (igual ao L5)
S = load(ficheiro_valid);
X = S.input;
t = X.time(:);
sigs = X.signals.values;
u_bruto = sigs(:, col_u);
thetae  = sigs(:, col_theta);
alphae  = sigs(:, col_alpha);

saltos = [0; diff(thetae)];
thetae = thetae - 10 * cumsum(round(saltos / 10));   % tira saltos do pot.
ytrend = thetae * Kp + alphae * Ke;                  % angulo total [graus]

Afilt = [1 -af];
Bfilt = (1 - af) * [1 -1];
yf = filter(Bfilt, Afilt, ytrend);                   % diferenciar + filtrar

idx = t >= t(1) + t_descarta;
tt  = t(idx);
yf_v = yf(idx);
u_v  = u_bruto(idx) - mean(u_bruto(idx));

%% 4. VALIDACAO (a): modelo ARMAX vs dados NAO usados na identificacao
fit = @(y, ys) 100 * (1 - norm(y - ys) / norm(y - mean(y)));
ys_v = filter(B1, A1, u_v);
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