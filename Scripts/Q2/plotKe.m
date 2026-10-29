load("/home/rodrigobrito/Documents/Rodrigo/MEEC/CSCib/LABS/CSCib/Datasets/KP_LAB3_4.mat")

figure;

plot(tensao_pot.time, tensao_pot.signals.values(:,2), ...
    'DisplayName', '\alpha_e [V]');

xlim([0 90]);

xlabel('Time [s]');
ylabel('\alpha_e [V]');

legend('show', 'Location', 'best');

grid on;

figure;
plot(input.time, input.signals.values(:,2), ...
    'DisplayName', 'u(t) [V] - input');
hold on; % Essencial para não apagar o primeiro gráfico

plot(input.time, input.signals.values(:,1), ...
    'DisplayName', '\theta_e [V]');

xlim([0 20]);
xlabel('Time [s]');
ylabel('Amplitude [V]'); % Ajusta o rótulo se necessário, pois misturas V e \theta_e
legend('show', 'Location', 'best');
grid on;
hold off; % Boa prática para repor o estado original dos eixos