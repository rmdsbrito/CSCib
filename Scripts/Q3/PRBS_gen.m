clear; clc; close all;

f = 0.4;             % signal frequency [Hz]
Period = 1/f;        % Period [s]

Fs = 50;            % Samp Frequency [Hz]
NumChannel = 1;
NumPeriod = 20;

Range = [-1, 1];
Band = [0, 0.1];

NumSamples = round(Period * Fs * NumPeriod);

u1 = idinput(NumSamples, 'prbs', Band, Range);

t = (0:NumSamples-1)/Fs;

plot(t, u1)
xlim([0 20]);
xlabel('Time [s]')
ylabel('Amplitude')
title('PRBS Signal')
grid on