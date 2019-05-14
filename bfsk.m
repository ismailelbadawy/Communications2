% Clear the memory of the executed file
clc
% Clear the console output
clear

% constants
twopi = 2 * pi;

% Values that can change, but kept constant for this project.
n0 = 4;
n1 = 1;
Ts = 1;
N = 40;

% (twopi / Tb) = (twopi * N / Ts)
Tb = N * Ts;

NumberOfBits = 100;

% in BFSK we have two frequencies so we need to calculate those two frequencies.
w1 = twopi * (n0 + n1) / Tb;
w2 = twopi * (n0 - n1) / Tb;

% omega 
wc = 4 * twopi / Tb;

% we need to generate a random stream of bits, this will be our message
bk = randi([0 1], NumberOfBits, 1);

% change the zeros into -1s
bk(bk == 0) = -1;

% change the bits into a rectangular pulse train
message =  rectpulse(bk, Tb);

% Plot the rectangular pulse train of the message.
figure;
plot(message);
ylim([-2 2]);
title('The raw rectangular pulse message');
xlabel('Time')
ylabel('Volts')


% we need to plot the PSD of the train of pulses
% This is needed to calculate the effective bandwidth
[psd, f] = periodogram(message, [], [], 1);
figure;
plot(f, psd);
title("Message power spectral density");

