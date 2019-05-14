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
% Add the noise No
No = 2;

% (twopi / Tb) = (twopi * N / Ts)
Tb = N * Ts;

NumberOfBits = 50;

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
% The plot is provided in a screenshot called "Original Message.JPG"
figure;
plot(message);
ylim([-2 2]);
title('The raw rectangular pulse message');
xlabel('Time')
ylabel('Volts')


% we need to plot the PSD of the train of pulses
% The theoretical bandwidth of this pulse Rb is equal to 1/Tb = 1/40 = 0.025
% This is needed to calculate the effective bandwidth
[psd, f] = periodogram(message, [], [], 1);
figure;
plot(f, psd);
title("Message power spectral density");
% On the plot the effective bandwidth would be the position of the first null
% That is provided through a screenshot called "PSD (Effective bandwidth).JPG"

% Note that : SNR changes from -4 to 4 dB.
% We need to generate the A for each SNR value
SNR = -4 : 1 : 4;

% We know from the description of the system that SNR = 10 log ((A^2) * Tb/ 2No)
% and we also know that Tb & No are givens so A would be the only dependent variable if we know the SNR
A = sqrt(10.^(SNR / 10) * 2 * No/ Tb);

% We have to modulate the message using BFSK for every different amplitude.
t = 1 : NumberOfBits * Tb;
FSK = t * 0;
for amplitude = A
    for bit = 1 : NumberOfBits * Tb
        if(message(bit) == 1)
            FSK(A == amplitude, bit) = amplitude * cos(w1 * bit);
        else
            FSK(A == amplitude, bit) = amplitude * cos(w2 * bit);
        end
    end
end


% Now we need to plot every transmitted signal.
figure;
plot(1:NumberOfBits * Tb, FSK(9, 1: NumberOfBits * Tb));
title('FSK for SNR = 4 and amplitude = 0.50');
xlim([1 400]);

