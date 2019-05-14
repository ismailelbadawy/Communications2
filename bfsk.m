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
xlim([1 400]);
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

% We have to modulate the message using BFSK and BPSK for every different amplitude.
t = 1 : NumberOfBits * Tb;
FSK = t * 0;
PSK = t * 0;
for amplitude = A
    for bit = 1 : NumberOfBits * Tb
        if(message(bit) == 1)
            % we need to set both the FSK and PSK
            FSK(A == amplitude, bit) = amplitude * cos(w1 * bit);
            PSK(A == amplitude, bit) = amplitude * cos(wc * bit);
        else
            FSK(A == amplitude, bit) = amplitude * cos(w2 * bit);
            PSK(A == amplitude, bit) = -1 * amplitude * cos(wc * bit);
        end
    end
end


% Simulate that we pass the signal through the channel meaning 
% we add the signal and the noise
for amplitude = A
    % Should be the same noise but we need to
    % make sure it corresponds to the length of each modulation technique
    fskNoise = wgn(1, length(FSK), No/2);
    pskNoise = wgn(1, length(PSK), No/2);

    receivedFSK(A == amplitude, :) = FSK(A == amplitude) + fskNoise;
    receivedPSK(A == amplitude, :) = PSK(A == amplitude) + pskNoise;
end


% we need to create a matched filter for each amplitude
% Note that a matched filter h(t) = s(Tb - t)
% We also need to set the value of the output of the filter,
% the output to the filter is the convolution between received signal and the matched filter
for amplitude = A
    for time = 1 : Tb
        matchedFilter1(A == amplitude, time) = amplitude * cos(w1 * (Tb - time));
        matchedFilter2(A == amplitude, time) = amplitude * cos(w2 * (Tb - time));
    end
    outputFilter1(A == amplitude, :) = conv(receivedFSK(A == amplitude, :), matchedFilter1(A == amplitude,  :));
    outputFilter2(A == amplitude, :) = conv(receivedFSK(A == amplitude, :), matchedFilter2(A == amplitude, :));
end

% Now we need to plot the transmitted signals and the received
% Plots starting with PSK or FSK are the transmitted and received signal plots.
% We also need to plot the output of the filters
for i = 1 : 9
    snrValue = (-5 + i);
    if(snrValue == 3 || snrValue == -1)
        figure;
        subplot(4, 1, 1);
        plot(1:NumberOfBits * Tb, FSK(i, 1: NumberOfBits * Tb));
        title("Transmitted FSK for SNR = " + snrValue + " and amplitude = " + A(i));
        xlim([1 1000]);
        ylim([-1 1]);
        subplot(4, 1, 2);
        plot(1: NumberOfBits * Tb, receivedFSK(i, 1:NumberOfBits * Tb));
        title("Recevied signal ")
        xlim([1 1000]);
        ylim([-10 10]);
        subplot(4, 1, 3);
        plot(1 : NumberOfBits * Tb, outputFilter1(i, 1: NumberOfBits * Tb));
        title("Filter 1 output");
        xlim([1 1000]);
        subplot(4, 1, 4);
        plot(1 : NumberOfBits * Tb, outputFilter2(i, 1: NumberOfBits * Tb));
        title("Filter 2 output");
        xlim([1 1000]);
        % figure;
        % subplot(2, 1, 1);
        % plot(1:NumberOfBits * Tb, PSK(i, 1: NumberOfBits * Tb));
        % title("Transmitted PSK for SNR = " + snrValue + " and amplitude = " + A(i));
        % xlim([1 400]);
        % ylim([-1 1]);
        % subplot(2,1, 2);
        % plot(1:NumberOfBits * Tb, receivedPSK(i, 1: NumberOfBits * Tb));
        % title("Received signal ");
        % ylim([-5 5]);
        % xlim([1 400]);
    end
end

for amplitude = A
    
end