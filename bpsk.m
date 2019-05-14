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
message = rectpulse(bk, Tb);

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

% We have to modulate the message using BFSK and BPSK for every different amplitude.
t = 1 : NumberOfBits * Tb;
FSK = t * 0;
PSK = t * 0;
for k = 1 : 9
    for bit = 1 : NumberOfBits * Tb
        if(message(bit) == 1)
            % we need to set both the FSK and PSK
            PSK(bit) = A(k) * cos(wc * bit);
        else
            PSK(bit) = -1 * A(k) * cos(wc * bit);
        end
    end

    % We need to try this 20 times and get different probability of error
    for trial = 1 : 20
        % Simulate that we pass the signal through the channel meaning 
        % we add the signal and the noise
        % Should be the same noise but we need to
        % make sure it corresponds to the length of each modulation technique
        pskNoise = wgn(1, length(PSK), No / 2);
        
        receivedPSK = PSK + pskNoise;

        % we need to create a matched filter for each amplitude
        % Note that a matched filter h(t) = s(Tb - t)
        % We also need to set the value of the output of the filter,
        % the output to the filter is the convolution between received signal and the matched filter
        for time = 1 : Tb
            matchedFilterPSK(time) = A(k) * cos(wc * (Tb - time));
        end
        outputFilterPSK = conv(receivedPSK, matchedFilterPSK);

        if (trial == 1)
            % Now we need to plot the transmitted signals and the received
            % Plots starting with PSK or FSK are the transmitted and received signal plots.
            % We also need to plot the output of the filters
            snrValue = (-5 + k);
            if(snrValue == 3 || snrValue == -1)
                
                % PSK Plotting
                figure;
                subplot(3, 1, 1);
                plot(1:NumberOfBits * Tb, PSK(1: NumberOfBits * Tb));
                title("Transmitted PSK for SNR = " + snrValue + " and amplitude = " + A(k));
                ylim([-1 1]);
                
                subplot(3, 1, 2);
                plot(1:NumberOfBits * Tb, receivedPSK(1: NumberOfBits * Tb));
                title("Received signal ");
                ylim([-5 5]);
                
                % plot the filter output
                subplot(3, 1, 3);
                plot(1 : NumberOfBits * Tb, outputFilterPSK( 1: NumberOfBits * Tb));
                title("Filter output ");
        
            end
        end

        index = 1;
        for time = Tb : Tb : length(outputFilterPSK)
            if(outputFilterPSK(time) > 0)
                decidedPSK(index) = 1;
            else
                decidedPSK(index) = -1;
            end
            index = index + 1;
        end

        erroneousBitsPSK = 0;
        bk = transpose(bk);
        for index = 1 : length(decidedPSK)
            if(decidedPSK(index) ~= bk(index))
                erroneousBitsPSK = erroneousBitsPSK + 1; 
            end
        end
        Pe(trial) = erroneousBitsPSK / length(decidedPSK);
    end
    Pavg(k) = mean(transpose(Pe));
end

Pexact = (1 / 2) * erfc(sqrt((A .^ 2) * Tb /4));
SNR = -4 : 4;
figure;
semilogy(SNR, Pavg, 'b');
hold on;
semilogy(SNR, Pexact, 'r');
title("BPSK Simulated Pe and Exact Pe");
xlabel("SNR in dB"); ylabel('Bit error rate');
legend("Bit error rate simulation", "Bit error rate calculated");

