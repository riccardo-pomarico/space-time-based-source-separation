%% DAAP Homework #2 Group composition:
% Marazzi Alice alice.marazzi@mail.polimi.it
% Pomarico Riccardo riccardo.pomarico@mail.polimi.it

close all
clearvars
clc

% Loading audio files

[y1, Fs] = audioread("y1.wav");
y2 = audioread("y2.wav");
s1 = audioread("s1.wav");
s2 = audioread("s2.wav");
s3 = audioread("s3.wav");

% Reducing peaks of the signals

y1 = 0.6*y1;
y2 = 0.6*y2;

% Defining the intervals

len = min(length(y1), length(y2));
y1 = y1(1:len);
y2 = y2(1:len);

% Since we have more sources than microphones we need to take advantage of 
% spatial information in order to get a good output

% Converting angles from degrees to radians

theta_1 = deg2rad(30);
theta_2 = deg2rad(85);
theta_3 = deg2rad(-40);

% Computing TDOA

TDOA_1 = (0.75/340) + 1000*(0.09 * sin(theta_1) / 340);
TDOA_2 = (0.75/340) + 1000*(0.09 * sin(theta_2) / 340);
TDOA_3 = (0.75/340) + 1000*(0.09 * sin(theta_3) / 340);

% The attenuation of the signal received from the relative source of the 
% second microphone with respect to the first one

attenuation_first_1 = sqrt(1/(1 + (1-cos(theta_1))^2));
attenuation_second_1 = sqrt(1/(1 + (1-cos(theta_2))^2));
attenuation_third_1 = sqrt(1/(1 + (1-cos(theta_3))^2));

% The attenuation of the signal received from the relative source of the 
% first microphone with respect to the second one

attenuation_first_2 = 1/attenuation_first_1;
attenuation_second_2 = 1/attenuation_second_1;
attenuation_third_2 = 1/attenuation_third_1;

% We use matrix parameter as a cluster centroid in the kmeans algorithm

matrix = [attenuation_first_1, attenuation_first_2, TDOA_1; 
    attenuation_second_1, attenuation_second_2, TDOA_2; 
    attenuation_third_1, attenuation_third_2, TDOA_3];

% Defining the block size

window_len = 1024;
hop = window_len/16;
blocks = floor((len-window_len)/hop)+1;

% Hanning Window

window = hann(window_len);

% Initializing output

output_s1 = zeros(len, 1);
output_s2 = zeros(len, 1);
output_s3 = zeros(len, 1);

A_1 = zeros(blocks-1, window_len/2);
A_2 = zeros(blocks-1, window_len/2);
P = zeros(blocks-1, window_len/2);

% Parameters in order to save all the values of the binary masks

binarymask1 = zeros(blocks-1, window_len);
binarymask2 = zeros(blocks-1, window_len);
binarymask3 = zeros(blocks-1, window_len);

for i = 1 : blocks-1

    % Windowing
    
    y1_block = window.*y1((i-1) * hop + 1 : (i-1) * hop + window_len);
    y2_block = window.*y2((i-1) * hop + 1 : (i-1) * hop + window_len);

    % FFT
    
    y1_fft = fft(y1_block);
    y2_fft = fft(y2_block);

    for omega = 1 : window_len/2
       
        % Computing 3-dimensional feature vector

        P(i, omega) = 1/(2*pi)*angle(y2_fft(omega)/y1_fft(omega));
        A_1(i, omega) = abs(y1_fft(omega)) / sqrt (abs(y1_fft(omega))^2 + abs(y2_fft(omega))^2);
        A_2(i, omega) = abs(y2_fft(omega)) / sqrt (abs(y1_fft(omega))^2 + abs(y2_fft(omega))^2);
        phi(omega, :) = [A_1(i, omega); A_2(i, omega); P(i, omega)];

    end
    
    cluster = zeros(size(y1_fft));
    cluster(1:window_len/2) = kmeans(phi, 3, 'Start', matrix);
    
    for omega = 1 : window_len/2
        cluster(end-omega+1) = cluster(omega+1);
    end

    binarymask1_block = zeros(size(y1_fft));
    binarymask2_block = zeros(size(y1_fft));
    binarymask3_block = zeros(size(y1_fft));

    binarymask1_block(cluster == 1) = 1;
    binarymask2_block(cluster == 2) = 1;
    binarymask3_block(cluster == 3) = 1;

    binarymask1(i, :) = binarymask1_block;
    binarymask2(i, :) = binarymask2_block;
    binarymask3(i, :) = binarymask3_block;

    S_1 = y1_fft.*binarymask1_block;
    S_2 = y1_fft.*binarymask2_block;
    S_3 = y1_fft.*binarymask3_block;

    s_1_block = ifft(S_1);
    s_2_block = ifft(S_2);
    s_3_block = ifft(S_3);

    % Computing overlap and add

    output_s1((i-1)*hop+1 : (i-1)*hop+window_len) = output_s1((i-1)*hop+1 : (i-1)*hop+window_len)+s_1_block;
    output_s2((i-1)*hop+1 : (i-1)*hop+window_len) = output_s2((i-1)*hop+1 : (i-1)*hop+window_len)+s_2_block;
    output_s3((i-1)*hop+1 : (i-1)*hop+window_len) = output_s3((i-1)*hop+1 : (i-1)*hop+window_len)+s_3_block;

end

% Plotting the log-amplitude spectograms for the mixture signals

[S, F, T] = spectrogram(y1, window, 0, window_len);
[S_2, F_2, T_2] = spectrogram(y2, window, 0, window_len);

figure(1)
set(gcf,'WindowState','maximized')
subplot(1, 2, 1)
surf(T, F, abs(S))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Logarithmic Amplitude Spectrogram';'Y_1'})

subplot(1, 2, 2)
surf(T_2, F_2, abs(S_2))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Logarithmic Amplitude Spectrogram';'Y_2'})

% Plotting the log-amplitude spectrograms of the true source signals

[S, F, T] = spectrogram(s1, window, 0, window_len);
[S_2, F_2, T_2] = spectrogram(s2, window, 0, window_len);
[S_3, F_3, T_3] = spectrogram(s3, window, 0, window_len);

figure(2)
set(gcf,'WindowState','maximized')
subplot(1, 3, 1)
surf(T, F, abs(S))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Real Signals Logarithmic Amplitude Spectrograms';'S_1'})

subplot(1, 3, 2)
surf(T_2, F_2, abs(S_2))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Real Signals Logarithmic Amplitude Spectrograms';'S_2'})

subplot(1, 3, 3)
surf(T_3, F_3, abs(S_3))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Real Signals Logarithmic Amplitude Spectrograms';'S_3'})

% Plotting the log-amplitude spectograms of the estimated source signals

[S, F, T] = spectrogram(output_s1, window, 0, window_len);
[S_2, F_2, T_2] = spectrogram(output_s2, window, 0, window_len);
[S_3, F_3, T_3] = spectrogram(output_s3, window, 0, window_len);

figure(3)
set(gcf,'WindowState','maximized')
subplot(1, 3, 1)
surf(T, F, abs(S))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Estimated Logarithmic Amplitude Spectrogram';'S_1'})

subplot(1, 3, 2)
surf(T_2, F_2, abs(S_2))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Estimated Logarithmic Amplitude Spectrogram';'S_2'})

subplot(1, 3, 3)
surf(T_3, F_3, abs(S_3))
view([0 90])
axis tight
xlabel('Time [s]')
ylabel('Frequency [Hz]')
set(gca, 'YScale', 'log')
title ({'Estimated Logarithmic Amplitude Spectrogram';'S_3'})

% Plotting binary masks

figure(4)
set(gcf,'WindowState','maximized')
subplot(1, 3, 1)                      
imagesc(binarymask1', [0 1]);          
colormap(gray);                              
axis tight                                  
xlabel('Time Frame')
ylabel('Frequency Bin')
title ({'Binary Mask';'M_1'})

subplot(1, 3, 2)                        
imagesc(binarymask2', [0 1]);          
colormap(gray);                              
axis tight                                   
xlabel('Time Frame')
ylabel('Frequency Bin')
title ({'Binary Mask';'M_2'})

subplot(1, 3, 3)                        
imagesc(binarymask3', [0 1]);          
colormap(gray);                              
axis tight                                   
xlabel('Time Frame')
ylabel('Frequency Bin')
title ({'Binary Mask';'M_3'})

% Plotting density plot each pair of features with the histograms of the
% individual features

figure(5)
set(gcf,'WindowState','maximized')
subplot(3, 3, 1)
histogram2(A_1, A_2, 'DisplayStyle', 'tile');
xlabel('A_1')
ylabel('A_2')
zlabel('Elements [int]')
title ({'Density Plot'; '[A_1(m,w), A_2(m,w)]'})

subplot(3, 3, 2)
histogram2(A_1, P, 'DisplayStyle', 'tile');
xlabel('A_1')
ylabel('P')
zlabel('Elements [int]')
title ({'Density Plot'; '[A_1(m,w), P(m,w)]'})

subplot(3, 3, 3)
histogram2(A_2, P, 'DisplayStyle','tile');
xlabel('A_2')
ylabel('P')
zlabel('Elements [int]')
title ({'Density Plot'; '[A_2(m,w), P(m,w)]'})

subplot(3, 3, 4)
histogram(A_1)
xlabel('A_1')
ylabel('Elements [int]')
title ({'Histogram'; 'A_1(m,w)'})

subplot(3, 3, 5)
histogram(A_2)
xlabel('A_2')
ylabel('Elements [int]')
title ({'Histogram'; 'A_2(m,w)'})

subplot(3, 3, 6)
histogram(P)
xlabel('P')
ylabel('Elements [int]')
title ({'Histogram'; 'P(m,w)'})

audiowrite('Marazzi_Pomarico_s1.wav', output_s1, Fs)
audiowrite('Marazzi_Pomarico_s2.wav', output_s2, Fs)
audiowrite('Marazzi_Pomarico_s3.wav', output_s3, Fs)