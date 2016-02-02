function [time, simEEG] = simulateEEG(N,fs)
% Create simulated EEG data with a 1/f frequency distribution.
%
% INPUTS:
% N = length of simulated EEG signal
% fs = sampling frequency
%
% OUTPUT:
% simEEG = simulated EEG data
time = (1:N)/fs; % Define time vector
x = randn(N,1); % signal starts as Gaussian noise
xfft = fft(x); % take FFT
fVec = linspace(0, fs/2, N/2+1); % Define frequency vector
fVec = fVec(2:end); % Remove zero because we divide by fVec
xfft2 = xfft(1:(N/2))./fVec'; % Multiply power spectrum by 1/f
xf_filtered = [xfft2; flipud(xfft2)]; % Add coefficients for negative frequencies
simEEG = real(ifft(xf_filtered')); % simEEG is real part of IFFT