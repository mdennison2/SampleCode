%% photoPlot
% This function loads in photocell data and performs and FFT
% to ensure that a stimulus is flickering at the intended frequency

% load the data
file = 'flicker15.cnt';
d = readantcnt(file);

% get channel labels
channelLabel = d.channellabel';
photoLabel = channelLabel(71:72);

% store the data and sampling rate and
% take exactly the desired number of samples
Fs = d.samplingrate;
testDuration = Fs*15;

% parse out the data
photoData = -d.data(1:testDuration,71:72);

% do FFT
windowSize = testDuration;
freqs = linspace(0,Fs/2,windowSize/2+1);
maxFreqBin = find(freqs == 80);
freqs = freqs(2:maxFreqBin);

fcoefs = fft(photoData,[],1)/length(photoData);
fcoefs = fcoefs(2:maxFreqBin,:);

% get power and plot
power = 2*abs(fcoefs).^2;
maxY = max(max(power));

figure;

subplot(1,2,1);
plot(freqs,power(:,1),'LineWidth',1.5);
ylim([0 maxY]);
set(gca,'FontSize',14);
title('Photocell 1','FontSize',18);
xlabel('frequency (Hz)','FontSize',18);
ylabel('power (mV^2/Hz)','FontSize',18);

subplot(1,2,2);
time = linspace(0,500,1024);
plot(time,photoData(1:1024,1));
set(gca,'FontSize',14);
xlabel('time (ms)','FontSize',18);
ylabel('amplitude (mV)','FontSize',18);
