%% TEST THE PCA CODE

% make fake data to test PCA
Fs = 2048; % Sampling frequency
N = Fs*10; % Length of signal
time = (1:N)/Fs; % Time vector
fVec = linspace(0, Fs/2, N/2+1); % Frequency Vector
maxFreq = fVec(end);
sinChannel = sin(2*pi*8*time);

load('ANTWAVE64');
hm = ANTWAVE64;
nChan = 64;

EEG = zeros(64,N);
% Set channel Fz and Pz to contain an 8hz sine wave
EEG(6,:) = sin(2*pi*8*time);
EEG(26,:) = sin(2*pi*8*time);

% PCA
% manual PCA of the data
nTimePoints = length(EEG);
pcaEEG = EEG;
pcaEEG = pcaEEG - repmat(mean(pcaEEG,2),1,nTimePoints);
cov = (pcaEEG*pcaEEG')./(nTimePoints-1);
[eigvecs,eigvals] = eig(cov);
varExplained = 100*diag(eigvals)./sum(diag(eigvals));
varExpected = 100/size(pcaEEG,1);

% Select only the significant component wights
sigPCs = find(varExplained > varExpected);

% create PC "channels" from the significant component weights
pcChan = zeros(length(sigPCs),nTimePoints);
for k = 1:length(sigPCs)
    pcChan(k,:) = sum(repmat(eigvecs(:,sigPCs(k)),1,length(pcaEEG)).*pcaEEG);
end

% Plot the components and their variance accounted for
plot(1:size(pcaEEG,1),varExplained,'k-o')
hold on;
line([1 size(pcaEEG,1)],[100/size(pcaEEG,1) 100/size(pcaEEG,1)])
xlabel('Component Number','FontSize',16);
legend({'% Variance Explained ','Chance-Level Variance'});
set(gca,'FontSize',16);
set(gcf,'color','w');

% Topo for each significant component
topoPCs = eigvecs;
corttopo(topoPCs(:,sigPCs(1)),hm,'badchans',[13 19]);

%% TEST THE SFFT CODE

% Create a changing 10 Hz and 15 Hz sin wave
Fs = 2048;
N = Fs*120; % Length of a signal chunk
time = (1:N)/Fs; % Time vector
comboSig = [sin(2*pi*8*time(1:length(time)/2)) ...
    sin(2*pi*15*time(length(time)/2+1:end))];
pcChan = repmat(comboSig,1,4);

windowSize = Fs*1;
stepSize = windowSize/2;

maxFreq = 30;
freqs = linspace(0,Fs/2,windowSize/2+1);
freqs = freqs(2:maxFreq+1);

startInds = 1:stepSize:length(pcChan);
startInds = startInds(1:end-1);

power = zeros(length(freqs),length(startInds));
for k = 1:length(startInds)
    temp = pcChan(startInds(k):startInds(k)+windowSize-1);
    fcoefs = fft(temp,[],2)/length(temp);
    fcoefs = fcoefs(2:maxFreq+1);
    power(:,k) = abs(fcoefs).^2;
end

% Plot the fake component spectrogram
figure;
time = linspace(1,(startInds(end)+windowSize-1)/Fs/60,length(startInds));
imagesc(time,freqs,power);
set(gca,'ydir','normal','FontSize',14);
title('Test Signal','FontSize',14);
colormap hot;
ylabel('freq (Hz)','FontSize',14);
xlabel('time (min)','FontSize',14);
set(gcf,'color','w');