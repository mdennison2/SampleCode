%% LOADING
% load the data
baselineFile = 'Data/Mark/marktest.cnt';
d1 = readantcnt(baselineFile);
runFile = 'Data/Mark/markvr.cnt';
d2 = readantcnt(runFile);

% load the head model
load('ANTWAVE64');
hm = ANTWAVE64;
nChan = 64;

% get channel labels
channelLabel = d1.channellabel';
EEGchanLabel = channelLabel(1:nChan);

% store the data and sampling rate and
% take exactly the desired number of samples
Fs = d1.samplingrate;
baselineDuration = Fs*60*2;
expDuration = Fs*60*8;
baseline = d1.data(1:baselineDuration,1:nChan);
EEG = d2.data(1:expDuration,1:nChan);

% clear the original structure
clear d1 d2;

%% PRE-PROCESSING

% Remove mean from each channel to fix DC offset
baseline = baseline - repmat(mean(baseline)', 1, length(baseline))';
EEG = EEG - repmat(mean(EEG)', 1, length(EEG))';

% Filter the data from 1 to 60Hz and detrend
baseline = ndetrend(filtereeg(baseline,Fs));
EEG = ndetrend(filtereeg(EEG,Fs));

EEG = EEG';
baseline = baseline';

% Normalize the EEG and covert to z-score
baselineMean = mean(baseline(:,1:Fs*10),2);
zEEG = (EEG - repmat(baselineMean,1,length(EEG)))...
    ./repmat(std(baseline,[],2),1,length(EEG));

% zero out the mastoids
EEG([13 19],:) = 0;
baseline([13 19],:) = 0;

%% PCA
nTimePoints = length(zEEG);
pcaEEG = zEEG;
pcaEEG([13 19],:) = [];
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
figure;
plot(1:size(pcaEEG,1),varExplained,'k-o')
hold on;
line([1 size(pcaEEG,1)],[100/size(pcaEEG,1) 100/size(pcaEEG,1)])
xlabel('Component Number','FontSize',16);
legend({'% Variance Explained ','Chance-Level Variance'});
set(gca,'FontSize',16);
set(gcf,'color','w');

% Topo for each significant component
figure;
topoPCs = zeros(64,64);
topoPCs([1:12 14:18 20:64],[1:12 14:18 20:64]) = eigvecs;
for k = 1:length(sigPCs)
    subplot(3,2,k);
    corttopo(topoPCs(:,sigPCs(k)),hm,'badchans',[13 19],'channumbers',0);
    title(sprintf('%.2f Var Explained',varExplained(sigPCs(k))),'FontSize',14);
end
set(gcf,'color','w');

%% Short-Time FFT with overlapping windows
windowSize = Fs*1;
stepSize = windowSize/2;

maxFreq = 30;
freqs = linspace(0,Fs/2,windowSize/2+1);
freqs = freqs(2:maxFreq+1);

startInds = 1:stepSize:length(pcChan);
startInds = startInds(1:end-1);

power = zeros(length(sigPCs),length(freqs),length(startInds));
for k = 1:length(startInds)
    temp = pcChan(:,startInds(k):startInds(k)+windowSize-1);
    fcoefs = fft(temp,[],2)/length(temp);
    fcoefs = fcoefs(:,2:maxFreq+1);
    power(:,:,k) = abs(fcoefs).^2;
end

% Plot the component spectrograms
figure;
time = linspace(1,(startInds(end)+windowSize-1)/Fs/60,length(startInds));
for k = 1:length(sigPCs)
    s = subplot(3,2,k);
    imagesc(time,freqs,squeeze(power(k,:,:)));
    set(gca,'ydir','normal','FontSize',14);
    ylim([1 15]);
    title(sprintf('PC %d - %.2f%% VarExp.',k,varExplained(sigPCs(k))),'FontSize',14);
    caxis([1.96 4]);
    colormap hot;
    sPos = get(s,'position');
    if k == 2 || k == 4 || k == 6    
        c = colorbar('location','eastoutside');
        ylabel(c, 'Z-Score','FontSize',14);
    end
    if k == 1 || k == 3 || k == 5
        ylabel('freq (Hz)','FontSize',14);
    end
    if k == 5 || k == 6
        xlabel('time (min)','FontSize',14);
    end
    set(s,'position',sPos);
end
set(gcf,'color','w');