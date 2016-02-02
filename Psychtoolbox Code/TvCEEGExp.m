function TvCEEGExp
% This function runs a threshold vs contrast psychophysical experiment
% All parameters are entered upon function launch
% Note that this experiment requires Psychtoolbox3 to run properly

clear;
if exist('onCleanup','class'), oC_Obj = onCleanup(@()sca);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GUI PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt={'Output Name:','Block Number:','% Gabor Contrasts:',...
    'Noise Contrasts:','Flicker Rates:'};
def={'S70','1','[0.1 0.1 0.1]','[0.1 0.1 0.1]','[8 12]'};
title='TvC Exp';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
trial.OutFileName = (char(answer(1)));
trial.BlockNumber = str2double(char(answer(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         TRIAL PROPERTIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1 for simulated responses, 0 for human subject
    simulation = 0;
    % 1 for on, 0 for off
    debug = 0;

    trial.randomseed = round(sum(clock)*10);
    rng(trial.randomseed,'twister');
    
    RefreshRate = 120;
    
    TotalTrials = 270;
    
    if debug == 1
        TotalTrials = 5;
    end
    
    % 1 = Low, 2 = Medium, 3 = High
    GaborContrastLevel = repmat([1 1 1 2 2 2 3 3 3],1,30);
    
    % 1 = Low, 2 = Medium, 3 = High
    NoiseContrastLevel = repmat([1 2 3 1 2 3 1 2 3],1,30);
    
    % Trial type marker
    TrialTypeMarker = repmat(1:9,1,30);
    
    trialmat = [GaborContrastLevel;NoiseContrastLevel;TrialTypeMarker];
    
    % Mix the trial properties matrix and split into "trial." structure
    shuffledtrialmat = shuffle2(trialmat,2);
    
    trial.GaborContrastLevel = zeros(1,TotalTrials);
    trial.NoiseContrastLevel = zeros(1,TotalTrials);
    trial.TrialTypeMarker = zeros(1,TotalTrials);
    
    for k = 1:TotalTrials
        trial.GaborContrastLevel(k) = shuffledtrialmat(1,k);
        trial.NoiseContrastLevel(k) = shuffledtrialmat(2,k);
        trial.TrialTypeMarker(k) = shuffledtrialmat(3,k);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   STIMULUS AND NOISE CONSTANTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    % Number of frames for each frequency cycle
    flickerframes = eval(char(answer(5)));
    
    trial.FlickerFrame(1) = flickerframes(1); % Gabor
    trial.FlickerFrame(2) = flickerframes(2); % Noise
    
    gaborcontrasts = eval(char(answer(3)));
    noisecontrasts = eval(char(answer(4)));
    if noisecontrasts > 1
        noisecontrasts = noisecontrasts/100;
    end
    if gaborcontrasts > 1
        gaborcontrasts = gaborcontrasts/100;
    end
    
    % GABOR INFO
    gaborangle = 0; % starting angle
    trial.GaborContrasts = gaborcontrasts;
    gabor_contrast(1) = gaborcontrasts(1);
    gabor_contrast(2) = gaborcontrasts(2);
    gabor_contrast(3) = gaborcontrasts(3);
    gabordim = 250; % Dimension of the gabor

    % NOISE INFO
    trial.NoiseContrasts = noisecontrasts;
    noise_contrast(1) = noisecontrasts(1);
    noise_contrast(2) = noisecontrasts(2);
    noise_contrast(3) = noisecontrasts(3);
    
    fixation_times = 700:50:1000;
    for n = 1:TotalTrials
        trial.fixation_frames(n) = fixation_times(randperm(7,1))/(1000/RefreshRate);
    end
    
    isi = 1000; %ms
    trial.isi = isi;
    shiftisi_frames = round(trial.isi/(1000/RefreshRate));
    postshift_frames = round(100/(1000/RefreshRate));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                      SETUP TRIAL STIM ARCHITECTURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    disp('Loading trial information...');
    for trialnum = 1:TotalTrials
        
        % Compute the number of frames the stimulus is on for each trial
        trial.NumStimFrames(trialnum) = shiftisi_frames+postshift_frames;
        
        % Determines on what frames the gabors need to rotate
        trial.shiftframe(trialnum) = shiftisi_frames+1;
        
        % Creates the vector for the reversing stimuli
        for m = 1:2
            temp = repmat([ones(1,trial.FlickerFrame(m)/2) ...
                zeros(1,trial.FlickerFrame(m)/2)],1,trial.NumStimFrames(trialnum)*2);
            trial.phaseshiftstream{trialnum}(m,:) = temp(1:trial.NumStimFrames(trialnum));
            % Creates the vector for the corresponding photocells
            if ismember(m,1:2)
                trial.photostream{trialnum}(m,:) = temp(1:trial.NumStimFrames(trialnum));
            end   
        end
        
        % Creates the angle matrix for each gabor frame stream
        targetStream = trial.phaseshiftstream{trialnum};
        anglestream = gaborangle*ones(1,length(targetStream));
        x = [-12 12];
        y = [-12 12];

        anglestream(trial.shiftframe(trialnum):end) = gaborangle+y(randi(length(y)));
 
        trial.TargetRotation(trialnum) = x(randi([1 2]));
        anglestream(trial.shiftframe(trialnum):end) = gaborangle+trial.TargetRotation(trialnum);
        
        % Gather all of the angle stream information into one cell struct
        trial.angleStream{trialnum} = anglestream;
        
        % Computes the indices when a new mesh is needed
        noiseshiftstream = trial.phaseshiftstream{trialnum}(2,:);
        noise_switch_ind{trialnum} = [1 strfind(noiseshiftstream,[0 1])+1];
    end
    disp('Loading complete. Press any key to begin...');
    
try
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     PSYCHTOOLBOX DEFAULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Color
    black = [0 0 0];
    white = [255 255 255];
    gray = [128 128 128];
    lightgray = [78 78 78];
   
    % Fix Java bs
    PsychJavaTrouble;
    
    % Display information
    displayGamma = 2.000; % Best to get the true value, rather than this
    ScreenDistance = 58; % cm
    ScreenHeight = 29.5;
    
    %Open a fullscreen window on the first screen with gray (~0.5) background
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UseVirtualFramebuffer');
    PsychImaging('AddTask','FinalFormatting','DisplayColorCorrection','SimpleGamma');
    whichscreen = max(Screen('screens'));
    [win,rect] = PsychImaging('OpenWindow', whichscreen,gray);
    PsychColorCorrection('SetEncodingGamma',win,1/displayGamma);
    PsychGPUControl('FullScreenWindowDisablesCompositor', 1);
    HideCursor;
    
    % Check that the refresh rate of the monitor is correct
    hz = Screen('NominalFrameRate', win);
    if hz ~= RefreshRate
        sca;
        disp('Refresh rate is not what you think it is!');
        return;
    end   
    
    % Display measurements
    ppd = pi/180 * ScreenDistance / ScreenHeight * rect(4); % pixels per degree
    StimCm = [7 7]; % stimulus size in cm
    m = 2*round(StimCm*ppd/2); % h and v stim size in pixel cm
    sf = 1/ppd; % cycles per pixel
    spatfreq = sf;
    
    trial.VisualAngle = rad2deg(2*atan(StimCm(2)/(2*ScreenDistance)));
    
    % Enable texture blending so the noise and gabors mix nicely
    
    % Unify key names across operating systems (just in case)
    KbName('UnifyKeyNames');
    
    % Load fonts
    myfont = '-bitstream-courier 10 pitch-bold-i-normal--0-0-0-0-m-0-ascii-0';
    
    % Set common key name values
    escapeKey = KbName('ESCAPE');
    SpaceBar = KbName('Space');
    left_shift = KbName('LEFTARROW');
    right_shift = KbName('RIGHTARROW');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   PREPARE STIM DISPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get stim screen resolution info
    [xcenter,ycenter] = RectCenterd(rect);
    
    % Set the stimulus and border rects
    srect = [0 0 gabordim gabordim];
    brect = [0 0 gabordim+5 gabordim+5];
    
    %%%%%%%%%%%%
    %          %
    %          %
    %        S %
    %%%%%%%%%%%%
    
    % Creates the rects for the photocells
    k = 0;
    photorect = [0 0 100 90];
    for m = 1:6
        pRect(m,:) = CenterRectOnPoint(photorect,50,50+k);
        k = k + 160;
    end
    
    % Offsets the rect to the appropriate corner
    stimrect = CenterRectOnPoint(srect, xcenter+200, ycenter+200);
    
    % Creates the rect for the box borders
    borderrect = CenterRectOnPoint(brect, xcenter+200, ycenter+200);
    
    %    Photocell x RGB
    %     0   1
    %     0   1
    %     0   1
    
    % Sets the photocells to black initially
    Screen('FillRect',win, black, pRect');
    
    % Flips the text and photocell rects to the screen
    Screen('Flip', win);
    
    % Draw experiment start text and black photocells
    Screen('TextSize', win, 28);
    Screen('TextFont', win, myfont);
    Screen('DrawText', win, 'Please wait for the experimenter to begin.', xcenter-(xcenter/2), ycenter, white);
    Screen('FillRect',win, black, pRect');
    Screen('Flip', win);
    
    % Waits until the space bar is pressed to begin the experiment
    while KbCheck;
    end;
    while true
        [resp, ~, keyCode] = KbCheck;
        
        if resp
            if keyCode(SpaceBar)
                break;
            elseif keyCode(escapeKey)
                sca;
                return;
            end
            while KbCheck;end
        end
    end
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                     PRIMARY EXPERIMENT TRIAL LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Priority(MaxPriority(win));
    
    for trialnum = 1:TotalTrials
        
        % Make each combined Gabor+Noise texture matrix for all frames of all trials
        ncont = noise_contrast(trial.NoiseContrastLevel(trialnum));
        gcont = gabor_contrast(trial.GaborContrastLevel(trialnum));
        for frame = 1:trial.NumStimFrames(trialnum)
            if ismember(frame,noise_switch_ind{trialnum})
                randmat{frame} = randn(64);
            else
                randmat{frame} = randmat{frame-1};
            end
            noisegabormat{frame} = MakeNoiseGabor(gcont,ncont,trial.angleStream{trialnum}(frame),...
                trial.phaseshiftstream{trialnum}(1,frame),...
                trial.phaseshiftstream{trialnum}(2,frame),...
                randmat{frame});
        end
        
        % Turn the Noise-Gabor matrix into an OpenGL texture
        for frame = 1:trial.NumStimFrames(trialnum)
            stimid(frame) = Screen('MakeTexture', win, noisegabormat{frame});
        end
        
        % Make a just noise texture to put up until the subject responds
        [~,endnoisemat] = MakeNoiseGabor(gcont,ncont,trial.angleStream{trialnum}(frame),...
                    trial.phaseshiftstream{trialnum}(1,frame),...
                    trial.phaseshiftstream{trialnum}(2,frame),...
                    randmat{frame});   

        endnoiseid = Screen('MakeTexture', win, endnoisemat+127);
                
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  DISPLAY SUBJECT REVIVAL TEXT  % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if trialnum == (TotalTrials/2)+1
                
            Screen('TextSize', win, 48);
            Screen('TextFont', win, myfont);
            Screen('DrawText', win, 'Yay you are 50% done! Wooo!', xcenter-250, ycenter-18, white);
            Screen('FillRect',win, black, pRect');
            Screen('Flip', win);
                    
            % Waits until the space bar is pressed to begin the experiment
            while KbCheck;
            end;
            while true
                [resp, ~, keyCode] = KbCheck;
                
                if resp
                    if keyCode(SpaceBar)
                        break;
                    elseif keyCode(escapeKey)
                        sca;
                        return;
                    end
                    while KbCheck;end
                end
            end
        
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  DISPLAY FIXATION CROSS AND BORDERS  % 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for pre_frame = 1:trial.fixation_frames(trialnum)
            Screen('TextSize', win, 48);
            Screen('TextFont', win, myfont);
            Screen('DrawText', win, '+', xcenter-16, ycenter-18, white);
            Screen('FrameRect',win, lightgray, borderrect',3);
            Screen('FillRect',win, black, pRect');
            Screen('Flip', win);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %              PRIMARY TRIAL STIMULUS MOVIE LOOP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        for frame = 1:trial.NumStimFrames(trialnum)
            
            Screen('FillRect',win, black, pRect(3,:));
            Screen('FillRect',win, black, pRect(4,:));
            Screen('FillRect',win, black, pRect(6,:));
            
            % Photocell spike for start of trial
            if ismember(frame,1)
                Screen('FillRect',win, white, pRect(5,:));
            else
                Screen('FillRect',win, black, pRect(5,:));
            end
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  DISPLAY FIXATION CROSS %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Screen('TextSize', win, 48);
            Screen('TextFont', win, myfont);
            Screen('DrawText', win, '+', xcenter-16, ycenter-18, white);
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %  DISPLAY BOX BORDERS %
            %%%%%%%%%%%%%%%%%%%%%%%%
            
            Screen('FrameRect',win, lightgray, borderrect',3);

            %%%%%%%%%%%%%%%%%%%%%
            %  DISPLAY STIMULUS %
            %%%%%%%%%%%%%%%%%%%%%
            
              Screen('DrawTexture',win,stimid(frame),[],stimrect,...
                  [],[],[],[],[],kPsychDontDoRotation);
              
            %%%%%%%%%%%%%%%%%%%%%%%
            %  DISPLAY PHOTOCELLS % 
            %%%%%%%%%%%%%%%%%%%%%%%
            
            % Gabor and Noise reversal spikes
            for m = 1:2
                if trial.photostream{trialnum}(m,frame) == 1
                    Screen('FillRect',win, white, pRect(m,:));
                else
                    Screen('FillRect',win, black, pRect(m,:));
                end
            end           
            
            % Photocell spike for the target rotation
            if ismember(frame,round(trial.shiftframe(trialnum))-1:round(trial.shiftframe(trialnum)))
                Screen('FillRect',win, white, pRect(6,:));
            else
                Screen('FillRect',win, black, pRect(6,:));
            end
            
            % Flip the noise and gabor texture layers to the screen
            Screen('Flip', win);
            
        end
        
        % Clear the texture to free up memory
        clear noisegabormat endnoisemat;
        
        Screen('DrawTextures',win,endnoiseid,[],stimrect,...
                  [],[],[],[],[],kPsychDontDoRotation);
        
        % Keep the border bars on until the subject responds
        Screen('FrameRect',win, lightgray, borderrect,3);
        
        % Trial Number Display
        Screen('TextSize', win, 48);
        Screen('TextFont', win, myfont);
        Screen('DrawText', win, num2str(trialnum), 20, 1000, black);
        
        % Keep photocells on until subject responds
        Screen('FillRect',win, black, pRect');

        % Keep fixation cross on until response
        Screen('TextSize', win, 48);
        Screen('TextFont', win, myfont);
        Screen('DrawText', win, '+', xcenter-16, ycenter-18, white);
        
        % Flip static textures to the screen
        vbl = Screen('Flip', win);
        
        % Records the trial duration
        trial.trialduration(trialnum) = toc;
        
        % Close open textures at the end of each trial
        Screen('Close');
        
        if simulation == 0;
            % Gets subject response for the target shift direction
            while KbCheck;end;
            
            % Waits until the space bar is pressed to continue
            % or the Esc key is pressed to close Screen
            while true
                [resp, ~, keyCode] = KbCheck;
                if resp
                    if keyCode(right_shift)
                        trial.response(trialnum) = -12;
                    elseif keyCode(left_shift)
                        trial.response(trialnum) = 12;
                    elseif keyCode(escapeKey)
                        % Save out the data
                        save(strcat(trial.OutFileName,'_TvCExp_',num2str(trial.BlockNumber)),'trial');
                        sca;
                        return;
                    end
                    
                    if keyCode(right_shift) || keyCode(left_shift)
                        if trial.response(trialnum) == trial.TargetRotation(trialnum)
                            trial.correct(trialnum) = 1;
                            break;
                        else
                            trial.correct(trialnum) = 0;
                            break;
                        end
                    end
                end
            end
            
            while KbCheck;end; % wait until all keys are released.
            
            % Waits until the space bar pressed to begin the experiment
            while true
                [resp, ~, keyCode] = KbCheck;
                
                if resp
                    if keyCode(SpaceBar)
                        break;
                    elseif keyCode(escapeKey)
                        sca;
                        return;
                    end
                    while KbCheck;end
                end
            end
        else
            trial.response(trialnum) = -12;
            if trial.response(trialnum) == trial.TargetRotation(trialnum)
                trial.correct(trialnum) = 1;
            else
                trial.correct(trialnum) = 0;
            end
            pause(1);
        end
        

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% Post-Movie Interval %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Draw experiment end text and black photocells
    Screen('TextSize', win, 28);
    Screen('TextFont', win, myfont);
    Screen('DrawText', win, 'Done! Please wait for the experimenter.', xcenter-(xcenter/2), ycenter, white);
    for m = 1:6
        Screen('FillRect',win, black, pRect(m,:));
    end
    Screen('Flip', win);
    
    % Save out the data
    save(strcat(trial.OutFileName,'_TvCExp_',num2str(trial.BlockNumber)),'trial');
    
    % Waits until the space bar is pressed to begin the experiment
    while KbCheck;
    end;
    while true
        [resp, ~, keyCode] = KbCheck;
        
        if resp
            if keyCode(SpaceBar)
                break;
            end
            while KbCheck;end
        end
    end
    
    sca
    clc
    
    disp(['Experiment block: ' num2str(trial.BlockNumber) ' complete. Behavioral data saved.']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
catch error
    ShowCursor;
    rethrow(error)
    save(strcat(trial.OutFileName,'_TvCExp_',num2str(trial.BlockNumber)),'trial');
    sca
end

end
