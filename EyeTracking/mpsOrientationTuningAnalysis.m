
function [r N bestfit] = mpsOrientationTuningAnalysis(subj,directory)
% Created: MPS - Aug 24 2010
%this function analyzes data from Orientation tuning. Give it two inputs:
%subj is the name of the file to be loaded, as a string, asks if no input
%directory is the name of the folder where the data is located, defaults to
%/desktop/OrientationTuning/psychophysics_data/ if no input
%% FYI

% Sequence = (randperm(40)+9); % set sequence order t for this trial %%%%
% runData.Sequence{RunNo,condition}(iOrder + sequencelength *(iTrial-1) + sequencelength*nTrialsEach*(iBlock-1)) = Sequence(iOrder);
% runData.response{RunNo,condition}(iSeq + sequencelength *(iTrial-1) + sequencelength*nTrialsEach*(iBlock-1)) = 1; 
%stimPhase = floor(Sequence(iSeq)/10)
%stimOrient = (Sequence(iSeq) - 10*stimPhase + 1)
%% Load data
if ~exist('subj','var'); subj = input('Please enter subject ID:  ','s'); end
if ~exist('directory','var'); directory = '~/Desktop/OrientationTuning/psychophysics_data'; end

dataFile = fullfile(directory,[subj]);
load(dataFile);
%% Set task variables if they don't exist
if ~isfield(runData,'nRuns')
    runData.nRuns = 1;
end
if ~isfield(runData,'nTrialsEach')
    runData.nTrialsEach = 45;
end
if ~isfield(runData,'nBlocks')
    runData.nBlocks = 15;
end
if ~isfield(runData,'targetOrients')
    runData.targetOrients = [0:pi/10:9*pi/10];
end
if ~isfield(runData,'conditionOrder')
    runData.conditionOrder = randperm(3);
end
if ~isfield(runData,'phases')
    runData.phases = 1:4;
end
if ~isfield(runData,'sequencelength')
    runData.sequencelength = 40;
end
if ~isfield(runData,'exOrient'); runData.exOrient = input('Enter the target orientation(s) as a vector: '); end
if ~isfield(runData,'speed'); runData.speed = input('Enter the stimulus presentation speed: '); end

Tsize = 60/runData.speed; %used to make the window of time below
%% Check out our data

for T = 1:(Tsize +1) %this will be the window of time (one second) that we will look back before the response, want 0-30 (@ 30 hz) but must start with a 1, so subtract later
    for aRun = 1:runData.nRuns
        for aCon = 1:length(runData.conditionOrder)
            iCount(aRun,aCon) = 0;
            for theta = 1:length(runData.targetOrients) %Create a set of counters for all 10 orientations (theta) with a range of 0-Tzise frames (1 sec) 
            r{aCon}(T,theta) = 0;
            end
            for aBlock = [1:4,6:9,11:14] %1:runData.nBlocks - don't use blocks in which Bullseye was presented
            for aTrial = 1:runData.nTrialsEach
                for aSeq = 1:runData.sequencelength
                    if (aSeq + runData.sequencelength *(aTrial-1) + runData.sequencelength*runData.nTrialsEach*(aBlock-1)) <= size(runData.response{aRun,runData.conditionOrder(aCon)}(:,:),2) 
                        Lookup = (aSeq + runData.sequencelength *(aTrial-1) + runData.sequencelength*runData.nTrialsEach*(aBlock-1));
                        if runData.response{aRun,runData.conditionOrder(aCon)}(Lookup) == 1;
                            
                            if aTrial == 1 && aBlock == 1
                                if aSeq - (T - 1) > 0 % can't go back before the beginning!
                                    r{aCon}(T,(runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1)) ...
                                    = (r{aCon}(T,((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1))) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1))+1); %this is a counter which is incremented for each response, looks up the orientation T frames back, and averages across phase
                                    iCount(aRun,aCon) = iCount(aRun,aCon) + 1; %keep track of the number of responses
                                end
                            else
                                r{aCon}(T,((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1))) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1)) ...
                                = (r{aCon}(T,((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1))) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1))+1);
                                iCount(aRun,aCon) = iCount(aRun,aCon) + 1; %keep track of the number of responses
                            end
                        end;end
                end
            end; end
        end
    end
end

for aRun = 1:runData.nRuns
    for aCon = 1:length(runData.conditionOrder)
        for T = 1:(Tsize +1)
            for theta = 1:length(runData.targetOrients)
               r{aCon}(T,theta) = (r{aCon}(T,theta)/iCount(aRun,aCon)); %normalize by total number of responses for that condition to obtain a probability that each orientation (theta) was present T images before the response
            end
        end
        Rmax(aCon) = max(r{aCon}(:,runData.exOrient(runData.conditionOrder(aCon)))); %find the max of r --> this is T'
        Imax(aCon) = find(r{aCon}(:,runData.exOrient(runData.conditionOrder(aCon))) == Rmax(aCon),1); %find the index
        if Imax(aCon) == 1
            Rchunk{aCon}(1,:) = r{aCon}(Imax(aCon),:);
            Rchunk{aCon}(2,:) = r{aCon}(Imax(aCon) +1,:);
        elseif Imax(aCon) == 31
            Rchunk{aCon}(1,:) = r{aCon}(Imax(aCon) - 1,:);
            Rchunk{aCon}(2,:) = r{aCon}(Imax(aCon),:);
        else
        Rchunk{aCon}(1,:) = r{aCon}(Imax(aCon) - 1,:);%these next three lines grab a chunk of the data corresponding to T'-1, T', T'+1 --> similar amount of time to +/- 1 SD in RT from Ringach
        Rchunk{aCon}(2,:) = r{aCon}(Imax(aCon),:);
        Rchunk{aCon}(3,:) = r{aCon}(Imax(aCon) +1,:);
        end
        Rhat{aCon} = mean(Rchunk{aCon},1);%this takes the average of the chunk
        for theta = 1:length(runData.targetOrients)
            N{aCon}(theta) = ((Rhat{aCon}(theta) - min(Rhat{aCon}))/(max(Rhat{aCon}) - min(Rhat{aCon}))); %this will normalize the data between 0 and 1
        end
        if runData.exOrient(runData.conditionOrder(aCon)) == 1  %this section re-orders data to display it as a centered "mexican hat"
            for theta = 1:length(runData.targetOrients)/2
                M(theta) = N{aCon}(theta + length(runData.targetOrients)/2);
            end
            for theta = (length(runData.targetOrients)/2 +1):length(runData.targetOrients)
                M(theta) = N{aCon}(theta - length(runData.targetOrients)/2);
            end
            N{aCon} = M;
        end
        N{aCon}(11) = N{aCon}(1);
    end
end

%% BE data analysis

for T = 1:(Tsize +1) %this will be the window of time (one second) that we will look back before the response, want 0-30 (@ 30 hz) but must start with a 1, so subtract later
    for aRun = 1:runData.nRuns
        for aCon = 1:length(runData.conditionOrder)
            iCountBE(aRun,aCon) = 0;
            for theta = 1:length(runData.targetOrients) %Create a set of counters for all 10 orientations (theta) with a range of 0-Tzise frames (1 sec) 
            rBE{aCon}(T,theta) = 0;
            end
            for aBlock = [5,10,15] %1:runData.nBlocks - don't use blocks in which Bullseye was presented
            for aTrial = 1:runData.nTrialsEach
                for aSeq = 1:runData.sequencelength
                    if (aSeq + runData.sequencelength *(aTrial-1) + runData.sequencelength*runData.nTrialsEach*(aBlock-1)) <= size(runData.BEresponse{aRun,runData.conditionOrder(aCon)}(:,:),2) 
                        Lookup = (aSeq + runData.sequencelength *(aTrial-1) + runData.sequencelength*runData.nTrialsEach*(aBlock-1));
                        if runData.BEresponse{aRun,runData.conditionOrder(aCon)}(Lookup) == 1;
                            
                            if aTrial == 1 && aBlock == 1
                                if aSeq - (T - 1) > 0 % can't go back before the beginning!
                                    rBE{aCon}(T,(runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1)) ...
                                    = (rBE{aCon}(T,((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1))) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1))+1); %this is a counter which is incremented for each response, looks up the orientation T frames back, and averages across phase
                                    iCountBE(aRun,aCon) = iCountBE(aRun,aCon) + 1; %keep track of the number of responses
                                end
                            else
                                rBE{aCon}(T,((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1))) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1)) ...
                                = (rBE{aCon}(T,((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1))) - (floor((runData.Sequence{aRun,runData.conditionOrder(aCon)}(Lookup - (T - 1)))/10)*10) + 1))+1);
                                iCountBE(aRun,aCon) = iCountBE(aRun,aCon) + 1; %keep track of the number of responses
                            end
                        end;end
                end
            end; end
        end
    end
end

for aRun = 1:runData.nRuns
    for aCon = 1:length(runData.conditionOrder)
        for T = 1:(Tsize +1)
            for theta = 1:length(runData.targetOrients)
               rBE{aCon}(T,theta) = (rBE{aCon}(T,theta)/iCountBE(aRun,aCon)); %normalize by total number of responses for that condition to obtain a probability that each orientation (theta) was present T images before the response
            end
        end
        RBEmax(aCon) = max(rBE{aCon}(:,3)); %find the max of r --> this is T'
        IBEmax(aCon) = find(rBE{aCon}(:,3) == RBEmax(aCon),1); %find the index
        if IBEmax(aCon) == 1
            RBEchunk{aCon}(1,:) = rBE{aCon}(IBEmax(aCon),:);
            RBEchunk{aCon}(2,:) = rBE{aCon}(IBEmax(aCon) +1,:);
        elseif IBEmax(aCon) == 31
            RBEchunk{aCon}(1,:) = rBE{aCon}(IBEmax(aCon) - 1,:);
            RBEchunk{aCon}(2,:) = rBE{aCon}(IBEmax(aCon),:);
        else
        RBEchunk{aCon}(1,:) = rBE{aCon}(IBEmax(aCon) - 1,:);%these next three lines grab a chunk of the data corresponding to T'-1, T', T'+1 --> similar amount of time to +/- 1 SD in RT from Ringach
        RBEchunk{aCon}(2,:) = rBE{aCon}(IBEmax(aCon),:);
        RBEchunk{aCon}(3,:) = rBE{aCon}(IBEmax(aCon) +1,:);
        end
        RBEhat{aCon} = mean(RBEchunk{aCon},1);%this takes the average of the chunk
        for theta = 1:length(runData.targetOrients)
            NBE{aCon}(theta) = ((RBEhat{aCon}(theta) - min(RBEhat{aCon}))/(max(RBEhat{aCon}) - min(RBEhat{aCon}))); %this will normalize the data between 0 and 1
        end
        %if runData.exOrient(runData.conditionOrder(aCon)) == 1  %this section re-orders data to display it as a centered "mexican hat"
%             for theta = 1:length(runData.targetOrients)/2
%                 MBE(theta) = NBE{aCon}(theta + length(runData.targetOrients)/2);
%             end
%             for theta = (length(runData.targetOrients)/2 +1):length(runData.targetOrients)
%                 MBE(theta) = NBE{aCon}(theta - length(runData.targetOrients)/2);
%             end
%             NBE{aCon} = MBE;
        %end
        NBE{aCon}(11) = NBE{aCon}(1);
    end
end
BE_SNR = RBEhat{1}(1,3)/mean(std(RBEchunk{1}(:,[1,2,4:10])));
%% fit and make some graphs
time = 0:1:10;

hold on
for aCon = 1:length(runData.conditionOrder)
%      Min1 = find(N{aCon} == min(N{aCon}(2:6)));
%      Min2 = find(N{aCon} == min(N{aCon}(6:10)));
%      MeanInh{aCon} = mean([N{aCon}(1:Min1) N{aCon}(Min2:10)]);
%      MeanEx{aCon} = mean(N{aCon}(Min1:Min2));
%      for theta1 = Min1:Min2
%          DiffEx{aCon}(theta1) = (N{aCon}(theta1) - MeanEx{aCon})^2; end
%      for theta2 = [1:Min1 Min2:10]
%          DiffInh{aCon}(theta2) = (N{aCon}(theta2) - MeanInh{aCon})^2; end
%      stdevEx(aCon) = sqrt(1/length(theta1)*sum(DiffEx{aCon}(1:end)));
%      stdevInh(aCon) = sqrt(1/(theta2)*sum(DiffInh{aCon}(1:end)));

if Rhat{aCon}(1,runData.exOrient(runData.conditionOrder(aCon))) < 2.5*mean(std(r{aCon}(:,find([1:10]~=runData.exOrient(runData.conditionOrder(aCon))))))+.1 ...
        || Rhat{aCon}(1,runData.exOrient(runData.conditionOrder(aCon))) < .95*max(max(r{aCon}(:,find([1:10]~=runData.exOrient(runData.conditionOrder(aCon)))))); % Rhat peak needs to be greater than 2 Std dev (of orientations other than target) higher than the mean, and greater than 90% of the peak value in r
no_fit = 1;
bestfit = 0;
break
else no_fit = 0;
end
     
%      minfit(aCon,:) = [1,stdevEx(aCon)-2,-.2,-1,stdevInh(aCon)-2,-.2];
%      maxfit(aCon,:) = [2,stdevEx(aCon)+2,.2,0,stdevInh(aCon)+2,.2];
     minfit(aCon,:) = [0 0 -.0001 -2 0 -.0001]; %make no assumptions about the shap of Gaussians
     maxfit(aCon,:) = [2 Inf .0001 0 Inf .0001];

    figure
    subplot(length(runData.conditionOrder),2,(2*aCon-1))
    plot(N{aCon}(1:11))
    eval(['titlename = runData.dataLabel' num2str(runData.conditionOrder(aCon)) ';']);
    title(titlename);
    xlabel('Orientation');
    ylabel('Normalized Response');
    axis tight;
    subplot(length(runData.conditionOrder),2,(2*aCon))
    imagesc(r{aCon}(:,:))
    xlabel('Orientation');
    ylabel('Time (frames)');
    title(subj(1:8));
    set(gcf,'POS',[5 513 560 420]);
    figure
    bestfit{aCon} = mpsFitDOG(N{aCon},time,minfit(aCon,:),maxfit(aCon,:),1);
    set(gcf,'POS',[5 47 560 420]);
end
hold off

%% BE plot
figure
hold on
    subplot(length(runData.conditionOrder),2,(2*aCon-1))
    bar(NBE{aCon}(1:11))
    bar(3,NBE{aCon}(3),'r')
    %eval(['titlename = runData.dataLabel' num2str(runData.conditionOrder(aCon)) ';']);
    title(['Bull''s Eye Control SNR = ' num2str(BE_SNR)]);
    xlabel('Orientation, 3 = Target');
    ylabel('Normalized Response');
    axis tight;
    subplot(length(runData.conditionOrder),2,(2*aCon))
    imagesc(rBE{aCon}(:,:))
    xlabel('Orientation, 3 = Target');
    ylabel('Time (frames)');
    %title(subj(1:8));
    set(gcf,'POS',[5 513 560 420]);
hold off
%% Save
if no_fit ==1
saveFile = fullfile(directory,[subj(1:end-4) '_Analyzed_noFit.mat']);
else
saveFile = fullfile(directory,[subj(1:end-4) '_Analyzed.mat']);
end
save(saveFile,'N','r','runData','bestfit');
