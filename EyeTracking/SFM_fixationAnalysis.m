% Funciton to load in edf converted files and perform fixation based analysis.

function [options,data] = SFM_fixationAnalysis(options)

%
%%
% addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode/'));
% addpath(genpath('E:/GitRepos/SFM.git'));
addpath(genpath('/home/shaw-raid1/data/psychophysics/SFM.git'));
addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode'));
% addpath(genpath('E:/GitRepos/edf-converter.git'));

%% opt
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'top_dir')
    options.top_dir = '/home/shaw-raid1/data/pHCP.git/subjectResponseFiles/';
end
if ~isfield(options,'dateCutoff')
    options.dateCutoff = 1;   % Only take participants after a given date
    options.dateCutoffVal = 20211102;   % Don't take ppt after this date
    %     options.dateCutoffVal = 20210815;
    options.dateCutoffValConverted = datenum(num2str(options.dateCutoffVal),'yyyymmdd');
end
if ~isfield(options,'excludeRedcap')
    options.excludeRedcap = 1; % exclude subjects based on control task performance, 0 = no, 1 = yes
    addpath(genpath('/home/shaw-raid1/matlab_tools/COP_analysis.git')) % Add path for redcap excludion function (in COP_analysis.git)
end
if ~isfield(options,'excludeTypeA')
    options.excludeTypeA = 1; % exclude subjects based on control task performance, 0 = no, 1 = yes
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 0; % 1 = on, 0 = off
end
if ~isfield(options,'displayETFigs')
    options.displayETFigs = 1; % 1 = on, 0 = off
end
if ~isfield('options','buttonEventTesting')
    options.buttonEventTesting = 0;
end

%% Load in the data
% Grab the behavioral data
options.curDur = '/home/shaw-raid1/data/psychophysics/SFM.git/';
cd(options.curDur)

% pull in SFM data
if options.excludeTypeA
    [dataAve] = SFM_Behav_Control_Group_Analysis(options);
    illusoryCutoff = dataAve.illusoryCutoff;
    [~,illusoryCutoffIdx] = unique(illusoryCutoff(:,1));
    illusoryCutoff = illusoryCutoff(illusoryCutoffIdx,:);
    options.illusoryCutoff = illusoryCutoff;
    SFMdata = analyze_SFM_data(options);
    clear dataAve illusoryCutoff illusoryCutoff illusoryCutoffIdx
else
    SFMdata = analyze_SFM_data(options);
end

counter = 0;
% Remove subjects who have no responses in their B visit A run or all 5
% blocks of B run
for iI=1:length(SFMdata.flip_time(:,1,1,1))
    if ~isempty(SFMdata.flip_time{iI,1,1,1}) &&...
            ~(isempty(SFMdata.flip_time{iI,1,2,1}) &&...
            isempty(SFMdata.flip_time{iI,1,2,2}) &&...
            isempty(SFMdata.flip_time{iI,1,2,3}) &&...
            isempty(SFMdata.flip_time{iI,1,2,4}) &&...
            isempty(SFMdata.flip_time{iI,1,2,5}))
        counter = counter+1;
        flipTimeHolder(counter,:,:,:) = SFMdata.flip_time(iI,:,:,:);
        behavSubjHolder(counter) = SFMdata.subj_number(iI);
        behavDateHolder(counter,:) = SFMdata.date_num(iI,:);
        switchRateHolder(counter,:,:,:) = SFMdata.Hz_flips(iI,:,:,:);
    end
end

% Grab the subj numbers, dates, and flip times for behavioral data
for iI = 1:size(flipTimeHolder,1)   % Num subjs
    % flip_times = (subj x B,Z x A,B x block)
    for iJ = 1:size(flipTimeHolder,2)   % B/Z
        data.behav.A.flipTimes{iI,iJ}{1} = flipTimeHolder{iI,iJ,1,1} .* 1000;
        data.behav.A.flipTimesIdx{iI,iJ}{1} = data.behav.A.flipTimes{iI,iJ}{1} ./ 2;
        for iK = 1:size(flipTimeHolder,4)   % Block num
            data.behav.B.flipTimes{iI,iJ}{iK} = flipTimeHolder{iI,iJ,2,iK} .* 1000;
            data.behav.B.flipTimesIdx{iI,iJ}{iK} = data.behav.B.flipTimes{iI,iJ}{iK} ./ 2;
        end
    end
end
data.behav.A.Hz_flips = squeeze(switchRateHolder(:,:,1,1));
data.behav.B.Hz_flips = squeeze(switchRateHolder(:,:,2,:));
data.behav.subjID = behavSubjHolder;
data.behav.dateNum = behavDateHolder;
data.behav.date(:,1) = [str2num(datestr(data.behav.dateNum(:,1),'yyyymmdd'))];
data.behav.date(~isnan(data.behav.dateNum(:,2)),2)=str2num(datestr(data.behav.dateNum(~isnan(data.behav.dateNum(:,2)),2),'yyyymmdd'));
data.behav.date(data.behav.date(:,2)==0,2) = NaN;
clear SFMdata counter behavSubjHolder behavDateHolder flipTimeHolder

% Load in the control switch values and response data
data.behav.controlSwitchTimesOrig = readtable('SFMControlSwitches.xlsx');
data.behav.controlSwitchTimesOrig = table2cell(data.behav.controlSwitchTimesOrig);
data.behav.controlSwitchNames = {data.behav.controlSwitchTimesOrig{4:end,1}};
data.behav.controlSwitchTimes = [data.behav.controlSwitchTimesOrig{4:end,2}] .* 1000;
data.behav.controlSwitchInd = data.behav.controlSwitchTimes ./ 2;

% Load in participant file names
cd ./EyeTracking/
options.curDur = pwd;
% options.etDataDur = 'Z:/eyetracking_data/psychophysics_eyetracking/';
options.etDataDur = '/home/shaw-raid1/data/psychophysics/eyetracking';
cd(options.etDataDur)
options.fileNames = dir('P*');

% Make new list that uses dates from behavioral data to determine B/Z
% visits for each ET data file for all subjects.
% Will want arrays of ET data files structured like:
% data.et.A.etDataOrig{subjNum x B/Z}
% Then for all other data:
% data.et.A.EXAMPLE_DATA_STRUCT{subjNum x B/Z}{runNum}(....)
% Where runNum will = 1 in A and = 5 in B

% For each participant load in the .mat files
counterA = 0;
counterB = 0;
for iI = 1:length(options.fileNames)
    
    fprintf('%s%d%s%s\n','Loading participant ',iI,': ',options.fileNames(iI).name)
    
    % CD in participant dir
    cd(['./' options.fileNames(iI).name])
    options.dirA{iI} = dir('*_A_*_Useable.mat');
    options.dirB{iI} = dir('*_B_*_Useable.mat');
    
    % Load A
    if isempty(options.dirA{iI})
        warning(['No ''A'' data for subj: ' options.fileNames(iI).name]);
    else
        % If this subj is not in the behavioral subjid list then exclude
        % Find subj Id for this et file
        % If there are 2 dirs grab the first one
        % (Since subjID will be the same for both)
        holderName = {options.dirA{iI}.name};
        startInd = strfind(holderName{1},'P')+1;   % char after the 'P'
        endInd = strfind(holderName{1},'r')-2;   % char before the '_run'
        subjNumHolder = str2num(holderName{1}(startInd:endInd));
        subjExcludeSwitch = find(data.behav.subjID == subjNumHolder);
        
        if ~isempty(subjExcludeSwitch)
            counterA = counterA+1;
            clear etDataUseable
            
            % Record the subj number for this et data file
            data.et.A.subjID(counterA,1) = subjNumHolder;
            clear holderName startInd endInd subjNumHolder subjExcludeSwitch
            
            for iJ=1:length(options.dirA{iI})
                % Record the date of this file
                startInd = strfind(options.dirA{iI}(iJ).name,'A')+2;   % char after the 'A_'
                endInd = strfind(options.dirA{iI}(iJ).name,'U')-7;   % char before the '_TIME_U'
                holderName = ...
                    str2num(options.dirA{iI}(iJ).name(startInd:endInd));
                
                % For this subj, determine if this is a B/Z visit (1st or 2nd)
                % using behavioral data. Find correct index or date, 1st or 2nd
                % position in the date array.
                % First find the index of the subject in the behavioral array (presumably,
                % if the subject has ET data, they also have behavioral data..)
                partIdx = find(data.behav.subjID ==...
                    data.et.A.subjID(counterA));
                % Then find the index of the correct date
                dateIdx = find(data.behav.date(partIdx,:) ==...
                    holderName);
                % Then assign the date to the correct place in the date array
                data.et.A.date(counterA,dateIdx) =...
                    holderName;
                
                % Load in the data in the correct position
                load(options.dirA{iI}(iJ).name)
                data.et.A.etDataOrig{counterA,dateIdx} = etDataUseable;
                
                clear holderName startInd endInd dateIdx partIdx etDataUseable
            end
        end
    end
    
    % Load B
    if isempty(options.dirB{iI})
        warning(['No ''B'' data for subj: ' options.fileNames(iI).name]);
    else
        % If this subj is not in the behavioral subjid list then exclude
        % Find subj Id for this et file
        % If there are 2 dirs grab the first one
        % (Since subjID will be the same for both)
        holderName = {options.dirB{iI}.name};
        startInd = strfind(holderName{1},'P')+1;   % char after the 'P'
        endInd = strfind(holderName{1},'r')-2;   % char before the '_run'
        subjNumHolder = str2num(holderName{1}(startInd:endInd));
        subjExcludeSwitch = find(data.behav.subjID == subjNumHolder);
        
        if ~isempty(subjExcludeSwitch)
            counterB = counterB+1;
            clear etDataUseable
            
            % Record the subj number for this et data file
            data.et.B.subjID(counterB,1) = subjNumHolder;
            clear holderName startInd endInd subjNumHolder subjExcludeSwitch
            
            for iJ=1:length(options.dirB{iI})
                % Record the date of this file
                startInd = strfind(options.dirB{iI}(iJ).name,'B')+2;   % char after the 'A_'
                endInd = strfind(options.dirB{iI}(iJ).name,'U')-7;   % char before the '_TIME_U'
                holderName = ...
                    str2num(options.dirB{iI}(iJ).name(startInd:endInd));
                
                % For this subj, determine if this is a B/Z visit (1st or 2nd)
                % using behavioral data. Find correct index or date, 1st or 2nd
                % position in the date array.
                % First find the index of the subject in the behavioral array (presumably,
                % if the subject has ET data, they also have behavioral data..)
                partIdx = find(data.behav.subjID ==...
                    data.et.B.subjID(counterB));
                % Then find the index of the correct date
                dateIdx = find(data.behav.date(partIdx,:) ==...
                    holderName);
                % Then assign the date to the correct place in the date array
                data.et.B.date(counterB,dateIdx) =...
                    holderName;
                
                % Load in the data in the correct position
                load(options.dirB{iI}(iJ).name)
                data.et.B.etDataOrig{counterB,dateIdx} = etDataUseable;
                
                clear holderName startInd endInd dateIdx partIdx etDataUseable
            end
        end
    end
    cd ../
end
clear counterA counterB

% CD back to ET analysis folder for saving later
cd(options.curDur)

% Change date 0's to NaNs
data.et.A.date(data.et.A.date==0) = NaN;
data.et.B.date(data.et.B.date==0) = NaN;

% Now update the behavioral data to exclude subjects not found in the ET
% data list.
counter = 0;
for iI=1:length(data.behav.subjID)
    if ~isempty(find(data.et.A.subjID == data.behav.subjID(iI)))
        counter = counter+1;
        
        subjIDHolder(counter) = data.behav.subjID(iI);
        dateNumHolder(counter,:) = data.behav.dateNum(iI,:);
        dateHolder(counter,:) = data.behav.date(iI,:);
        
        switchRateHolderA(counter,:) = data.behav.A.Hz_flips(iI,:);
        switchRateHolderB(counter,:,:) = data.behav.B.Hz_flips(iI,:,:);
        
        AflipTimes(counter,:) = data.behav.A.flipTimes(iI,:);
        BflipTimes(counter,:) = data.behav.B.flipTimes(iI,:);
        AflipTimesIdx(counter,:) = data.behav.A.flipTimesIdx(iI,:);
        BflipTimesIdx(counter,:) = data.behav.B.flipTimesIdx(iI,:);
    else
    end
end
data.behav.subjID = subjIDHolder';
data.behav.dateNum = dateNumHolder;
data.behav.date = dateHolder;

data.behav.A.flipTimes = AflipTimes;
data.behav.B.flipTimes = BflipTimes;
data.behav.A.flipTimesIdx = AflipTimesIdx;
data.behav.B.flipTimesIdx = BflipTimesIdx;

data.behav.A.Hz_flips = switchRateHolderA;
data.behav.B.Hz_flips = switchRateHolderB;

clear subjIDHolder dateNumHolder dateHolder AflipTimes BflipTimes...
    AflipTimesIdx BflipTimesIdx switchRateHolderA switchRateHolderB

% Find the total number of subjects in each group
data.behav.grpIdx{1} = find(data.behav.subjID < 2000000);
data.behav.grpIdx{2} = find(data.behav.subjID > 2000000 & data.behav.subjID < 6000000);
data.behav.grpIdx{3} = find(data.behav.subjID > 6000000);

data.et.A.grpIdx{1} = find(data.et.A.subjID < 2000000);
data.et.A.grpIdx{2} = find(data.et.A.subjID > 2000000 & data.et.A.subjID < 6000000);
data.et.A.grpIdx{3} = find(data.et.A.subjID > 6000000);
data.et.B.grpIdx{1} = find(data.et.B.subjID < 2000000);
data.et.B.grpIdx{2} = find(data.et.B.subjID > 2000000 & data.et.B.subjID < 6000000);
data.et.B.grpIdx{3} = find(data.et.B.subjID > 6000000);

%% Pull out and organize the relevant data from the raw ET data
% Params
% Screen
options.screenSize = [1 1 1920 1080];

% Run
options.runType = {'A','B'};
options.blockNumA = 1;
options.blockNumB = 5;

% Calculat PPD value
options.mon_width_cm = 53;   % Width of the monitor (cm)
options.mon_dist_cm = 70;   % Viewing distance (cm)
options.mon_width_deg = 2 * (180/pi) * atan((options.mon_width_cm/2)/options.mon_dist_cm);   % Monitor width in DoVA
options.PPD = (options.screenSize(3)/options.mon_width_deg);   % pixels per degree

% First grab the ET data for each block for each participant
for iB=1:2   % A/B run
    for iI=1:size(data.et.(options.runType{iB}).etDataOrig,1)   % Num participants
        for iZ=1:size(data.et.(options.runType{iB}).etDataOrig,2)   % B/Z day
            % Only look at participats/runs that have ET data
            if isnan(data.et.(options.runType{iB}).date(iI,iZ))
                warning('%s%s%s%d%s%d','No fixation data for ',(options.runType{iB}),' run for participant: ',...
                    data.et.(options.runType{iB}).subjID(iI,1),', day: ', iZ)
            else
                fprintf('%s%s%s%d%s%d\n','Grabbing fixation data for ',(options.runType{iB}),' run for participant: ',...
                    data.et.(options.runType{iB}).subjID(iI,1),', day: ', iZ)
                
                %% Find events
                % First find all start and end events
                data.et.(options.runType{iB}).blockStartIdx{iI,iZ} = find(~cellfun(@isempty,...
                    strfind(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.info,...
                    'start block')));
                data.et.(options.runType{iB}).blockStartType{iI,iZ} = ...
                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.info(...
                    data.et.(options.runType{iB}).blockStartIdx{iI,iZ});
                
                data.et.(options.runType{iB}).blockEndIdx{iI,iZ} = find(~cellfun(@isempty,...
                    strfind(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.info,...
                    'end block')));
                data.et.(options.runType{iB}).blockEndType{iI,iZ} = ...
                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.info(...
                    data.et.(options.runType{iB}).blockEndIdx{iI,iZ});
                
                % Make sure the number of start/end event codes match.
                % If there's an extra start event w/ no corresponding end,
                % find the end time point, and set values based on that. - KWK 20211109
                % Create array of found block start and end events
                data.et.(options.runType{iB}).blockEventsFound{iI,iZ} =...
                    zeros([2,options.(['blockNum', options.runType{iB}])]);
                % Grab the event numbers from the event strings
                data.et.(options.runType{iB}).blockStartNum{iI,iZ} =...
                    cellfun(@str2num,cellfun(@(x) x(end),{data.et.(options.runType{iB}).blockStartType{iI,iZ}{:}},...
                    'UniformOutput',false));
                data.et.(options.runType{iB}).blockEndNum{iI,iZ} =...
                    cellfun(@str2num,cellfun(@(x) x(end),{data.et.(options.runType{iB}).blockEndType{iI,iZ}{:}},...
                    'UniformOutput',false));
                
                % First check whether or not there are the CORRECT number
                % of events. Then check if therer is one start/end event
                % missing. If so we can find it's matching event using
                % timing. If not, and we're missing an end and start, throw
                % error, b/c we are missing an entire block.
                if ~((numel(data.et.(options.runType{iB}).blockStartNum{iI,iZ}) ==...
                        options.(['blockNum', options.runType{iB}])) && ...
                        (numel(data.et.(options.runType{iB}).blockEndNum{iI,iZ}) ==...
                        options.(['blockNum', options.runType{iB}])))
                    
                    % Check that there are the same number of start and end events
                    if ~(numel(data.et.(options.runType{iB}).blockStartNum{iI,iZ}) ==...
                            numel(data.et.(options.runType{iB}).blockEndNum{iI,iZ}))
                        % Keep track of runs that don't have matching start/end events
                        % Find the missing event
                        counter1=0;
                        counter2=0;
                        clear blockIdxHolder
                        for iJ=1:options.(['blockNum', options.runType{iB}])
                            % For start events
                            if sum(data.et.(options.runType{iB}).blockStartNum{iI,iZ}+1==iJ)==1
                                counter1=counter1+1;
                                data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(1,iJ) = 1;
                                blockIdxHolder(1,iJ) = ...
                                    data.et.(options.runType{iB}).blockStartIdx{iI,iZ}(counter1);
                            else
                                data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(1,iJ) = 0;
                                blockIdxHolder(1,iJ) = 0;
                            end
                            % For end events
                            if sum(data.et.(options.runType{iB}).blockEndNum{iI,iZ}+1==iJ)==1
                                counter2=counter2+1;
                                data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(2,iJ) = 1;
                                blockIdxHolder(2,iJ) = ...
                                    data.et.(options.runType{iB}).blockEndIdx{iI,iZ}(counter2);
                            else
                                data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(2,iJ) = 0;
                                blockIdxHolder(2,iJ) = 0;
                            end
                        end
                    else
                        % If blocks are the same size, but not = correct
                        % num elements for that block, throw an error.
                        error('%s%s%s%d%s%d','Incorrect number of blocks found for ',(options.runType{iB}),' run for participant: ',...
                            data.et.(options.runType{iB}).subjID(iI,1),', day: ', iZ)
                    end
                else
                    % If there are the corr number of blocks in each run
                    data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(:,:) = 1;
                end
                clear counter1 counter2
                
                % For each block in the run
                for iJ=1:length(data.et.(options.runType{iB}).blockStartIdx{iI,iZ})
                    
                    fprintf('%s%d\n','     Processing block: ',iJ)
                    
                    %% Grab start end points using event codes
                    % If we have the correct number of start/end events for
                    % this run, iZ, load event times/data in normally.
                    if sum(data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(1,:))==...
                            options.(['blockNum', options.runType{iB}]) &&...
                            sum(data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(2,:))==...
                            options.(['blockNum', options.runType{iB}])
                        % Sometimes start event doesn't match up with any sample event time
                        % point. Find the closest event to the start event. KWK - 20211026
                        [~, closest] = min(...
                            abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                            data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.time(...
                            data.et.(options.runType{iB}).blockStartIdx{iI,iZ}(iJ))));
                        data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) =...
                            data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                        clear closest
                        % Seems like end event time is one time point (4 ms) beyond the last
                        % sample so fdin the closest time point to get the idx for last
                        % event - KWK 20211026
                        [~, closest] = min(...
                            abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                            data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.time(...
                            data.et.(options.runType{iB}).blockEndIdx{iI,iZ}(iJ))));
                        data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ) =...
                            data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                        clear closest
                    else
                        if data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(1,iJ)==1 &&...
                                data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(2,iJ)==1
                            % Sometimes start event doesn't match up with any sample event time
                            % point. Find the closest event to the start event. KWK - 20211026
                            [~, closest] = min(...
                                abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                                data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.time(...
                                blockIdxHolder(1,iJ))));
                            data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) =...
                                data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                            clear closest
                            % Seems like end event time is one time point (4 ms) beyond the last
                            % sample so fdin the closest time point to get the idx for last
                            % event - KWK 20211026
                            [~, closest] = min(...
                                abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                                data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.time(...
                                blockIdxHolder(2,iJ))));
                            data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ) =...
                                data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                            clear closest
                        else
                            % If we have a missing start or end for this run,
                            % then find the corresponding time value.
                            blockLengthHolder = 120050;
                            if ~data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(1,iJ)==1
                                % If start is missing first find end
                                [~, closest] = min(...
                                    abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.time(...
                                    data.et.(options.runType{iB}).blockEndIdx{iI,iZ}(iJ))));
                                data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                                clear closest
                                
                                % Subtract block length value from end time.
                                data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) = ...
                                    data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ)-...
                                    blockLengthHolder;
                                % B/c there may not be a time point at that
                                % exact value, find the nearest time point.
                                [~, closest] = min(...
                                    abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ)));
                                data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                                clear closest
                                
                            elseif ~data.et.(options.runType{iB}).blockEventsFound{iI,iZ}(2,iJ)==1
                                % If end is missing find start first
                                [~, closest] = min(...
                                    abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Messages.time(...
                                    data.et.(options.runType{iB}).blockStartIdx{iI,iZ}(iJ))));
                                data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                                clear closest
                                
                                % Add block length value to start time.
                                data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ) = ...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ)+...
                                    blockLengthHolder;
                                % B/c there may not be a time point at that
                                % exact value, find the nearest time point.
                                [~, closest] = min(...
                                    abs(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time-...
                                    data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ)));
                                data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(closest);
                                clear closest
                            end
                        end
                    end
                    
                    %% Find eyeblink events
                    leftCounter = 0;
                    rightCounter = 0;
                    for iE = 1:length(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.eye)
                        % Ignore events that are before or after start/end of block
                        if data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.start(iE) > ...
                                data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) && ...
                                data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.start(iE) <...
                                data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ)
                            if strcmp(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.eye{iE},'LEFT')
                                leftCounter = leftCounter + 1;
                                % Subtract the eyeblink value from block start time, to
                                % get relative eyeblink time.
                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(leftCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.start(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(leftCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.end(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                            elseif strcmp(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.eye{iE},'RIGHT')
                                rightCounter = rightCounter + 1;
                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(rightCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.start(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(rightCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Eblink.end(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                            end
                        end
                    end
                    % If no blinks during this block make variable empty
                    if leftCounter == 0
                        data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left = [];
                        data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left = [];
                    end
                    if rightCounter == 0
                        data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right = [];
                        data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right = [];
                    end
                    clear leftCounter rightCounter
                    
                    % Find blinks that occurred in both eyes 
                    % Loop through all blinks in Right eye and if there is
                    % a corresponding blink in the Left eye record as blink
                    % Look for blink times w/in +/-25ms 
                    % Find and loop through length of the smaller number of
                    % blink events, since the one w/ extra events will have
                    % the leftover events. 
                    [~, smallerIdx] = min([length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left),...
                        length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)]);
                    if smallerIdx == 1  % Left eye has less events
                        smallerArray = data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left;
                        largerArray = data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right;
                    elseif smallerIdx == 2   % Right eye has less events
                        smallerArray = data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right;
                        largerArray = data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left;
                    end
                    counter = 0;
                    for iE=1:length(smallerArray)
                        diffHolder = smallerArray(iE) - largerArray;
                        diffHolderIdx = abs(diffHolder)<50;
                        % if there is a corresponding time point
                        if sum(diffHolderIdx)>=1
                            counter=counter+1;
                            % If there are more than 1 difs <50ms find the
                            % smallest diff to pair
                            minDiffHolder = find(min(diffHolder(diffHolderIdx))==diffHolder);
                            data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.both(1,counter)=...
                                smallerArray(iE);
                            data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.both(2,counter)=...
                                largerArray(minDiffHolder);
                        end
                        clear diffHolder diffHolderIdx minDiffHolder 
                    end
                    % If no eye events that match between left and right
                    % for this block set array to empty
                    if counter==0
                        data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.both = [];
                    end
                    clear smallerIdx smallerArray largerArray counter
                    
                    % Count the total number of blinks
                    if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ})
                        if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                            data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).left = ...
                                numel(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left);
                        else   % If no eye this block make NaN (no eye recorded)
                            data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).left = ...
                                NaN;
                        end
                        if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                            data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).right = ...
                                numel(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right);
                        else   % If no eye this block make NaN (no eye recorded)
                            data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).right = ...
                                NaN;
                        end
                        if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.both)
                            data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).both = ...
                                size(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.both,2);
                        else   % If no eye this block make NaN (no eye recorded)
                            data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).both = ...
                                NaN;
                        end
                    else   % If no B/Z day for this participant make total NaN
                        data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).left = ...
                            NaN;
                        data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).right = ...
                            NaN;
                        data.et.(options.runType{iB}).blinkNumber(iI,iZ,iJ).both = ...
                            NaN;
                    end
                    
                    
                    %% Find eyemovement events
                    leftCounter = 0;
                    rightCounter = 0;
                    for iE = 1:length(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.eye)
                        % Ignore events that are before or after start/end of block
                        if data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.start(iE) > ...
                                data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ) && ...
                                data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.start(iE) <...
                                data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ)
                            if strcmp(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.eye{iE},'LEFT')
                                leftCounter = leftCounter + 1;
                                data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(leftCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.start(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                                data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.left(leftCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.end(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                            elseif strcmp(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.eye{iE},'RIGHT')
                                rightCounter = rightCounter + 1;
                                data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(rightCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.start(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                                data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.right(rightCounter) =...
                                    data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Events.eventCode.Esacc.end(iE)-...
                                    data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ);
                            end
                        end
                    end
                    % If no movements during this block make variable empty
                    if leftCounter == 0
                        data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left = [];
                        data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.left = [];
                    end
                    if rightCounter == 0
                        data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right = [];
                        data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.right = [];
                    end
                    clear leftCounter rightCounter
                    
                    % Find saccs that occurred in both eyes 
                    % Loop through all saccs in Right eye and if there is
                    % a corresponding sacc in the Left eye record as sacc
                    % Look for sacc times w/in +/-25ms 
                    % Find and loop through length of the smaller number of
                    % sacc events, since the one w/ extra events will have
                    % the leftover events. 
                    [~, smallerIdx] = min([length(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left),...
                        length(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right)]);
                    if smallerIdx == 1  % Left eye has less events
                        smallerArray = data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left;
                        largerArray = data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right;
                    elseif smallerIdx == 2   % Right eye has less events
                        smallerArray = data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right;
                        largerArray = data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left;
                    end
                    counter = 0;
                    for iE=1:length(smallerArray)
                        diffHolder = smallerArray(iE) - largerArray;
                        diffHolderIdx = abs(diffHolder)<50;
                        % if there is a corresponding time point
                        if sum(diffHolderIdx)>=1
                            counter=counter+1;
                            % If there are more than 1 difs <50ms find the
                            % smallest diff to pair
                            minDiffHolder = find(min(diffHolder(diffHolderIdx))==diffHolder);
                            data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.both(1,counter)=...
                                smallerArray(iE);
                            data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.both(2,counter)=...
                                largerArray(minDiffHolder);
                        end
                        clear diffHolder diffHolderIdx minDiffHolder 
                    end
                    % If no eye events that match between left and right
                    % for this block set array to empty
                    if counter==0
                        data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.both = [];
                    end
                    clear smallerIdx smallerArray largerArray counter
                    
                    % Count the total number of movements
                    if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ})
                        if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left)
                            data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).left = ...
                                numel(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left);
                        else   % If no eye this block make NaN (no eye recorded)
                            data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).left = ...
                                NaN;
                        end
                        if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right)
                            data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).right = ...
                                numel(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right);
                        else   % If no eye this block make NaN (no eye recorded)
                            data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).right = ...
                                NaN;
                        end
                        if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.both)
                            data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).both = ...
                                size(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.both,2);
                        else   % If no eye this block make NaN (no eye recorded)
                            data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).both = ...
                                NaN;
                        end
                    else   % If no B/Z day for this participant make total NaN
                        data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).left = ...
                            NaN;
                        data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).right = ...
                            NaN;
                        data.et.(options.runType{iB}).saccNumber(iI,iZ,iJ).both = ...
                            NaN;
                    end
                    
                    %% Pull out the relevant data
                    % Find the index for first/last data point in the samples array
                    data.et.(options.runType{iB}).fixDataStartIdx{iI,iZ}(iJ) =...
                        find(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time==...
                        data.et.(options.runType{iB}).blockStartTime{iI,iZ}(iJ));
                    data.et.(options.runType{iB}).fixDataEndIdx{iI,iZ}(iJ) =...
                        find(data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time==...
                        data.et.(options.runType{iB}).blockEndTime{iI,iZ}(iJ));
                    
                    % Grab the first and last time point
                    data.et.(options.runType{iB}).fixDataStartTime{iI,iZ}(iJ) =...
                        data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(...
                        data.et.(options.runType{iB}).fixDataStartIdx{iI,iZ}(iJ),:);
                    data.et.(options.runType{iB}).fixDataEndTime{iI,iZ}(iJ) =...
                        data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.time(...
                        data.et.(options.runType{iB}).fixDataEndIdx{iI,iZ}(iJ),:);
                    data.et.(options.runType{iB}).fixDataTotalTime{iI,iZ}(iJ) =...
                        data.et.(options.runType{iB}).fixDataEndTime{iI,iZ}(iJ) -...
                        data.et.(options.runType{iB}).fixDataStartTime{iI,iZ}(iJ);
                    
                    % Grab the data
                    % There's 'posX' and 'posY' which I thought were position, but
                    % they seem to match 'gx' and 'gy' (gaze?).
                    data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ} = ...
                        data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.gx(...
                        data.et.(options.runType{iB}).fixDataStartIdx{iI,iZ}(iJ):...
                        data.et.(options.runType{iB}).fixDataEndIdx{iI,iZ}(iJ),:);
                    data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ} = ...
                        data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.gy(...
                        data.et.(options.runType{iB}).fixDataStartIdx{iI,iZ}(iJ):...
                        data.et.(options.runType{iB}).fixDataEndIdx{iI,iZ}(iJ),:);
                    % Grab pupil size data as well
                    data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ} = ...
                        data.et.(options.runType{iB}).etDataOrig{iI,iZ}.Samples.pupilSize(...
                        data.et.(options.runType{iB}).fixDataStartIdx{iI,iZ}(iJ):...
                        data.et.(options.runType{iB}).fixDataEndIdx{iI,iZ}(iJ),:);
                    
                    % Sometimes the pos data is > total screen  area...if so make
                    % it NaN
                    data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ}(...
                        data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ}>options.screenSize(3)) =...
                        NaN;
                    data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ}(...
                        data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ}<options.screenSize(1)) =...
                        NaN;
                    data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ}(...
                        data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ}>options.screenSize(4)) =...
                        NaN;
                    data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ}(...
                        data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ}<options.screenSize(2)) =...
                        NaN;
                    
                    % Get rid of the first 10ms (make nans) before and after a blink
                    blinkBuffer = 10;
                    eyeLabel = {'left','right'};
                    for iE=1:2   % For left and right eyes
                        if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},eyeLabel(iE))
                            if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE}))
                                for iK=1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE}))
                                    % If the blink starts or ends beyond +/- the blink
                                    % buffer, than only extend to beginning/end.
                                    if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2 <= blinkBuffer
                                        blinkBufferHolder =...
                                            (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2) - 1;
                                    else
                                        blinkBufferHolder = blinkBuffer;
                                    end
                                    data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ}(...
                                        (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2)-blinkBufferHolder:...
                                        (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2),:) = NaN;
                                    data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ}(...
                                        (data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2):...
                                        (data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2)+blinkBufferHolder,:) = NaN;
                                    data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ}(...
                                        (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2)-blinkBufferHolder:...
                                        (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2),:) = NaN;
                                    data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ}(...
                                        (data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2):...
                                        (data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2)+blinkBufferHolder,:) = NaN;
                                    % Do the same for pupil size
                                    data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ}(...
                                        (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2)-blinkBufferHolder:...
                                        (data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2),:) = NaN;
                                    data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ}(...
                                        (data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2):...
                                        (data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.(eyeLabel{iE})(iK)/2)+blinkBufferHolder,:) = NaN;
                                end
                            end
                        end
                    end
                    
                    % Take the average between the left and right eyes
                    data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ} = nanmean(data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ},2);
                    data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ} = nanmean(data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ},2);
                    data.et.(options.runType{iB}).pupilSize_ave{iI,iZ}{iJ} = nanmean(data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ},2);
                    
                    % If the position is greater than or less than the
                    % screen size make it NaN.
                    data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}>options.screenSize(3)) =...
                        NaN;
                    data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}<options.screenSize(1)) =...
                        NaN;
                    data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}>options.screenSize(4)) =...
                        NaN;
                    data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}<options.screenSize(2)) =...
                        NaN;         
                    
                    % Take the difference between fixation center and eye position
                    data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_posX{iI,iZ}{iJ} - ...
                        (options.screenSize(3)/2);
                    data.et.(options.runType{iB}).fixData_diffPosY{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_posY{iI,iZ}{iJ} - ...
                        (options.screenSize(4)/2);
                    data.et.(options.runType{iB}).fixData_aveDiffPosX{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ} - ...
                        (options.screenSize(3)/2);
                    data.et.(options.runType{iB}).fixData_aveDiffPosY{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ} - ...
                        (options.screenSize(4)/2);
                    
                    % Calculate the average distance in position from fixation in x and y
                    data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ}(iJ) = ...
                        nanmean(data.et.(options.runType{iB}).fixData_aveDiffPosX{iI,iZ}{iJ});
                    data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ}(iJ) = ...
                        nanmean(data.et.(options.runType{iB}).fixData_aveDiffPosY{iI,iZ}{iJ});
                    
                    % Calculate the Euclidian distance from  fixation of
                    % the average position.
                    data.et.(options.runType{iB}).eucDistFromFix{iI,iZ}{iJ} = sqrt(...
                        (data.et.(options.runType{iB}).fixData_aveDiffPosX{iI,iZ}{iJ} .^ 2) + ...
                        (data.et.(options.runType{iB}).fixData_aveDiffPosY{iI,iZ}{iJ} .^ 2));
                    
                    % Calculate mean/std distance uncorrected
                    data.et.(options.runType{iB}).stats.aveEucDisFromFix{iI,iZ}(iJ) =...
                        nanmean(data.et.(options.runType{iB}).eucDistFromFix{iI,iZ}{iJ});
                    data.et.(options.runType{iB}).stats.stdEucDisFromFix{iI,iZ}(iJ) =...
                        nanstd(data.et.(options.runType{iB}).eucDistFromFix{iI,iZ}{iJ});
                    
                    %% Make corrections to position
                    % Correct the position data based on calc'd correction value
                    data.et.(options.runType{iB}).corFixData_avePosX{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ} -...
                        data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ}(iJ);
                    data.et.(options.runType{iB}).corFixData_avePosY{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ} -...
                        data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ}(iJ);
                    
                    % Find the new distance from fixation for the corrected data
                    data.et.(options.runType{iB}).corFixData_aveDiffPosX{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).corFixData_avePosX{iI,iZ}{iJ} - ...
                        (options.screenSize(3)/2);
                    data.et.(options.runType{iB}).corFixData_aveDiffPosY{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).corFixData_avePosY{iI,iZ}{iJ} - ...
                        (options.screenSize(4)/2);
                    
                    % Calculate corrected Euclidian distance
                    data.et.(options.runType{iB}).corrEucDistFromFix{iI,iZ}{iJ} = sqrt(...
                        (data.et.(options.runType{iB}).corFixData_aveDiffPosX{iI,iZ}{iJ} .^ 2) + ...
                        (data.et.(options.runType{iB}).corFixData_aveDiffPosY{iI,iZ}{iJ} .^ 2));
                    
                    % Calculate mean/std distance
                    data.et.(options.runType{iB}).stats.corrAveEucDisFromFix{iI,iZ}(iJ) =...
                        nanmean(data.et.(options.runType{iB}).corrEucDistFromFix{iI,iZ}{iJ});
                    data.et.(options.runType{iB}).stats.corrStdEucDisFromFix{iI,iZ}(iJ) =...
                        nanstd(data.et.(options.runType{iB}).corrEucDistFromFix{iI,iZ}{iJ});
                    
                    
                    %% Make the heat map for this block
                    % Make for uncorrected data
                    % First make a grid (array) that is the size of the screen
                    data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{iJ} =...
                        zeros([options.screenSize(3) options.screenSize(4)]);
                    
                    % Find the heat map array indices (fixation position) that need to be counted.
                    clear finderArray1 finderArray2
                    finderArray1(1,:) = round(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}(:))-1;
                    finderArray1(2,:) = round(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}(:));
                    finderArray1(3,:) = round(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ}(:))+1;
                    finderArray2(1,:) = round(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}(:))-1;
                    finderArray2(2,:) = round(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}(:));
                    finderArray2(3,:) = round(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ}(:))+1;
                    % B/c we are expanding area around fixation center,
                    % If there are values > max(screen) or < min(screen), erase
                    % them. (Can't count up for indices that don't exist).
                    finderArray1...
                        ((finderArray1 >= (options.screenSize(3)+1) |...
                        finderArray1 <= (options.screenSize(1)-1))) = NaN;
                    finderArray2...
                        ((finderArray2 >= (options.screenSize(4)+1) |...
                        finderArray2 <= (options.screenSize(2)-1))) = NaN;
                    
                    % For each fixation, determine it's location w/in the array +/-1 and
                    % count up 1.
                    for iK=1:3
                        for iL=1:length(finderArray1)
                            if ~isnan(finderArray1(iK,iL)) & ~isnan(finderArray2(iK,iL))
                                data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{iJ}(...
                                    finderArray1(iK,iL),...
                                    finderArray2(iK,iL)) = ...
                                    data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{iJ}(...
                                    finderArray1(iK,iL),...
                                    finderArray2(iK,iL)) + 1;
                            end
                        end
                    end
                    
                    % Count the total number of non nan fixations
                    data.et.(options.runType{iB}).nonNanCounter{iI,iZ}(iJ) =...
                        sum(~isnan(round(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ})));
                    
                    % Now divide by the total number of non-nan fixations
                    data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{iJ} ./...
                        data.et.(options.runType{iB}).nonNanCounter{iI,iZ}(iJ);
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    % Make for corrected data
                    % First make a grid (array) that is the size of the screen
                    data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{iJ} =...
                        zeros([options.screenSize(3) options.screenSize(4)]);
                    
                    % Find the heat map array indices (fixation position) that need to be counted.
                    clear finderArray1 finderArray2
                    finderArray1(1,:) = round(data.et.(options.runType{iB}).corFixData_avePosX{iI,iZ}{iJ}(:))-1;
                    finderArray1(2,:) = round(data.et.(options.runType{iB}).corFixData_avePosX{iI,iZ}{iJ}(:));
                    finderArray1(3,:) = round(data.et.(options.runType{iB}).corFixData_avePosX{iI,iZ}{iJ}(:))+1;
                    finderArray2(1,:) = round(data.et.(options.runType{iB}).corFixData_avePosY{iI,iZ}{iJ}(:))-1;
                    finderArray2(2,:) = round(data.et.(options.runType{iB}).corFixData_avePosY{iI,iZ}{iJ}(:));
                    finderArray2(3,:) = round(data.et.(options.runType{iB}).corFixData_avePosY{iI,iZ}{iJ}(:))+1;
                    % B/c we are expanding area around fixation center,
                    % If there are values > max(screen) or < min(screen), erase
                    % them. (Can't count up for indices that don't exist).
                    finderArray1...
                        ((finderArray1 >= (options.screenSize(3)+1) |...
                        finderArray1 <= (options.screenSize(1)-1))) = NaN;
                    finderArray2...
                        ((finderArray2 >= (options.screenSize(4)+1) |...
                        finderArray2 <= (options.screenSize(2)-1))) = NaN;
                    
                    % For each fixation, determine it's location w/in the array +/-1 and
                    % count up 1.
                    counter = 0;
                    for iK=1:3
                        for iL=1:length(finderArray1)
                            if ~isnan(finderArray1(iK,iL)) & ~isnan(finderArray2(iK,iL))
                                counter=counter+1;
                                data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{iJ}(...
                                    finderArray1(iK,iL),...
                                    finderArray2(iK,iL)) = ...
                                    data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{iJ}(...
                                    finderArray1(iK,iL),...
                                    finderArray2(iK,iL)) + 1;
                            end
                        end
                    end
                    
                    % Count the total number of non nan fixations
                    data.et.(options.runType{iB}).nonNanCounterCorr{iI,iZ}(iJ) =...
                        sum(~isnan(round(data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ})));
                    
                    % Now divide by the total number of non-nan fixations
                    data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{iJ} =...
                        data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{iJ} ./...
                        data.et.(options.runType{iB}).nonNanCounterCorr{iI,iZ}(iJ);
                end
            end
        end
    end
end

%% Create demographics table for participants used here

% Create index to only include parts/dates that are not nans
demogIdx = ~isnan(data.et.B.date);
subjHolder = repmat(data.et.B.subjID, [1 2]);
groupHolder = repmat([ones([length(data.et.B.grpIdx{1}) 1]);...
    ones([length(data.et.B.grpIdx{2}) 1])+1;...
    ones([length(data.et.B.grpIdx{3}) 1])+2],[1 2]);
% Convert date to proper format
for iI = 1:size(data.et.B.date,1)
    for iJ = 1:size(data.et.B.date,2)
        if ~isnan(data.et.B.date(iI,iJ))
            dateHolder(iI,iJ) = datenum(num2str(data.et.B.date(iI,iJ)),'yyyymmdd');
        else
           dateHolder(iI,iJ) = NaN; 
        end
    end
end

counter = 1;
for iI = 1:length(demogIdx)
    demogOptions.subj_number(counter:(counter-1)+sum(demogIdx(iI,:))) = ...
        subjHolder(iI,demogIdx(iI,:));
    demogOptions.partGroup(counter:(counter-1)+sum(demogIdx(iI,:))) = ...
        groupHolder(iI,demogIdx(iI,:));
    demogOptions.date_number(counter:(counter-1)+sum(demogIdx(iI,:))) = ...
        dateHolder(iI,demogIdx(iI,:));
    counter = counter+sum(demogIdx(iI,:));
end
demogOptions.subj_group_def = 1;

% Only look at unique subjects
[~, uniqueIdx, repeatSubjs] = unique(demogOptions.subj_number');
demogOptions.subj_number = demogOptions.subj_number(uniqueIdx);
demogOptions.partGroup = demogOptions.partGroup(uniqueIdx);
demogOptions.date_number = demogOptions.date_number(uniqueIdx);
clear uniqueIdx

% Call MPS phcp_demographics function
cd /home/shaw-raid1/matlab_tools/mpsCode/
options.demoData = phcp_demographics(demogOptions);
cd(options.curDur);
clear demogOptions demogIdx groupHolder subjHolder dateHolder

% Make table
rowNames = {'Age','Gender','Education','Estimated_IQ','Visual_Acuity'};
columnNames = {sprintf('%s%d%s','g1_',options.demoData.group_n_unique(1)),...
    sprintf('%s%d%s','g2_',options.demoData.group_n_unique(2)),...
    sprintf('%s%d%s','g3_',options.demoData.group_n_unique(3)),...
    'stats'};

g1 = {sprintf('%.3f%s%.3f%s',options.demoData.Age.('g1').mean,'(',options.demoData.Age.('g1').sd,')'),...
    sprintf('%2.0f%s',options.demoData.Gender.('g1').pct_female*100,'%'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Education.('g1').mean,'(',options.demoData.Education.('g1').sd,')'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Estimated_IQ.('g1').mean,'(',options.demoData.Estimated_IQ.('g1').sd,')'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Visual_Acuity.('g1').mean,'(',options.demoData.Visual_Acuity.('g1').sd,')')}';
g2 = {sprintf('%.3f%s%.3f%s',options.demoData.Age.('g2').mean,'(',options.demoData.Age.('g2').sd,')'),...
    sprintf('%2.0f%s',options.demoData.Gender.('g2').pct_female*100,'%'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Education.('g2').mean,'(',options.demoData.Education.('g2').sd,')'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Estimated_IQ.('g2').mean,'(',options.demoData.Estimated_IQ.('g2').sd,')'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Visual_Acuity.('g2').mean,'(',options.demoData.Visual_Acuity.('g2').sd,')')}';
g3 = {sprintf('%.3f%s%.3f%s',options.demoData.Age.('g3').mean,'(',options.demoData.Age.('g3').sd,')'),...
    sprintf('%2.0f%s',options.demoData.Gender.('g3').pct_female*100,'%'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Education.('g3').mean,'(',options.demoData.Education.('g3').sd,')'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Estimated_IQ.('g3').mean,'(',options.demoData.Estimated_IQ.('g3').sd,')'),...
    sprintf('%.3f%s%.3f%s',options.demoData.Visual_Acuity.('g3').mean,'(',options.demoData.Visual_Acuity.('g3').sd,')')}';
stats = {sprintf('%s%d%s%d%s%1.3f%s%1.3f','F(',options.demoData.Age.('stats').all3{3,3},',',...
    options.demoData.Age.('stats').all3{2,3},')=',options.demoData.Age.('stats').all3{3,6},...
    '; p=',options.demoData.Age.('stats').all3{3,7}),...
    
    sprintf('%s%d%s%1.3f%s%1.3f','X2(',options.demoData.Gender.('stats').all3.df,')=',...
    options.demoData.Gender.('stats').all3.Chi_squared,'; p=',...
    options.demoData.Gender.('stats').all3.p_value),...
    
    sprintf('%s%d%s%d%s%1.3f%s%1.3f','F(',options.demoData.Education.('stats').all3{3,3},',',...
    options.demoData.Education.('stats').all3{2,3},')=',options.demoData.Education.('stats').all3{3,6},...
    '; p=',options.demoData.Education.('stats').all3{3,7}),...
    
    sprintf('%s%d%s%d%s%1.3f%s%1.3f','F(',options.demoData.Estimated_IQ.('stats').all3{3,3},',',...
    options.demoData.Estimated_IQ.('stats').all3{2,3},')=',options.demoData.Estimated_IQ.('stats').all3{3,6},...
    '; p=',options.demoData.Estimated_IQ.('stats').all3{3,7}),...
    
    sprintf('%s%d%s%d%s%1.3f%s%1.3f','F(',options.demoData.Visual_Acuity.('stats').all3{3,3},',',...
    options.demoData.Visual_Acuity.('stats').all3{2,3},')=',options.demoData.Visual_Acuity.('stats').all3{3,6},...
    '; p=',options.demoData.Visual_Acuity.('stats').all3{3,7})};

% Convert to table format
options.demoData.Table = table(g1,g2,g3,stats,'RowNames',rowNames);
options.demoData.Table.Properties.VariableNames = columnNames;
clear columnNames rowNames

%% Exclude participants based on poor fixation
% If participant has a corrected euc average > 1.5 degs and std > 1.5 then cut them.
% KWK 20230720
% for iB=1:2   % A/B runs
for iI=1:size(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix,1)   % Subj
    for iZ=1:size(data.et.B.stats.corrAveEucDisFromFix,2)   % B/Z days
        if isempty(data.et.B.stats.corrAveEucDisFromFix{iI,iZ})
            partExclusionHolder(iI,iZ) = 0;
        elseif nanmean(data.et.B.stats.corrAveEucDisFromFix{iI,iZ} ./ options.PPD)>1.5 && ...
                nanmean(data.et.B.stats.corrAveEucDisFromFix{iI,iZ} ./ options.PPD)>1.5
            partExclusionHolder(iI,iZ) = 1;
        else
            partExclusionHolder(iI,iZ) = 0;
        end
    end
end
% % end
data.et.B.partExclusion = partExclusionHolder(:,1) | partExclusionHolder(:,2);
data.et.B.partExclusion = boolean(data.et.B.partExclusion);
clear partExclusionHolder

% Record the excluded subjIDs
data.et.B.partExclusionSubjID = data.et.B.subjID(data.et.B.partExclusion);

% Clear excluded subjects
runLabel = {'A' 'B'};
for iR=1:2   % A/B runs
    fieldNameHolder = fieldnames(data.et.(runLabel{iR}));
    for iF=1:size(fieldNameHolder)
        if strcmp(fieldNameHolder(iF),'stats')
            fieldNameHolderStats = fieldnames(data.et.(runLabel{iR}).stats);
            for iFF=1:size(fieldNameHolderStats)
                data.et.(runLabel{iR}).stats.(fieldNameHolderStats{iFF})(data.et.B.partExclusion(:,1),:) = [];
            end
        elseif strcmp(fieldNameHolder(iF),'grpIdx')   % Redo the grping index
        elseif strcmp(fieldNameHolder(iF),'partExclusionSubjID')   % Don't exclude from the exclusion subj lists
        elseif strcmp(fieldNameHolder(iF),'partExclusion')
        else
            data.et.(runLabel{iR}).(fieldNameHolder{iF})(data.et.B.partExclusion(:,1),:) = [];
        end
    end
    clear fieldNameHolder fieldNameHolderStats
end

fieldNameHolder = fieldnames(data.behav);
for iF=1:size(fieldNameHolder)
    if strcmp(fieldNameHolder(iF),'A') || strcmp(fieldNameHolder(iF),'B')
        fieldNameHolderRun = fieldnames(data.behav.(runLabel{iR}));
        for iFF=1:size(fieldNameHolderRun)
            data.behav.(fieldNameHolder{iF}).(fieldNameHolderRun{iFF})(data.et.B.partExclusion(:,1),:) = [];
        end
        clear fieldNameHolderRun
    elseif strcmp(fieldNameHolder(iF),'controlSwitchTimesOrig') || ...
            strcmp(fieldNameHolder(iF),'controlSwitchNames') || ...
            strcmp(fieldNameHolder(iF),'controlSwitchTimes') || ...
            strcmp(fieldNameHolder(iF),'controlSwitchInd')
    elseif strcmp(fieldNameHolder(iF),'grpIdx')   % Redo the grping index
    else
        data.behav.(fieldNameHolder{iF})(data.et.B.partExclusion(:,1),:) = [];
    end
end
clear fieldNameHolder

% Redo the grouping index
data.behav.grpIdx{1} = find(data.behav.subjID < 2000000);
data.behav.grpIdx{2} = find(data.behav.subjID > 2000000 & data.behav.subjID < 6000000);
data.behav.grpIdx{3} = find(data.behav.subjID > 6000000);

data.et.A.grpIdx{1} = find(data.et.A.subjID < 2000000);
data.et.A.grpIdx{2} = find(data.et.A.subjID > 2000000 & data.et.A.subjID < 6000000);
data.et.A.grpIdx{3} = find(data.et.A.subjID > 6000000);
data.et.B.grpIdx{1} = find(data.et.B.subjID < 2000000);
data.et.B.grpIdx{2} = find(data.et.B.subjID > 2000000 & data.et.B.subjID < 6000000);
data.et.B.grpIdx{3} = find(data.et.B.subjID > 6000000);

% Find the number of subjects with only one eye
for iI=1:size(data.et.B.fixData_posX,1)   % Subject
    for iZ = 1:size(data.et.B.fixData_posX,2)   % B/Z
        % Only need to look at 1 block to see if there's data for both eyes
        if isempty(data.et.B.fixData_posX{iI,iZ})
            data.et.B.eyesCollected(iI,iZ) = NaN;
        elseif isnan(nanmean(data.et.B.fixData_posX{iI,iZ}{1}(:,1))) && ...
                ~isnan(nanmean(data.et.B.fixData_posX{iI,iZ}{1}(:,2)))   % right eye nan
            data.et.B.eyesCollected(iI,iZ) = 1;
        elseif isnan(nanmean(data.et.B.fixData_posX{iI,iZ}{1}(:,2))) && ...
                ~isnan(nanmean(data.et.B.fixData_posX{iI,iZ}{1}(:,1)))   % left eye nan
            data.et.B.eyesCollected(iI,iZ) = 2;
        elseif ~isnan(nanmean(data.et.B.fixData_posX{iI,iZ}{1}(:,1))) && ...
                ~isnan(nanmean(data.et.B.fixData_posX{iI,iZ}{1}(:,2)))   % neither eye nan
            data.et.B.eyesCollected(iI,iZ) = 3;
        end
    end
end
holder = [data.et.B.eyesCollected==1 | data.et.B.eyesCollected==2];
holder = holder(:,1) | holder(:,2);
data.et.B.eyesCollected_subjID = data.et.B.subjID(holder);
clear holder


%% Run contingency table analysis to look for differences between number of excluded participants
% Make new array to exclude repeats in counts
n_cat = 2;
samples = [sum(data.et.B.subjID<2000000)...
    sum(data.et.B.partExclusionSubjID<2000000);...
    sum(data.et.B.subjID>=2000000 & data.et.B.subjID<6000000)...
    sum(data.et.B.partExclusionSubjID>=2000000 & data.et.B.partExclusionSubjID<6000000);...
    sum(data.et.B.subjID>6000000)...
    sum(data.et.B.partExclusionSubjID>6000000)];
yates = 0;
[statsHolder,dataHolder] = mpsContingencyTable(n_cat, samples, yates);
data.et.B.partExclusionContingency.data = dataHolder;
data.et.B.partExclusionContingency.stats = statsHolder;
clear statsHolder dataHolder

%% Save data structure
options.dateRun = datestr(now,'YYYYMMDD');
% data is super big... maybe don't save on the server?
% save(['SFM_ET_Analysis_Output_' options.dateRun],'data','options','-v7.3')

%% Plot the fixation data for A over time for each block for each subject
% Graph params
if options.displayETFigs
    % Group colors
    col_list{1} = [0 1 0];
    col_list{2} = [0 0 1];
    col_list{3} = [1 0 0];
    
    % Group labels
    for iB=1:2   % A/B runs
        options.(['group_labels_', options.runType{iB}]) = {...
            ['C, n=' num2str(length(data.et.(options.runType{iB}).grpIdx{1}))],...
            ['R, n=' num2str(length(data.et.(options.runType{iB}).grpIdx{2}))],...
            ['P, n=' num2str(length(data.et.(options.runType{iB}).grpIdx{3}))]};
    end
    
    % Heat map resolution
    % Set 'zoom' resolution for the heat maps, so its not the entire 1920x180 screen
    heatMapZoom = 6; % in degrees
    % Only plot heatmap values from center +/- zoom res
    heatMapIdx = round([(options.screenSize(3)/2) - (heatMapZoom*options.PPD),...
        (options.screenSize(4)/2) - (heatMapZoom*options.PPD),...
        (options.screenSize(3)/2) + (heatMapZoom*options.PPD),...
        (options.screenSize(4)/2) + (heatMapZoom*options.PPD)]);
    heatMapSize = [0 0 ...
        heatMapIdx(3)-heatMapIdx(1) heatMapIdx(4)-heatMapIdx(2)];
    heatMapXTicks = round(linspace(1,heatMapSize(3),heatMapZoom*2+1));
    heatMapXLabels = -heatMapZoom:1:heatMapZoom;
end
if options.displayETFigs
    plotNumRows = 5;
    
    % Make x-axis markers for time plots
    % (Used to draw x-axis for each of the 5 plots in
    % the figure)
    clear xAxisMarkers xAxisMarkers
    % Maximum range above/below each x axis
    maxYVal = 200/options.PPD;   % Convert to degrees
    xAxisMarkers = flip(maxYVal:maxYVal*2:maxYVal*10);
    
    for iI=1
        % for iI=1:size(data.et.(options.runType{iB}).etDataOrig,1)
        for iB=1:2   % A/B runs
            for iZ=1:size(data.et.(options.runType{iB}).etDataOrig,2)   % B/Z days
                
                % Make size values for individual plots
                clear timePlotSize scatterPlotSize timePupPlotSize
                if iB==1
                    timePlotSize(1,:) = [.33 .01 .33 .98];
                    timePupPlotSize = [.05 .01 .94 .94];
                    scatterPlotSize(1,:) = [0.03 .01 .28 .98];
                elseif iB==2
                    for iJ=1:length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ})
                        timePlotSize(iJ,:) = [.33 1.01-(.2*iJ) .33 .19];
                        timePupPlotSize(iJ,:) = [.05 .95-((.94*.2)*iJ) .94 (.94*.2)];
                        scatterPlotSize(iJ,:) = [.03 1.01-(.2*iJ) .28 .19];
                    end
                end
                
                % if there is a date for the fix data for this part for this run,
                % then plot the data
                if ~isnan(data.et.(options.runType{iB}).date(iI,iZ))

                    fig1=figure();
                    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
                    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
                    figSize.figSize = [0 0 ...
                        figSize.baseSize(3)...
                        figSize.baseSize(4)];   % Size/postion of fig
                    set(gcf,'color','w')
                    set(gcf,'Position', figSize.figSize)
                    
                    for iJ=1:length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ})   % Blocks
                        %% Plot position data in screen coords
                        % subplot(5,3,1+((iJ-1)*3))
                        
                        scatterPlot(iJ) = axes('Parent',fig1);
                        scatter(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ},...
                            data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ});
                        hold on
                        plot([options.screenSize(3)/2 options.screenSize(3)/2],...
                            [options.screenSize(4)/2 options.screenSize(4)/2],'.r');
                        th = 0:pi/50:2*pi;
                        xDegCirc = (options.PPD/2) * cos(th) + (options.screenSize(3)/2);
                        yDegCirc = (options.PPD/2) * sin(th) + (options.screenSize(4)/2);
                        plot(xDegCirc, yDegCirc,'k');
                        text(options.screenSize(3)*.05,options.screenSize(4)*.4,...
                            sprintf('%s\n%.2f','X Correction Used:',...
                            data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ}(iJ)/...
                            options.PPD,'°'));
                        text(options.screenSize(3)*.05,options.screenSize(4)*.15,...
                            sprintf('%s\n%.2f%s','Y Correction Used:',...
                            data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ}(iJ)/...
                            options.PPD,'°'));
                        xlim([options.screenSize(1) options.screenSize(3)])
                        ylim([options.screenSize(2) options.screenSize(4)])
                        hold off
                        
                        box off
%                         title(['Fixation points (block: ', ...
%                             num2str(iJ),')'],...
%                             'fontsize',10)
                        
                        set(scatterPlot(iJ),'xtick',[],...
                            'ytick',[],...
                            'Position',scatterPlotSize(iJ,:))
                        ylabel(['Block: ', num2str(iJ)])
                        xlabel([])
                    end
                        
                    %% Plot the heat map
                    % subplot(5,3,[3,6])
                    heatPlot(1) = axes('Parent',fig1);
                    % Average data across blocks
                    holderArray = cellfun(@(x) x(heatMapIdx(1):heatMapIdx(3),...
                        heatMapIdx(2):heatMapIdx(4)),...
                        data.et.(options.runType{iB}).fixData_heatMap{iI,iZ},'UniformOutput',false);
                    imagesc(nanmean(cat(3,holderArray{:}),3)');
                    clear holderArray
                    hold on
                    plot([heatMapSize(3)/2 heatMapSize(3)/2],...
                        [heatMapSize(4)/2 heatMapSize(4)/2],'.r');
                    set(gca,'YDir','normal')
                    th = 0:pi/50:2*pi;
                    xDegCirc = (options.PPD/2) * cos(th) + (heatMapSize(3)/2);
                    yDegCirc = (options.PPD/2) * sin(th) + (heatMapSize(4)/2);
                    plot(xDegCirc, yDegCirc,'k');
                    text(heatMapSize(3)*.1,heatMapSize(4)*.25,...
                            sprintf('%s\n%.2f%s','Ave Euc. Distance From Fix:',...
                            nanmean(data.et.B.stats.aveEucDisFromFix{iI,iZ})/...
                            options.PPD,'°'));
                    text(heatMapSize(3)*.1,heatMapSize(4)*.1,...
                        sprintf('%s\n%d','Num Valid Fixations:',...
                        data.et.(options.runType{iB}).nonNanCounter{iI,iZ}(iJ)))
                    title(['Proportion of eye gaze position' ...
                        sprintf('\n%s%d%s%d%s',' (P',data.et.(options.runType{iB}).subjID(iI),', ',iZ,' run)')],...
                        'fontsize',10)
                    set(heatPlot(1),'xtick',heatMapXTicks,...
                        'xticklabels',heatMapXLabels,...
                        'ytick',[],...
                        'OuterPosition',[.66 .5 .33 .5])
                    
                    % subplot(5,3,[12,15])
                    heatPlot(2) = axes('Parent',fig1);
                    % Average data across blocks
                    holderArray = cellfun(@(x) x(heatMapIdx(1):heatMapIdx(3),...
                        heatMapIdx(2):heatMapIdx(4)),...
                        data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ},'UniformOutput',false);
                    imagesc(nanmean(cat(3,holderArray{:}),3)');
                    clear holderArray
                    hold on
                    plot([heatMapSize(3)/2 heatMapSize(3)/2],...
                        [heatMapSize(4)/2 heatMapSize(4)/2],'.r');
                    set(gca,'YDir','normal')
                    th = 0:pi/50:2*pi;
                    xDegCirc = (options.PPD/2) * cos(th) + (heatMapSize(3)/2);
                    yDegCirc = (options.PPD/2) * sin(th) + (heatMapSize(4)/2);
                    plot(xDegCirc, yDegCirc,'k');
                    text(heatMapSize(3)*.1,heatMapSize(4)*.25,...
                            sprintf('%s\n%.2f','Ave Euc. Distance From Fix (cor.):',...
                            nanmean(data.et.B.stats.corrAveEucDisFromFix{iI,iZ})/...
                            options.PPD,'°'));
                    text(heatMapSize(3)*.1,heatMapSize(4)*.1,...
                        sprintf('%s\n%d','Num Valid Fixations:',...
                        data.et.(options.runType{iB}).nonNanCounterCorr{iI,iZ}(iJ)))
                    set(gca,'XColor','k','YColor','k')
                    title(['Proportion of eye gaze position (corr.)' ...
                        sprintf('\n%s%d%s%d%s',' (P',data.et.(options.runType{iB}).subjID(iI),', ',iZ,' run)')],...
                        'fontsize',10)
                    set(heatPlot(2),'xtick',heatMapXTicks,...
                        'xticklabels',heatMapXLabels,...
                        'ytick',[],...
                        'OuterPosition',[.66 0 .33 .5])
                        
                        
                        
                    %% Plot relative eye position in lx, rx, ly, ry
                    % Find length and slpit evenly into 5 rows to plot
                    for iJ=1:length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ})   % Blocks
                        
                        % Make x-axis
                        % (Will be same x vals for all plots)
                        xAxisLims = linspace(1,length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ}),plotNumRows+1);
                        
                        % Find x-axis labels based on time
                        xAxisTime = linspace(1,data.et.(options.runType{iB}).fixDataTotalTime{iI,iZ}(iJ),plotNumRows+1);
                        
                        % Plot 5 sets of lines w/ full eye position for
                        % each block
                        timePlot(iJ) = axes('Parent',fig1);
                        for iK=1:plotNumRows
                            
                            % Center
                            h1=plot([xAxisLims(1) xAxisLims(2)],[xAxisMarkers(iK) xAxisMarkers(iK)],'k');
                            h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            hold on
                            % X Right eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                (data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),1) / options.PPD) + xAxisMarkers(iK),...
                                'Color',[0 .9 0]);
                            % Y Right eye
                            h2 = plot(xAxisLims(1):xAxisLims(2),...
                                (data.et.(options.runType{iB}).fixData_diffPosY{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),1) / options.PPD) + xAxisMarkers(iK),...
                                'Color',[0 .6 0]);
                            h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            % X Left eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                (data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),2) / options.PPD) + xAxisMarkers(iK),...
                                'Color',[0 0 .9])
                            % Y Left eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                (data.et.(options.runType{iB}).fixData_diffPosY{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),2) / options.PPD) + xAxisMarkers(iK),...
                                'Color',[0 0 .6])
                        
                            
                            %% Plot blink marker
                            % Left eye blink
                            % If there are left blinks for this trial
                            if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},'left')
                                % If there are blinks for this trial
                                if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                                    for iL = 1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                                        % If the blink time falls w/in the plot time limit plot it
                                        if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL) > xAxisTime(iK) &&...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL) < xAxisTime(iK+1)
                                            
                                            % Scale x values (eyeblink times) based on xAxisTime(iK) to ensure each
                                            % point is plotted relative to x axis '0'. Scale y values relative to 
                                            % individual x axis markers.
                                            patch([data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL)/2]-...
                                                xAxisTime(iK)/2,...
                                                [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
                                                [0 1 0],...
                                                'FaceAlpha',.5,...
                                                'EdgeColor','none')
                                        end
                                    end
                                end
                            end
                            % Right eye blink
                            % If there are rightblinks for this trial
                            if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},'right')
                                % If there are blinks for this trial
                                if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                                    for iL = 1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                                        % If the blink time falls w/in the plot time limit plot it
                                        if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL) > xAxisTime(iK) &&...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL) < xAxisTime(iK+1)
                                            patch([data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL)/2]-...
                                                xAxisTime(iK)/2,...
                                                [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
                                                [0 0 1],...
                                                'FaceAlpha',.5,...
                                                'EdgeColor','none')
                                        end
                                    end
                                end
                            end
                            
%                             %% Plot sacc marker
%                             % Left eye sacc
%                             % If there are left saccs for this trial
%                             if isfield(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ},'left')
%                                 % If there are blinks for this trial
%                                 if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left)
%                                     for iL = 1:length(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left)
%                                         % If the sacc time falls w/in the plot time limit plot it
%                                         if data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL) > xAxisTime(iK) &&...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL) < xAxisTime(iK+1)
%                                             patch([data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.left(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.left(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL)/2]-...
%                                                 xAxisTime(iK)/2,...
%                                                 [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
%                                                 [0 1 0],...
%                                                 'FaceAlpha',.1,...
%                                                 'EdgeColor','none')
%                                         end
%                                     end
%                                 end
%                             end
%                             % Right sacc
%                             % If there are right saccs for this trial
%                             if isfield(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ},'right')
%                                 % If there are saccs for this trial
%                                 if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right)
%                                     for iL = 1:length(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right)
%                                         % If the sacc time falls w/in the plot time limit plot it
%                                         if data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL) > xAxisTime(iK) &&...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL) < xAxisTime(iK+1)
%                                             patch([data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.right(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.right(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL)/2]-...
%                                                 xAxisTime(iK)/2,...
%                                                 [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
%                                                 [0 0 1],...
%                                                 'FaceAlpha',.1,...
%                                                 'EdgeColor','none')
%                                         end
%                                     end
%                                 end
%                             end
%                             
                            %% Plot behavioral events
                            % Actual switch times
                            if iB==1
                                for iL = 1:length(data.behav.controlSwitchTimes)
                                    % If the event time falls w/in the plot time limit plot it
                                    if data.behav.controlSwitchTimes(iL) > xAxisTime(iK) &&...
                                            data.behav.controlSwitchTimes(iL) < xAxisTime(iK+1)
                                        plot([data.behav.controlSwitchInd(iL)...
                                            data.behav.controlSwitchInd(iL)]-...
                                            xAxisTime(iK)/2,...
                                            [maxYVal -maxYVal]+xAxisMarkers(iK),...
                                            'Color',[0 0 0])
                                    end
                                end
                            end
                            
                            % Responses
                            % If there are events for this trial
                            if ~isempty(data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ})
                                for iL = 1:length(data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ})
                                    % If the event time falls w/in the plot time limit plot it
                                    if data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ}(iL) > xAxisTime(iK) &&...
                                            data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ}(iL) < xAxisTime(iK+1)
                                        plot([data.behav.(options.runType{iB}).flipTimesIdx{iI,iZ}{iJ}(iL)...
                                            data.behav.(options.runType{iB}).flipTimesIdx{iI,iZ}{iJ}(iL)]-...
                                                xAxisTime(iK)/2,...
                                            [maxYVal -maxYVal]+xAxisMarkers(iK),...
                                            'Color',[1 0 0])
                                    end
                                end
                            end
                            
                            xlim([xAxisLims(1) xAxisLims(2)])
                            ylim([0 xAxisMarkers(1)+maxYVal]);
                            set(timePlot(iJ),'XColor','k','YColor','k')
                        end
                        set(timePlot(iJ),'xtick',[],...
                            'ytick',[],...
                            'Position',timePlotSize(iJ,:))
                        if iJ==1
                            legend({'Right','Left'},'Location','NorthWest')
                        end
                    end
                    
                    %% Plot pupil size across all blocks
                    % For pupsize, we want to scale each plot based on the
                    % max value for that block, where the plot centers around max
                    % value (max-200 to max+200). 
                    % (Since, unlike for relative position, pupil size 
                    % doesn't center on 0). 
                    startYValPup = max(max(data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ}))-250;  % Min value of y axis
                    maxYValPup = 500;   % Raw 'area' values for pup size
                    xAxisMarkersPup = flip(0:maxYValPup:maxYValPup*4);
                    xAxisRangePup = [maxYValPup maxYValPup*2];
                    
                    fig2=figure();
                    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
                    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
                    figSize.figSize = [0 0 ...
                        figSize.baseSize(3)...
                        figSize.baseSize(4)];   % Size/postion of fig
                    set(gcf,'color','w')
                    set(gcf,'Position', figSize.figSize)
                    annotation('textbox',[0.42, .99, .175, 0],'edgecolor','none','string',...
                        ['Pupil Size For ', options.runType{iB}, ' Run (All Blocks)'])
                    
                    % Find length and slpit evenly into 5 rows to plot
                    for iJ=1:length(data.et.(options.runType{iB}).pupilSize{iI,iZ})   % Blocks
                        
                        % Make x-axis
                        % (Will be same x vals for all plots)
                        xAxisLims = linspace(1,length(data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ}),plotNumRows+1);
                        
                        % Find x-axis labels based on time
                        xAxisTime = linspace(1,data.et.(options.runType{iB}).fixDataTotalTime{iI,iZ}(iJ),plotNumRows+1);
                                                
                        % Plot 5 sets of lines w/ full eye position for
                        % each block
                        timePlotPup(iJ) = axes('Parent',fig2);
                        for iK=1:plotNumRows
                            
                            % Center
                            h3=plot([xAxisLims(1) xAxisLims(2)],[xAxisMarkersPup(iK) xAxisMarkersPup(iK)],'k');
                            h3.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            hold on
                            % Right eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                (data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),1) - startYValPup) + xAxisMarkersPup(iK),...
                                'Color',[0 .9 0]);
                            % Left eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                (data.et.(options.runType{iB}).pupilSize{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),2) - startYValPup) + xAxisMarkersPup(iK),...
                                'Color',[0 0 .9])   
                            
                            %% Plot blink marker
                            % Left eye blink
                            % If there are left blinks for this trial
                            if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},'left')
                                % If there are blinks for this trial
                                if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                                    for iL = 1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                                        % If the blink time falls w/in the plot time limit plot it
                                        if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL) > xAxisTime(iK) &&...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL) < xAxisTime(iK+1)
                                            
                                            % Scale x values (eyeblink times) based on xAxisTime(iK) to ensure each
                                            % point is plotted relative to x axis '0'. Scale y values relative to 
                                            % individual x axis markers.
                                            patch([data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL)/2]-...
                                                xAxisTime(iK)/2,...
                                                [maxYValPup maxYValPup 0 0]+xAxisMarkersPup(iK),...
                                                [0 1 0],...
                                                'FaceAlpha',.5,...
                                                'EdgeColor','none')
                                        end
                                    end
                                end
                            end
                            % Right eye blink
                            % If there are rightblinks for this trial
                            if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},'right')
                                % If there are blinks for this trial
                                if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                                    for iL = 1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                                        % If the blink time falls w/in the plot time limit plot it
                                        if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL) > xAxisTime(iK) &&...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL) < xAxisTime(iK+1)
                                            patch([data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL)/2]-...
                                                xAxisTime(iK)/2,...
                                                [maxYValPup maxYValPup 0 0]+xAxisMarkersPup(iK),...
                                                [0 0 1],...
                                                'FaceAlpha',.5,...
                                                'EdgeColor','none')
                                        end
                                    end
                                end
                            end
                            
                            %% Plot behavioral events
                            % Actual switch times
                            if iB==1
                                for iL = 1:length(data.behav.controlSwitchTimes)
                                    % If the event time falls w/in the plot time limit plot it
                                    if data.behav.controlSwitchTimes(iL) > xAxisTime(iK) &&...
                                            data.behav.controlSwitchTimes(iL) < xAxisTime(iK+1)
                                        plot([data.behav.controlSwitchInd(iL)...
                                            data.behav.controlSwitchInd(iL)]-...
                                            xAxisTime(iK)/2,...
                                            [maxYValPup 0]+xAxisMarkersPup(iK),...
                                            'Color',[0 0 0])
                                    end
                                end
                            end
                            
                            % Responses
                            % If there are events for this trial
                            if ~isempty(data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ})
                                for iL = 1:length(data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ})
                                    % If the event time falls w/in the plot time limit plot it
                                    if data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ}(iL) > xAxisTime(iK) &&...
                                            data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ}(iL) < xAxisTime(iK+1)
                                        plot([data.behav.(options.runType{iB}).flipTimesIdx{iI,iZ}{iJ}(iL)...
                                            data.behav.(options.runType{iB}).flipTimesIdx{iI,iZ}{iJ}(iL)]-...
                                                xAxisTime(iK)/2,...
                                            [maxYValPup 0]+xAxisMarkersPup(iK),...
                                            'Color',[1 0 0])
                                    end
                                end
                            end
                            
                            xlim([xAxisLims(1) xAxisLims(2)])
                            ylim([0 xAxisMarkersPup(1)+maxYValPup]);
                            set(timePlotPup(iJ),'XColor','k','YColor','k')
                        end
                        set(timePlotPup(iJ),'xtick',[],...
                            'ytick',[],...
                            'Position',timePupPlotSize(iJ,:))
                        if iJ==1
                            legend({'Right','Left'},'Location','NorthWest')
                        end
                    end
                    
                end
            end
        end
    end
end


%% Group analysis and plotting
% Make empty cells NaNs
for iB=1:2   % A/B runs
    for iI=1:size(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix,1)   % Subj
        for iZ=1:size(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix,2)   % B/Z days
            if isempty(data.et.(options.runType{iB}).stats.aveEucDisFromFix{iI,iZ})
                data.et.(options.runType{iB}).stats.aveEucDisFromFix{iI,iZ} = NaN;
            end
            if isempty(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix{iI,iZ})
                data.et.(options.runType{iB}).stats.corrAveEucDisFromFix{iI,iZ} = NaN;
            end
            if isempty(data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ})
                data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ} = NaN;
            end
            if isempty(data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ})
                data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ} = NaN;
            end
            %             if isempty(data.et.(options.runType{iB}).fixData_heatMap{iI,iZ})
            %                 data.et.(options.runType{iB}).fixData_heatMap{iI,iZ} = NaN;
            %             end
            %             if isempty(data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ})
            %                 data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ} = NaN;
            %             end
        end
    end
end

%% Average correction value (Distance in x/y data are corrected by)
% Take absolute value, to get more accurate average.
for iB=1:2   % A/B runs
    for iI = 1:length(data.et.(options.runType{iB}).grpIdx)
        data.et.(options.runType{iB}).stats.fullAveDistFromFix_PosX(iI) =...
            nanmean(nanmean(abs(...
            cellfun(@nanmean,data.et.(options.runType{iB}).aveDistFromFix_PosX(...
            data.et.(options.runType{iB}).grpIdx{iI},:))),2));
        data.et.(options.runType{iB}).stats.fullSteDistFromFix_PosX(iI) =...
            nanstd(nanmean(abs(...
            cellfun(@nanmean,data.et.(options.runType{iB}).aveDistFromFix_PosX(...
            data.et.(options.runType{iB}).grpIdx{iI},:))),2)) /...
            sqrt(length(data.et.(options.runType{iB}).grpIdx{iI}));
        data.et.(options.runType{iB}).stats.fullAveDistFromFix_PosY(iI) =...
            nanmean(nanmean(abs(...
            cellfun(@nanmean,data.et.(options.runType{iB}).aveDistFromFix_PosY(...
            data.et.(options.runType{iB}).grpIdx{iI},:))),2));
        data.et.(options.runType{iB}).stats.fullSteDistFromFix_PosY(iI) =...
            nanstd(nanmean(abs(...
            cellfun(@nanmean,data.et.(options.runType{iB}).aveDistFromFix_PosY(...
            data.et.(options.runType{iB}).grpIdx{iI},:))),2)) /...
            sqrt(length(data.et.(options.runType{iB}).grpIdx{iI}));
    end
end

%% Run stats on correction
% X distance
% Create grouping variables for ANOVA
% Combine the blocks into a subjxdayxblock array
clear all_data all_subj all_group all_block grouping
nest = zeros(3,3);
nest(1,2) = 1;
if options.displayETFigs==1
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
else
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
end
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).aveDistFromFix_PosX,1)
        for iZ=1:size(data.et.(options.runType{iB}).aveDistFromFix_PosX,2)
            for iJ=1:5
                if ~isempty(data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ})
                    all_data(iI,iZ,1:5) = data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ}(:);
                else
                    all_data(iI,iZ,1:5) = NaN;
                end
                all_block(iI,iZ,iJ) = iJ;
                all_group(data.et.(options.runType{iB}).grpIdx{1},iZ,iJ) = 1;
                all_group(data.et.(options.runType{iB}).grpIdx{2},iZ,iJ) = 2;
                all_group(data.et.(options.runType{iB}).grpIdx{3},iZ,iJ) = 3;
                all_subj(iI,iZ,iJ) = iI;
            end
        end
    end
    
    grouping = all_group(:,1,1);
    
    %% ANOVA
    [data.et.(options.runType{iB}).stats.AnovaP_distFromFix_posX,...
        data.et.(options.runType{iB}).stats.AnovaTable_distFromFix_posX] =...
        anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig_anova);
    
    %% K-W
    [data.et.(options.runType{iB}).stats.KW_P_distFromFix_posX,...
        data.et.(options.runType{iB}).stats.KW_Table_distFromFix_posX, ...
        data.et.(options.runType{iB}).stats.KW_Stats_distFromFix_posX] = ...
        kruskalwallis(nanmean(squeeze(nanmean(all_data,2)),2), grouping, show_stats_fig_kw);
    
    clear all_data all_subj all_group all_block grouping
end

% Y distance
% Create grouping variables for ANOVA
% Combine the blocks into a subjxdayxblock array
clear all_data all_subj all_group all_block
nest = zeros(3,3);
nest(1,2) = 1;
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).aveDistFromFix_PosY,1)
        for iZ=1:size(data.et.(options.runType{iB}).aveDistFromFix_PosY,2)
            for iJ=1:5
                if ~isempty(data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ})
                    all_data(iI,iZ,1:5) = data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ}(:);
                else
                    all_data(iI,iZ,1:5) = NaN;
                end
                all_block(iI,iZ,iJ) = iJ;
                all_group(data.et.(options.runType{iB}).grpIdx{1},iZ,iJ) = 1;
                all_group(data.et.(options.runType{iB}).grpIdx{2},iZ,iJ) = 2;
                all_group(data.et.(options.runType{iB}).grpIdx{3},iZ,iJ) = 3;
                all_subj(iI,iZ,iJ) = iI;
            end
        end
    end
    
    grouping = all_group(:,1,1);
    
    %% ANOVA
    [data.et.(options.runType{iB}).stats.AnovaP_distFromFix_posY,...
        data.et.(options.runType{iB}).stats.AnovaTable_distFromFix_posY] =...
        anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig_anova);
    
    %% K-W
    [data.et.(options.runType{iB}).stats.KW_P_distFromFix_posY,...
        data.et.(options.runType{iB}).stats.KW_Table_distFromFix_posY, ...
        data.et.(options.runType{iB}).stats.KW_Stats_distFromFix_posY] = ...
        kruskalwallis(nanmean(squeeze(nanmean(all_data,2)),2), grouping, show_stats_fig_kw);
    
    clear all_data all_subj all_group all_block grouping
end

%% Plot correction value
if options.displayETFigs
    figure()
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(4)*...
        (figSize.aspectRatio(1)/figSize.aspectRatio(2))...
        figSize.baseSize(4)];   % Size/postion of fig
    
    %% Bar graphs
    %% Plot for X
    subplot(2,1,1)
    hb = bar([data.et.A.stats.fullAveDistFromFix_PosX;...
        data.et.B.stats.fullAveDistFromFix_PosX] ./ options.PPD);
    hold on
    
    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(abs(...
            cellfun(@nanmean,data.et.A.aveDistFromFix_PosX(...
            data.et.A.grpIdx{iI},:))),2) ./ options.PPD,...
            nanmean(abs(...
            cellfun(@nanmean,data.et.B.aveDistFromFix_PosX(...
            data.et.B.grpIdx{iI},:))),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.A.stats.fullAveDistFromFix_PosX(iI);...
            data.et.B.stats.fullAveDistFromFix_PosX(iI)] ./ options.PPD,...
            [data.et.A.stats.fullSteDistFromFix_PosX(iI);...
            data.et.B.stats.fullSteDistFromFix_PosX(iI)] ./ options.PPD,'.k');
    end
    
    % Plot the significance
    text(.75,4,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.A.stats.KW_Table_distFromFix_posX{2,3},')=',...
        data.et.A.stats.KW_Table_distFromFix_posX{2,5},...
        'p=',data.et.A.stats.KW_Table_distFromFix_posX{2,6}));
    text(1.75,4,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_distFromFix_posX{2,3},')=',...
        data.et.B.stats.KW_Table_distFromFix_posX{2,5},...
        'p=',data.et.B.stats.KW_Table_distFromFix_posX{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Real Switch Task','Bi-stable Switch Task'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 10],'ytick',[0.1 0.5 1 2.5 5 10])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};

    box off
    ylabel('Offset (°)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Average X correction value used (in absolute °)','fontsize',20)
    legend(options.(['group_labels_', options.runType{iB}]))
    
    %% Plot for Y
    subplot(2,1,2)
    hb2 = bar([data.et.A.stats.fullAveDistFromFix_PosY;...
        data.et.B.stats.fullAveDistFromFix_PosY] ./ options.PPD);
    hold on
    
    % Beeswarm
    for iI=1:numel(hb2)
        x_val = hb2(iI).XData + hb2(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(abs(...
            cellfun(@nanmean,data.et.A.aveDistFromFix_PosY(...
            data.et.A.grpIdx{iI},:))),2) ./ options.PPD,...
            nanmean(abs(...
            cellfun(@nanmean,data.et.B.aveDistFromFix_PosY(...
            data.et.B.grpIdx{iI},:))),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.A.stats.fullAveDistFromFix_PosY(iI);...
            data.et.B.stats.fullAveDistFromFix_PosY(iI)] ./ options.PPD,...
            [data.et.A.stats.fullSteDistFromFix_PosY(iI);...
            data.et.B.stats.fullSteDistFromFix_PosY(iI)] ./ options.PPD,'.k');
    end
    
    % Plot the significance
    text(.75,10,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.A.stats.KW_Table_distFromFix_posY{2,3},')=',...
        data.et.A.stats.KW_Table_distFromFix_posY{2,5},...
        'p=',data.et.A.stats.KW_Table_distFromFix_posY{2,6}));
    text(1.75,10,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_distFromFix_posY{2,3},')=',...
        data.et.B.stats.KW_Table_distFromFix_posY{2,5},...
        'p=',data.et.B.stats.KW_Table_distFromFix_posY{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Real Switch Task','Bi-stable Switch Task'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 15],'ytick',[0.1 0.5 1 2.5 5 10])
    set(hb2,'linewidth',2)
    hb2(1).FaceColor = col_list{1};
    hb2(2).FaceColor = col_list{2};
    hb2(3).FaceColor = col_list{3};
    
    box off
    ylabel('Offset (°)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Average Y correction value used (in absolute °)','fontsize',20)
end

%% Average and variance of corrected/uncorrected euclidian distance away from fixation
% Average distance away from fixation across groups
% Average across B/Z
% Average distance within groups
for iB=1:2   % A/B runs
    for iI = 1:length(data.et.(options.runType{iB}).grpIdx)   % Group
        %% Average distance
        % Uncorrected
        data.et.(options.runType{iB}).stats.grpAveDistFromFix(iI) =...
            nanmean(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.aveEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2));
        data.et.(options.runType{iB}).stats.grpSteDistFromFix(iI) =...
            nanstd(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.aveEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2)) ./...
            sqrt(length(data.et.(options.runType{iB}).grpIdx{iI}));
        
        % Corrected
        data.et.(options.runType{iB}).stats.grpAveCorrDistFromFix(iI) =...
            nanmean(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.corrAveEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2));
        data.et.(options.runType{iB}).stats.grpSteCorrDistFromFix(iI) =...
            nanstd(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.corrAveEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2)) ./...
            sqrt(length(data.et.(options.runType{iB}).grpIdx{iI}));
    
        %% Std of distance
        data.et.(options.runType{iB}).stats.grpAveVarFromFix(iI) =...
            nanmean(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.stdEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2));
        data.et.(options.runType{iB}).stats.grpSteVarFromFix(iI) =...
            nanstd(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.stdEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2)) ./...
            sqrt(length(data.et.(options.runType{iB}).grpIdx{iI}));
        
        % Corrected
        data.et.(options.runType{iB}).stats.grpAveCorrVarFromFix(iI) =...
            nanmean(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.corrStdEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2));
        data.et.(options.runType{iB}).stats.grpSteCorrVarFromFix(iI) =...
            nanstd(nanmean(cellfun(@nanmean,data.et.(options.runType{iB}).stats.corrStdEucDisFromFix(...
            data.et.(options.runType{iB}).grpIdx{iI},:)),2)) ./...
            sqrt(length(data.et.(options.runType{iB}).grpIdx{iI}));
    
    end
end


%% Run stats on average Euc distance
% Uncorrected Euc distance
% Create grouping variables for ANOVA
% Combine the blocks into a subjxdayxblock array
clear all_data all_subj all_group all_block grouping
nest = zeros(3,3);
nest(1,2) = 1;
if options.displayETFigs==1
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
else
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
end
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).stats.aveEucDisFromFix,1)
        for iZ=1:size(data.et.(options.runType{iB}).stats.aveEucDisFromFix,2)
            for iJ=1:5
                if ~isempty(data.et.(options.runType{iB}).stats.aveEucDisFromFix{iI,iZ})
                    all_data(iI,iZ,1:5) = data.et.(options.runType{iB}).stats.aveEucDisFromFix{iI,iZ}(:);
                else
                    all_data(iI,iZ,1:5) = NaN;
                end
                all_block(iI,iZ,iJ) = iJ;
                all_group(data.et.(options.runType{iB}).grpIdx{1},iZ,iJ) = 1;
                all_group(data.et.(options.runType{iB}).grpIdx{2},iZ,iJ) = 2;
                all_group(data.et.(options.runType{iB}).grpIdx{3},iZ,iJ) = 3;
                all_subj(iI,iZ,iJ) = iI;
            end
        end
    end
    
    grouping = all_group(:,1,1);
    
    %% ANOVA
    [data.et.(options.runType{iB}).stats.AnovaP_eucDisAve,...
        data.et.(options.runType{iB}).stats.AnovaTable_eucDisAve] =...
        anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig_anova);
        
    %% K-W
    [data.et.(options.runType{iB}).stats.KW_P_eucDisAve,...
        data.et.(options.runType{iB}).stats.KW_Table_eucDisAve, ...
        data.et.(options.runType{iB}).stats.KW_Stats_eucDisAve] = ...
        kruskalwallis(nanmean(squeeze(nanmean(all_data,2)),2), grouping, show_stats_fig_kw);
    
    clear all_data all_subj all_group all_block grouping
end

% Corrected Euc distance
% Create grouping variables for ANOVA
% Combine the blocks into a subjxdayxblock array
clear all_data all_subj all_group all_block grouping
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix,1)
        for iZ=1:size(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix,2)
            for iJ=1:5
                if ~isempty(data.et.(options.runType{iB}).stats.corrAveEucDisFromFix{iI,iZ})
                    all_data(iI,iZ,1:5) = data.et.(options.runType{iB}).stats.corrAveEucDisFromFix{iI,iZ}(:);
                else
                    all_data(iI,iZ,1:5) = NaN;
                end
                all_block(iI,iZ,iJ) = iJ;
                all_group(data.et.(options.runType{iB}).grpIdx{1},iZ,iJ) = 1;
                all_group(data.et.(options.runType{iB}).grpIdx{2},iZ,iJ) = 2;
                all_group(data.et.(options.runType{iB}).grpIdx{3},iZ,iJ) = 3;
                all_subj(iI,iZ,iJ) = iI;
            end
        end
    end
    
    grouping = all_group(:,1,1);
    
    %% ANOVA
    [data.et.(options.runType{iB}).stats.AnovaP_corrEucDisAve,...
        data.et.(options.runType{iB}).stats.AnovaTable_corrEucDisAve] =...
        anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig_anova);
    
    %% K-W
    [data.et.(options.runType{iB}).stats.KW_P_corrEucDisAve,...
        data.et.(options.runType{iB}).stats.KW_Table_corrEucDisAve, ...
        data.et.(options.runType{iB}).stats.KW_Stats_corrEucDis] = ...
        kruskalwallis(nanmean(squeeze(nanmean(all_data,2)),2), grouping, show_stats_fig_kw);
    
    clear all_data all_subj all_group all_block grouping
end

%% Run stats on average of Euc variance
% Uncorrected Euc variance
% Create grouping variables for ANOVA
% Combine the blocks into a subjxdayxblock array
clear all_data all_subj all_group all_block grouping
nest = zeros(3,3);
nest(1,2) = 1;
if options.displayETFigs==1
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
else
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
end
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).stats.stdEucDisFromFix,1)
        for iZ=1:size(data.et.(options.runType{iB}).stats.stdEucDisFromFix,2)
            for iJ=1:5
                if ~isempty(data.et.(options.runType{iB}).stats.stdEucDisFromFix{iI,iZ})
                    all_data(iI,iZ,1:5) = data.et.(options.runType{iB}).stats.stdEucDisFromFix{iI,iZ}(:);
                else
                    all_data(iI,iZ,1:5) = NaN;
                end
                all_block(iI,iZ,iJ) = iJ;
                all_group(data.et.(options.runType{iB}).grpIdx{1},iZ,iJ) = 1;
                all_group(data.et.(options.runType{iB}).grpIdx{2},iZ,iJ) = 2;
                all_group(data.et.(options.runType{iB}).grpIdx{3},iZ,iJ) = 3;
                all_subj(iI,iZ,iJ) = iI;
            end
        end
    end
    
    grouping = all_group(:,1,1);
    
    %% ANOVA
    [data.et.(options.runType{iB}).stats.AnovaP_eucVarAve,...
        data.et.(options.runType{iB}).stats.AnovaTable_eucVarAve] =...
        anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig_anova);
        
    %% K-W
    [data.et.(options.runType{iB}).stats.KW_P_eucVarAve,...
        data.et.(options.runType{iB}).stats.KW_Table_eucVarAve, ...
        data.et.(options.runType{iB}).stats.KW_Stats_eucVarAve] = ...
        kruskalwallis(nanmean(squeeze(nanmean(all_data,2)),2), grouping, show_stats_fig_kw);
    
    clear all_data all_subj all_group all_block grouping
end

% Corrected Euc variance
% Create grouping variables for ANOVA
% Combine the blocks into a subjxdayxblock array
clear all_data all_subj all_group all_block grouping
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).stats.corrStdEucDisFromFix,1)
        for iZ=1:size(data.et.(options.runType{iB}).stats.corrStdEucDisFromFix,2)
            for iJ=1:5
                if ~isempty(data.et.(options.runType{iB}).stats.corrStdEucDisFromFix{iI,iZ})
                    all_data(iI,iZ,1:5) = data.et.(options.runType{iB}).stats.corrStdEucDisFromFix{iI,iZ}(:);
                else
                    all_data(iI,iZ,1:5) = NaN;
                end
                all_block(iI,iZ,iJ) = iJ;
                all_group(data.et.(options.runType{iB}).grpIdx{1},iZ,iJ) = 1;
                all_group(data.et.(options.runType{iB}).grpIdx{2},iZ,iJ) = 2;
                all_group(data.et.(options.runType{iB}).grpIdx{3},iZ,iJ) = 3;
                all_subj(iI,iZ,iJ) = iI;
            end
        end
    end
    
    grouping = all_group(:,1,1);
    
    %% ANOVA
    [data.et.(options.runType{iB}).stats.AnovaP_corrEucVarAve,...
        data.et.(options.runType{iB}).stats.AnovaTable_corrEucVarAve] =...
        anovan(all_data(:),{all_subj(:),...
        all_group(:),all_block(:)},'random',1,'continuous',[3],...
        'nested',nest,'model','full','varnames',{'subj','group','block'},...
        'display',show_stats_fig_anova);
    
    %% K-W
    [data.et.(options.runType{iB}).stats.KW_P_corrEucVarAve,...
        data.et.(options.runType{iB}).stats.KW_Table_corrEucVarAve, ...
        data.et.(options.runType{iB}).stats.KW_Stats_corrEucVarAve] = ...
        kruskalwallis(nanmean(squeeze(nanmean(all_data,2)),2), grouping, show_stats_fig_kw);
    
    clear all_data all_subj all_group all_block grouping
end


%% Plot average of euc distance
if options.displayETFigs
    figure()
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(4)*...
        (figSize.aspectRatio(1)/figSize.aspectRatio(2))...
        figSize.baseSize(4)];   % Size/postion of fig
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    
    %% Bar graphs
    %% Uncorrected euc distance
    subplot(2,1,1)
    hb = bar([data.et.A.stats.grpAveDistFromFix;...
        data.et.B.stats.grpAveDistFromFix] ./ options.PPD);
    hold on
    
    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(cellfun(@nanmean,data.et.A.stats.aveEucDisFromFix(...
            data.et.A.grpIdx{iI},:)),2) ./ options.PPD,...
            nanmean(cellfun(@nanmean,data.et.B.stats.aveEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errobars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.A.stats.grpAveDistFromFix(iI);...
            data.et.B.stats.grpAveDistFromFix(iI)] ./ options.PPD,...
            [data.et.A.stats.grpSteDistFromFix(iI);...
            data.et.B.stats.grpSteDistFromFix(iI)] ./ options.PPD,...
            '.k','LineWidth',2);
    end
    
    % Plot the significance
    text(.75,10,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.A.stats.KW_Table_eucDisAve{2,3},')=',...
        data.et.A.stats.KW_Table_eucDisAve{2,5},...
        'p=',data.et.A.stats.KW_Table_eucDisAve{2,6}));
    text(1.75,10,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_eucDisAve{2,3},')=',...
        data.et.B.stats.KW_Table_eucDisAve{2,5},...
        'p=',data.et.B.stats.KW_Table_eucDisAve{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Real Switch Task','Bi-stable Switch Task'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 15],'ytick',[0.1 0.5 1 2.5 5 10])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};
    set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Offset (°)','fontsize',15)
    set(gca,'XColor','k','YColor','k')
    title('Average of Euclidian distance (in absolute °)','fontsize',20)
    legend(options.(['group_labels_', options.runType{iB}]))
    
    %% Corrected euc distance
    subplot(2,1,2)
    hb = bar([data.et.A.stats.grpAveCorrDistFromFix;...
        data.et.B.stats.grpAveCorrDistFromFix] ./ options.PPD);
    hold on

    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(cellfun(@nanmean,data.et.A.stats.corrAveEucDisFromFix(...
            data.et.A.grpIdx{iI},:)),2) ./ options.PPD,...
            nanmean(cellfun(@nanmean,data.et.B.stats.corrAveEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.A.stats.grpAveCorrDistFromFix(iI);...
            data.et.B.stats.grpAveCorrDistFromFix(iI)] ./ options.PPD,...
            [data.et.A.stats.grpSteCorrDistFromFix(iI);...
            data.et.B.stats.grpSteCorrDistFromFix(iI)] ./ options.PPD,...
            '.k','LineWidth',2);
    end
    
    % Plot the significance
    text(.75,6,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.A.stats.KW_Table_corrEucDisAve{2,3},')=',...
        data.et.A.stats.KW_Table_corrEucDisAve{2,5},...
        'p=',data.et.A.stats.KW_Table_corrEucDisAve{2,6}));
    text(1.75,6,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_corrEucDisAve{2,3},')=',...
        data.et.B.stats.KW_Table_corrEucDisAve{2,5},...
        'p=',data.et.B.stats.KW_Table_corrEucDisAve{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Real Switch Task','Bi-stable Switch Task'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 15],'ytick',[0.1 0.5 1 2.5 5 10])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};
    set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Offset (°)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Average of corrected Euclidian distance (in absolute °)','fontsize',20)
end


%% Plot average of euc variance
if options.displayETFigs
    figure()
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(4)*...
        (figSize.aspectRatio(1)/figSize.aspectRatio(2))...
        figSize.baseSize(4)];   % Size/postion of fig
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    
    %% Bar graphs
    %% Uncorrected euc variance
    subplot(2,1,1)
    hb = bar([data.et.A.stats.grpAveVarFromFix;...
        data.et.B.stats.grpAveVarFromFix] ./ options.PPD);
    hold on
    
    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(cellfun(@nanmean,data.et.A.stats.stdEucDisFromFix(...
            data.et.A.grpIdx{iI},:)),2) ./ options.PPD,...
            nanmean(cellfun(@nanmean,data.et.B.stats.stdEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errobars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.A.stats.grpAveVarFromFix(iI);...
            data.et.B.stats.grpAveVarFromFix(iI)] ./ options.PPD,...
            [data.et.A.stats.grpSteVarFromFix(iI);...
            data.et.B.stats.grpSteVarFromFix(iI)] ./ options.PPD,...
            '.k','LineWidth',2);
    end
    
    % Plot the significance
    text(.75,4,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.A.stats.KW_Table_eucVarAve{2,3},')=',...
        data.et.A.stats.KW_Table_eucVarAve{2,5},...
        'p=',data.et.A.stats.KW_Table_eucVarAve{2,6}));
    text(1.75,4,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_eucVarAve{2,3},')=',...
        data.et.B.stats.KW_Table_eucVarAve{2,5},...
        'p=',data.et.B.stats.KW_Table_eucVarAve{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Real Switch Task','Bi-stable Switch Task'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 5],'ytick',[0.1 0.5 1 2.5 5])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};
    set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Stability (SD in °)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Average of Euclidian variance (in absolute °)','fontsize',20)
    legend(options.(['group_labels_', options.runType{iB}]))
    
    %% Corrected euc variance
    subplot(2,1,2)
    hb = bar([data.et.A.stats.grpAveCorrVarFromFix;...
        data.et.B.stats.grpAveCorrVarFromFix] ./ options.PPD);
    hold on

    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(cellfun(@nanmean,data.et.A.stats.corrStdEucDisFromFix(...
            data.et.A.grpIdx{iI},:)),2) ./ options.PPD,...
            nanmean(cellfun(@nanmean,data.et.B.stats.corrStdEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.A.stats.grpAveCorrVarFromFix(iI);...
            data.et.B.stats.grpAveCorrVarFromFix(iI)] ./ options.PPD,...
            [data.et.A.stats.grpSteCorrVarFromFix(iI);...
            data.et.B.stats.grpSteCorrVarFromFix(iI)] ./ options.PPD,...
            '.k','LineWidth',2);
    end
    
    % Plot the significance
    text(.75,3,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.A.stats.KW_Table_corrEucVarAve{2,3},')=',...
        data.et.A.stats.KW_Table_corrEucVarAve{2,5},...
        'p=',data.et.A.stats.KW_Table_corrEucVarAve{2,6}));
    text(1.75,3,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_corrEucVarAve{2,3},')=',...
        data.et.B.stats.KW_Table_corrEucVarAve{2,5},...
        'p=',data.et.B.stats.KW_Table_corrEucVarAve{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Real Switch Task','Bi-stable Switch Task'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 5],'ytick',[0.1 0.5 1 2.5 5])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};
    set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Stability (SD in °)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Average of corrected Euclidian variance (in absolute °)','fontsize',20)
end


%% Plot Euc var/distance for both uncorr/corr bistable task only (for paper)
if options.displayETFigs
    figure()
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(4)*...
        (figSize.aspectRatio(1)/figSize.aspectRatio(2))...
        figSize.baseSize(4)];   % Size/postion of fig
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    
    %% Bar graphs
    %% Corr/Unc Euc distance (bistable)
    subplot(2,1,1)
    hb = bar([data.et.B.stats.grpAveDistFromFix;...
        data.et.B.stats.grpAveCorrDistFromFix] ./ options.PPD);
    hold on
    
    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(cellfun(@nanmean,data.et.B.stats.aveEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD,...
            nanmean(cellfun(@nanmean,data.et.B.stats.corrAveEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errobars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.B.stats.grpAveDistFromFix(iI);...
            data.et.B.stats.grpAveCorrDistFromFix(iI)] ./ options.PPD,...
            [data.et.B.stats.grpSteDistFromFix(iI);...
            data.et.B.stats.grpSteCorrDistFromFix(iI)] ./ options.PPD,...
            '.k','LineWidth',2);
    end
    
%     % Plot the significance
%     text(.75,4,...
%         sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_eucDisAve{2,3},')=',...
%         data.et.B.stats.KW_Table_eucDisAve{2,5},...
%         'p=',data.et.B.stats.KW_Table_eucDisAve{2,6}));
%     text(1.75,4,...
%         sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_corrEucDisAve{2,3},')=',...
%         data.et.B.stats.KW_Table_corrEucDisAve{2,5},...
%         'p=',data.et.B.stats.KW_Table_corrEucDisAve{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Uncorrected','Corrected'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 15],'ytick',[0.1 0.5 1 2.5 5 10])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};
    set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Offset (°)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Euclidian Distance for Bi-stable Task (in absolute °)','fontsize',20)
    
    %% Corr/Unc Euc variance (bistable)
    subplot(2,1,2)
    hb = bar([data.et.B.stats.grpAveVarFromFix;...
        data.et.B.stats.grpAveCorrVarFromFix] ./ options.PPD);
    hold on

    % Beeswarm
    for iI=1:numel(hb)
        x_val = hb(iI).XData + hb(iI).XOffset;
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({nanmean(cellfun(@nanmean,data.et.B.stats.stdEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD,...
            nanmean(cellfun(@nanmean,data.et.B.stats.corrStdEucDisFromFix(...
            data.et.B.grpIdx{iI},:)),2) ./ options.PPD},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    ngroups = 2;
    nbars = 3;
    groupwidth = min(nbars/(nbars+1.5));
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,[data.et.B.stats.grpAveVarFromFix(iI);...
            data.et.B.stats.grpAveCorrVarFromFix(iI)] ./ options.PPD,...
            [data.et.B.stats.grpSteVarFromFix(iI);...
            data.et.B.stats.grpSteCorrVarFromFix(iI)] ./ options.PPD,...
            '.k','LineWidth',2);
    end
    
    % Plot the significance
%     text(.75,3,...
%         sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_eucVarAve{2,3},')=',...
%         data.et.B.stats.KW_Table_eucVarAve{2,5},...
%         'p=',data.et.B.stats.KW_Table_eucVarAve{2,6}));
%     text(1.75,3,...
%         sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.KW_Table_corrEucVarAve{2,3},')=',...
%         data.et.B.stats.KW_Table_corrEucVarAve{2,5},...
%         'p=',data.et.B.stats.KW_Table_corrEucVarAve{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Uncorrected','Corrected'},'fontsize',15)
    set(gca,'YScale','log')
    set(gca,'ylim',[0 5],'ytick',[0.1 0.5 1 2.5 5])
    set(hb,'linewidth',2)
    hb(1).FaceColor = col_list{1};
    hb(2).FaceColor = col_list{2};
    hb(3).FaceColor = col_list{3};
    set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Offset (°)','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    title('Euclidian Variance for Bi-stable Task (in absolute °)','fontsize',20)
end


%% Average total number of eye blinks/saccs
% Average across blocks and days for each particpant for each run
for iB=1:2   % A/B runs
    for iI=1:size(data.et.(options.runType{iB}).blinkNumber,1)
        % Blinks
        data.et.(options.runType{iB}).stats.numBlinks_subj(1,iI) =...
            nanmean([data.et.(options.runType{iB}).blinkNumber(iI,:,:).left]);
        data.et.(options.runType{iB}).stats.numBlinks_subj(2,iI) =...
            nanmean([data.et.(options.runType{iB}).blinkNumber(iI,:,:).right]);
        data.et.(options.runType{iB}).stats.numBlinks_subj(3,iI) =...
            nanmean([data.et.(options.runType{iB}).blinkNumber(iI,:,:).both]);
        
        % Saccs
        data.et.(options.runType{iB}).stats.numSaccs_subj(1,iI) =...
            nanmean([data.et.(options.runType{iB}).saccNumber(iI,:,:).left]);
        data.et.(options.runType{iB}).stats.numSaccs_subj(2,iI) =...
            nanmean([data.et.(options.runType{iB}).saccNumber(iI,:,:).right]);
        data.et.(options.runType{iB}).stats.numSaccs_subj(3,iI) =...
            nanmean([data.et.(options.runType{iB}).saccNumber(iI,:,:).both]);
    end
   
    
    % For blinks
    for iE = 1:3   % for either eye events and combined eye events
        % Average across all participants
        % Blinks
        data.et.(options.runType{iB}).stats.numBlinks_fullAve(iE) = ...
            squeeze(nanmean(data.et.(options.runType{iB}).stats.numBlinks_subj(iE,:)));
        data.et.(options.runType{iB}).stats.numBlinks_fullSte(iE) = ...
            squeeze(nanstd(data.et.(options.runType{iB}).stats.numBlinks_subj(iE,:))) ./ ...
            sqrt(size(data.et.(options.runType{iB}).blinkNumber,1));
        for iG = 1:3   % Across the 3 groups
            % Blinks
            data.et.(options.runType{iB}).stats.numBlinks_grpAve(iE,iG) = ...
                squeeze(nanmean(data.et.(options.runType{iB}).stats.numBlinks_subj(iE,...
                data.et.(options.runType{iB}).grpIdx{iG})));
            data.et.(options.runType{iB}).stats.numBlinks_grpSte(iE,iG) = ...
                squeeze(nanstd(data.et.(options.runType{iB}).stats.numBlinks_subj(iE,...
                data.et.(options.runType{iB}).grpIdx{iG}))) ./ ...
                sqrt(length(data.et.(options.runType{iB}).grpIdx{iG}));
        end
    end
    
    
    % For saccs
    for iE = 1:3   % for both eyes
        % Average across all participants
        % Saccs
        data.et.(options.runType{iB}).stats.numSaccs_fullAve(iE) = ...
            squeeze(nanmean(data.et.(options.runType{iB}).stats.numSaccs_subj(iE,:)));
        data.et.(options.runType{iB}).stats.numSaccs_fullSte(iE) = ...
            squeeze(nanstd(data.et.(options.runType{iB}).stats.numSaccs_subj(iE,:))) ./ ...
            sqrt(size(data.et.(options.runType{iB}).blinkNumber,1));
        for iG = 1:3   % Across the 3 groups
            % Saccs
            data.et.(options.runType{iB}).stats.numSaccs_grpAve(iE,iG) = ...
                squeeze(nanmean(data.et.(options.runType{iB}).stats.numSaccs_subj(iE,...
                data.et.(options.runType{iB}).grpIdx{iG})));
            data.et.(options.runType{iB}).stats.numSaccs_grpSte(iE,iG) = ...
                squeeze(nanstd(data.et.(options.runType{iB}).stats.numSaccs_subj(iE,...
                data.et.(options.runType{iB}).grpIdx{iG}))) ./ ...
                sqrt(length(data.et.(options.runType{iB}).grpIdx{iG}));
        end
    end
end

%% Run stats on ave num blinks/saccades
% Look at KW test between number of blinks in each group
if options.displayETFigs==1
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
else
    show_stats_fig_kw = 'off';
    show_stats_fig_anova = 'off';
end

clear grouping
grouping = zeros([length(data.et.(options.runType{iB}).subjID) 1]);
for iG=1:3
    grouping(data.et.(options.runType{iB}).grpIdx{iG}) = iG;
end

% 3KW test between number of blinks that occured in both eyes
[data.et.B.stats.numBlinks_kw_P,...
    data.et.B.stats.numBlinks_kw_Table,...
    data.et.B.stats.numBlinks_kw_Stats] = ...
    kruskalwallis(data.et.(options.runType{iB}).stats.numBlinks_subj(3,:),...
    grouping,show_stats_fig_kw);

% 3KW test between number of saccs that occured in both eyes
[data.et.B.stats.numSaccs_kw_P,...
    data.et.B.stats.numSaccs_kw_Table,...
    data.et.B.stats.numSaccs_kw_Stats] = ...
    kruskalwallis(data.et.(options.runType{iB}).stats.numSaccs_subj(3,:),...
    grouping,show_stats_fig_kw);

% Look at correlation between blink number and switch rate
% First average switch rate across block, then days
data.behav.A.Hz_flips_aveFull(:) = squeeze(nanmean(squeeze(nanmean(data.behav.A.Hz_flips,3)),2));
data.behav.B.Hz_flips_aveFull(:) = squeeze(nanmean(squeeze(nanmean(data.behav.B.Hz_flips,3)),2));

% Correlate between switch rate and blink number
[data.et.B.stats.corr_switch_blinkLeft_r, data.et.B.stats.corr_switch_blinkLeft_p] = ...
    corr(data.behav.B.Hz_flips_aveFull', data.et.B.stats.numBlinks_subj(1,:)',...
    'rows','complete','type','Spearman');
[data.et.B.stats.corr_switch_blinkRight_r, data.et.B.stats.corr_switch_blinkRight_p] = ...
    corr(data.behav.B.Hz_flips_aveFull', data.et.B.stats.numBlinks_subj(2,:)',...
    'rows','complete','type','Spearman');
[data.et.B.stats.corr_switch_blinkBoth_r, data.et.B.stats.corr_switch_blinkBoth_p] = ...
    corr(data.behav.B.Hz_flips_aveFull', data.et.B.stats.numBlinks_subj(3,:)',...
    'rows','complete','type','Spearman');

% Correlate between switch rate and sacc number
[data.et.B.stats.corr_switch_saccLeft_r, data.et.B.stats.corr_switch_saccLeft_p] = ...
    corr(log10(data.behav.B.Hz_flips_aveFull)', data.et.B.stats.numSaccs_subj(1,:)',...
    'rows','complete','type','Spearman');
[data.et.B.stats.corr_switch_saccRight_r, data.et.B.stats.corr_switch_saccRight_p] = ...
    corr(log10(data.behav.B.Hz_flips_aveFull)', data.et.B.stats.numSaccs_subj(2,:)',...
    'rows','complete','type','Spearman');
[data.et.B.stats.corr_switch_saccBoth_r, data.et.B.stats.corr_switch_saccBoth_p] = ...
    corr(log10(data.behav.B.Hz_flips_aveFull)', data.et.B.stats.numSaccs_subj(3,:)',...
    'rows','complete','type','Spearman');


%% Plot average number eye movements and blinks for left/right/both
if options.displayETFigs
    figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.75 ...
        figSize.baseSize(4)];   % Size/postion of fig
    set(gcf,'color','w')
    set(gcf,'Position', figSize.figSize)
    
    %% Blinks
    subplot(2,1,1)
    hb = bar([[data.et.B.stats.numBlinks_fullAve(1) data.et.B.stats.numBlinks_grpAve(1,:)];...
        [data.et.B.stats.numBlinks_fullAve(2) data.et.B.stats.numBlinks_grpAve(2,:)];
        [data.et.B.stats.numBlinks_fullAve(3) data.et.B.stats.numBlinks_grpAve(3,:)]]);
    hold on
    
    % Beeswarm
    % For average
    x_val = (hb(1).XData + hb(1).XOffset);
    bee_bin_width = .1;
    bee_spread_width = .2;
    beePlot = plotSpread({data.et.B.stats.numBlinks_subj(1,:),...
        data.et.B.stats.numBlinks_subj(2,:),...
        data.et.B.stats.numBlinks_subj(3,:)},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.5 .5 .5]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    % For groups
    for iI=1:numel(hb)-1
        x_val = (hb(iI+1).XData + hb(iI+1).XOffset);
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{iI}),...
            data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{iI}),...
            data.et.B.stats.numBlinks_subj(3,data.et.B.grpIdx{iI})},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    ngroups = 3;
    nbars = 4;
    groupwidth = min(nbars/(nbars+1.5));
    holderMean =[[data.et.B.stats.numBlinks_fullAve(1) data.et.B.stats.numBlinks_grpAve(1,:)];...
        [data.et.B.stats.numBlinks_fullAve(2) data.et.B.stats.numBlinks_grpAve(2,:)];...
        [data.et.B.stats.numBlinks_fullAve(3) data.et.B.stats.numBlinks_grpAve(3,:)]];
    holderSte = [[data.et.B.stats.numBlinks_fullSte(1) data.et.B.stats.numBlinks_grpSte(1,:)];...
        [data.et.B.stats.numBlinks_fullSte(2) data.et.B.stats.numBlinks_grpSte(2,:)];...
        [data.et.B.stats.numBlinks_fullSte(3) data.et.B.stats.numBlinks_grpSte(3,:)]];
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,holderMean(:,iI),holderSte(:,iI),'.k','LineWidth',2);
    end
    
    % Plot the significance
    text(2.9,275,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.numBlinks_kw_Table{2,3},')=',...
        data.et.B.stats.numBlinks_kw_Table{2,5},...
        'p=',data.et.B.stats.numBlinks_kw_Table{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Left Eye','Right Eye','Comb. Eyes'},'fontsize',15)
    set(hb,'linewidth',2)
    hb(1).FaceColor = [0 0 0];
    hb(2).FaceColor = col_list{1};
    hb(3).FaceColor = col_list{2};
    hb(4).FaceColor = col_list{3};
%     set(gca,'XTick',1:3,'fontsize',15)
    
    box off
    ylabel('Num Blinks','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
        options.(['group_labels_', options.runType{iB}])],'FontSize',8)
    title('Average Number of Blinks','fontsize',20)
    
    %% Saccades
    subplot(2,1,2)
    hb = bar([[data.et.B.stats.numSaccs_fullAve(1) data.et.B.stats.numSaccs_grpAve(1,:)];...
        [data.et.B.stats.numSaccs_fullAve(2) data.et.B.stats.numSaccs_grpAve(2,:)];
        [data.et.B.stats.numSaccs_fullAve(3) data.et.B.stats.numSaccs_grpAve(3,:)]]);
    hold on
    
    % Beeswarm
    % For average
    x_val = (hb(1).XData + hb(1).XOffset);
    bee_bin_width = .1;
    bee_spread_width = .2;
    beePlot = plotSpread({data.et.B.stats.numSaccs_subj(1,:),...
        data.et.B.stats.numSaccs_subj(2,:),...
        data.et.B.stats.numSaccs_subj(3,:)},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.5 .5 .5]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    % For groups
    for iI=1:numel(hb)-1
        x_val = (hb(iI+1).XData + hb(iI+1).XOffset);
        bee_bin_width = .1;
        bee_spread_width = .2;
        beePlot = plotSpread({data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{iI}),...
            data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{iI}),...
            data.et.B.stats.numSaccs_subj(3,data.et.B.grpIdx{iI})},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.5 .5 .5]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
    end
    
    % Errorbars
    ngroups = 3;
    nbars = 4;
    groupwidth = min(nbars/(nbars+1.5));
    holderMean =[[data.et.B.stats.numSaccs_fullAve(1) data.et.B.stats.numSaccs_grpAve(1,:)];...
        [data.et.B.stats.numSaccs_fullAve(2) data.et.B.stats.numSaccs_grpAve(2,:)];...
        [data.et.B.stats.numSaccs_fullAve(3) data.et.B.stats.numSaccs_grpAve(3,:)]];
    holderSte = [[data.et.B.stats.numSaccs_fullSte(1) data.et.B.stats.numSaccs_grpSte(1,:)];...
        [data.et.B.stats.numSaccs_fullSte(2) data.et.B.stats.numSaccs_grpSte(2,:)];...
        [data.et.B.stats.numSaccs_fullSte(3) data.et.B.stats.numSaccs_grpSte(3,:)]];
    for iI=1:nbars
        x = (1:ngroups) - groupwidth/2 + (2*iI-1) * groupwidth / (2*nbars);
        errorbar(x,holderMean(:,iI),holderSte(:,iI),'.k','LineWidth',2);
    end
    
    % Plot the significance
    text(2.9,500,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.numSaccs_kw_Table{2,3},')=',...
        data.et.B.stats.numSaccs_kw_Table{2,5},...
        'p=',data.et.B.stats.numSaccs_kw_Table{2,6}));
    
    set(gca,'XTick',1:3,'XTickLabel',{'Left Eye','Right Eye'},'fontsize',15)
    set(hb,'linewidth',2)
    hb(1).FaceColor = [0 0 0];
    hb(2).FaceColor = col_list{1};
    hb(3).FaceColor = col_list{2};
    hb(4).FaceColor = col_list{3};
    set(gca,'XTick',1:2,'fontsize',15)
    
    box off
    ylabel('Num Saccades','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Average Number of Saccades','fontsize',20)
    
end

%% Plot average number blinks for combined events only
if options.displayETFigs
    figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [6.1906 3.1432];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.75 ...
        (figSize.baseSize(3)*.75)*(figSize.aspectRatio(2)/figSize.aspectRatio(1))];   % Size/postion of fig
    set(gcf,'color','w')
    set(gcf,'Position', figSize.figSize)
    
    %% Blinks
    hb1 = bar(1,[data.et.B.stats.numBlinks_fullAve(3)]);
    hold on
    hb2 = bar(2,[data.et.B.stats.numBlinks_grpAve(3,1)]);
    hb3 = bar(3,[data.et.B.stats.numBlinks_grpAve(3,2)]);
    hb4 = bar(4,[data.et.B.stats.numBlinks_grpAve(3,3)]);
    
    % Beeswarm
    % For average
    x_val = 1:4;
    bee_bin_width = .1;
    bee_spread_width = .75;
    beePlot = plotSpread({data.et.B.stats.numBlinks_subj(3,:),...
        data.et.B.stats.numBlinks_subj(3,data.et.B.grpIdx{1}),...
        data.et.B.stats.numBlinks_subj(3,data.et.B.grpIdx{2}),...
        data.et.B.stats.numBlinks_subj(3,data.et.B.grpIdx{3})},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.5 .5 .5]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    
    % Errorbars
    errorbar([data.et.B.stats.numBlinks_fullAve(3) data.et.B.stats.numBlinks_grpAve(3,1) ...
        data.et.B.stats.numBlinks_grpAve(3,2) data.et.B.stats.numBlinks_grpAve(3,3)],...
        [data.et.B.stats.numBlinks_fullSte(3) data.et.B.stats.numBlinks_grpSte(3,1) ...
        data.et.B.stats.numBlinks_grpSte(3,2) data.et.B.stats.numBlinks_grpSte(3,3)],...
        '.k','LineWidth',2);
    
    % Plot the significance
    text(3.9,75,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.numBlinks_kw_Table{2,3},')=',...
        data.et.B.stats.numBlinks_kw_Table{2,5},...
        'p=',data.et.B.stats.numBlinks_kw_Table{2,6}));
    
    set(gca,'XTick',1:4,'XTickLabel',{'All','Controls','Relatives','Psychosis'},...
        'ylim',[0 80],'fontsize',15)
    set(hb1,'linewidth',2)
    set(hb2,'linewidth',2)
    set(hb3,'linewidth',2)
    set(hb4,'linewidth',2)
    hb1.FaceColor = [0 0 0];
    hb2.FaceColor = col_list{1};
    hb3.FaceColor = col_list{2};
    hb4.FaceColor = col_list{3};
    
    box off
    ylabel('Num Blinks','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])],'FontSize',8)
    title('Average Number of Blinks','fontsize',20)
    
end

%% Plot average number eye movements for combined events only
if options.displayETFigs
    figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [6.1906 3.1432];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.75 ...
        (figSize.baseSize(3)*.75)*(figSize.aspectRatio(2)/figSize.aspectRatio(1))];   % Size/postion of fig
    set(gcf,'color','w')
    set(gcf,'Position', figSize.figSize)
    %% Saccades
    hb1 = bar(1,[data.et.B.stats.numSaccs_fullAve(3)]);
    hold on
    hb2 = bar(2,[data.et.B.stats.numSaccs_grpAve(3,1)]);
    hb3 = bar(3,[data.et.B.stats.numSaccs_grpAve(3,2)]);
    hb4 = bar(4,[data.et.B.stats.numSaccs_grpAve(3,3)]);
    
    % Beeswarm
    % For average
    x_val = 1:4;
    bee_bin_width = .1;
    bee_spread_width = .75;
    beePlot = plotSpread({data.et.B.stats.numSaccs_subj(3,:),...
        data.et.B.stats.numSaccs_subj(3,data.et.B.grpIdx{1}),...
        data.et.B.stats.numSaccs_subj(3,data.et.B.grpIdx{2}),...
        data.et.B.stats.numSaccs_subj(3,data.et.B.grpIdx{3})},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.5 .5 .5]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    
    % Errorbars
    errorbar([data.et.B.stats.numSaccs_fullAve(3) data.et.B.stats.numSaccs_grpAve(3,1) ...
        data.et.B.stats.numSaccs_grpAve(3,2) data.et.B.stats.numSaccs_grpAve(3,3)],...
        [data.et.B.stats.numSaccs_fullSte(3) data.et.B.stats.numSaccs_grpSte(3,1) ...
        data.et.B.stats.numSaccs_grpSte(3,2) data.et.B.stats.numSaccs_grpSte(3,3)],...
        '.k','LineWidth',2);
    
    % Plot the significance
    text(3.9,235,...
        sprintf('%s%d%s%.3f\n%s%.3f','X²(',data.et.B.stats.numSaccs_kw_Table{2,3},')=',...
        data.et.B.stats.numSaccs_kw_Table{2,5},...
        'p=',data.et.B.stats.numSaccs_kw_Table{2,6}));
    
    set(gca,'XTick',1:4,'XTickLabel',{'All','Controls','Relatives','Psychosis'},...
        'ylim',[0 250],'fontsize',15)
    set(hb1,'linewidth',2)
    set(hb2,'linewidth',2)
    set(hb3,'linewidth',2)
    set(hb4,'linewidth',2)
    hb1.FaceColor = [0 0 0];
    hb2.FaceColor = col_list{1};
    hb3.FaceColor = col_list{2};
    hb4.FaceColor = col_list{3};
    
    box off
    ylabel('Num Saccades','fontsize',15)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Average Number of Saccades','fontsize',20)
    
end

%% Plot blink # and switch rate correlations
if options.displayETFigs
    figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.5 ...
        figSize.baseSize(4)];   % Size/postion of fig
    set(gcf,'color','w')
    set(gcf,'Position', figSize.figSize)
    
    % Left eye
    subplot(3,1,1)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{1}))))',...
        data.et.B.stats.numBlinks_subj(1,...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{1}))))',...
        'og','MarkerFaceColor','w','MarkerSize',8)
    hold on
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{2}))))',...
        data.et.B.stats.numBlinks_subj(1,...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{2}))))',...
        'ob','MarkerFaceColor','w','MarkerSize',8)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{3}))))',...
        data.et.B.stats.numBlinks_subj(1,...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numBlinks_subj(1,data.et.B.grpIdx{3}))))',...
        'or','MarkerFaceColor','w','MarkerSize',8)
    
    % Plot the correlation
    text(.25,400,...
        sprintf('%s%.3f\n%s%.3f','r=',data.et.B.stats.corr_switch_blinkLeft_r,...
        'p=',data.et.B.stats.corr_switch_blinkLeft_p));
        
    box off
    ylabel('Num Blinks','fontsize',15)
    xlabel('Switch Rate (Hz)')
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gca,'XScale','log','xtick',[0.01 0.025 0.05 0.1 0.25 0.5])
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Correlation Between Blink Number and Switch Rate (Left Eye)','fontsize',12)
    
    % Right eye
    subplot(3,1,2)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{1}))))',...
        data.et.B.stats.numBlinks_subj(2,...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{1}))))',...
        'og','MarkerFaceColor','w','MarkerSize',8)
    hold on
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{2}))))',...
        data.et.B.stats.numBlinks_subj(2,...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{2}))))',...
        'ob','MarkerFaceColor','w','MarkerSize',8)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{3}))))',...
        data.et.B.stats.numBlinks_subj(2,...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{3}))))',...
        'or','MarkerFaceColor','w','MarkerSize',8)
    
    % Plot the correlation
    text(.25,400,...
        sprintf('%s%.3f\n%s%.3f','r=',data.et.B.stats.corr_switch_blinkRight_r,...
        'p=',data.et.B.stats.corr_switch_blinkRight_p));
    
    box off
    ylabel('Num Blinks','fontsize',15)
    xlabel('Switch Rate (Hz)')
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gca,'XScale','log','xtick',[0.01 0.025 0.05 0.1 0.25 0.5])
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Correlation Between Blink Number and Switch Rate (Right Eye)','fontsize',12)
    
    
    % Combined eyes
    subplot(3,1,3)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{1}))))',...
        data.et.B.stats.numBlinks_subj(3,...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{1}))))',...
        'og','MarkerFaceColor','w','MarkerSize',8)
    hold on
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{2}))))',...
        data.et.B.stats.numBlinks_subj(3,...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{2}))))',...
        'ob','MarkerFaceColor','w','MarkerSize',8)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{3}))))',...
        data.et.B.stats.numBlinks_subj(3,...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numBlinks_subj(2,data.et.B.grpIdx{3}))))',...
        'or','MarkerFaceColor','w','MarkerSize',8)
    
    % Plot the correlation
    text(.25,100,...
        sprintf('%s%.3f\n%s%.3f','r=',data.et.B.stats.corr_switch_blinkBoth_r,...
        'p=',data.et.B.stats.corr_switch_blinkBoth_p));
    
    box off
    ylabel('Num Blinks','fontsize',12)
    xlabel('Switch Rate (Hz)')
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gca,'XScale','log','xtick',[0.01 0.025 0.05 0.1 0.25 0.5])
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Correlation Between Blink Number and Switch Rate (Comb. Eyes)','fontsize',12)
end

%% Plot sacc # and switch rate correlations
if options.displayETFigs
    figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.5 ...
        figSize.baseSize(4)];   % Size/postion of fig
    set(gcf,'color','w')
    set(gcf,'Position', figSize.figSize)
    
    % Left eye
    subplot(3,1,1)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{1}))))',...
        data.et.B.stats.numSaccs_subj(1,...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{1}))))',...
        'og','MarkerFaceColor','w','MarkerSize',8)
    hold on
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{2}))))',...
        data.et.B.stats.numSaccs_subj(1,...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{2}))))',...
        'ob','MarkerFaceColor','w','MarkerSize',8)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{3}))))',...
        data.et.B.stats.numSaccs_subj(1,...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numSaccs_subj(1,data.et.B.grpIdx{3}))))',...
        'or','MarkerFaceColor','w','MarkerSize',8)
    
    % Plot the correlation
    text(.25,400,...
        sprintf('%s%.3f\n%s%.3f','r=',data.et.B.stats.corr_switch_saccLeft_r,...
        'p=',data.et.B.stats.corr_switch_saccLeft_p));
        
    box off
    ylabel('Num Saccs','fontsize',15)
    xlabel('Switch Rate (Hz)')
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gca,'XScale','log','xtick',[0.01 0.025 0.05 0.1 0.25 0.5])
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Correlation Between Sacc Number and Switch Rate (Left Eye)','fontsize',12)
    
    % Right eye
    subplot(3,1,2)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{1}))))',...
        data.et.B.stats.numSaccs_subj(2,...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{1}))))',...
        'og','MarkerFaceColor','w','MarkerSize',8)
    hold on
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{2}))))',...
        data.et.B.stats.numSaccs_subj(2,...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{2}))))',...
        'ob','MarkerFaceColor','w','MarkerSize',8)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{3}))))',...
        data.et.B.stats.numSaccs_subj(2,...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{3}))))',...
        'or','MarkerFaceColor','w','MarkerSize',8)
    
    % Plot the correlation
    text(.25,400,...
        sprintf('%s%.3f\n%s%.3f','r=',data.et.B.stats.corr_switch_saccRight_r,...
        'p=',data.et.B.stats.corr_switch_saccRight_p));
    
    box off
    ylabel('Num Saccs','fontsize',15)
    xlabel('Switch Rate (Hz)')
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gca,'XScale','log','xtick',[0.01 0.025 0.05 0.1 0.25 0.5])
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Correlation Between Sacc Number and Switch Rate (Right Eye)','fontsize',12)
    
    
    % Combined eyes
    subplot(3,1,3)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{1}))))',...
        data.et.B.stats.numSaccs_subj(3,...
        data.et.B.grpIdx{1}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{1}))))',...
        'og','MarkerFaceColor','w','MarkerSize',8)
    hold on
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{2}))))',...
        data.et.B.stats.numSaccs_subj(3,...
        data.et.B.grpIdx{2}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{2}))))',...
        'ob','MarkerFaceColor','w','MarkerSize',8)
    plot(data.behav.B.Hz_flips_aveFull(...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{3}))))',...
        data.et.B.stats.numSaccs_subj(3,...
        data.et.B.grpIdx{3}(~isnan(data.et.B.stats.numSaccs_subj(2,data.et.B.grpIdx{3}))))',...
        'or','MarkerFaceColor','w','MarkerSize',8)
    
    % Plot the correlation
    text(.25,100,...
        sprintf('%s%.3f\n%s%.3f','r=',data.et.B.stats.corr_switch_saccBoth_r,...
        'p=',data.et.B.stats.corr_switch_saccBoth_p));
    
    box off
    ylabel('Num Saccs','fontsize',12)
    xlabel('Switch Rate (Hz)')
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    set(gca,'XScale','log','xtick',[0.01 0.025 0.05 0.1 0.25 0.5])
    set(gcf,'Position', figSize.figSize)
%     legend([sprintf('%s%d','All, n=',length(data.et.B.subjID)) ...
%         options.(['group_labels_', options.runType{iB}])])
    title('Correlation Between Sacc Number and Switch Rate (Comb. Eyes)','fontsize',12)
end


%% Average heat maps across participants
% Take an average of the heat maps
% First average across blocks (should only be 1 for A runs)
for iB=1:2   % A/B runs
    for iI=1:size(data.et.(options.runType{iB}).fixData_heatMap,1)   % Subjs
        for iZ=1:size(data.et.(options.runType{iB}).fixData_heatMap,2)   % B/Z day
            % Uncorrected heat maps
            if isempty(data.et.(options.runType{iB}).fixData_heatMap{iI,iZ})
                holderArray.(options.runType{iB}){iI,iZ} = [];
            else
                holderArray.(options.runType{iB}){iI,iZ} = nanmean(cat(3,data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{:,:}),3);
            end
            % Corrected heat maps
            if isempty(data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ})
                holderArrayCor.(options.runType{iB}){iI,iZ} = [];
            else
                holderArrayCor.(options.runType{iB}){iI,iZ} = nanmean(cat(3,data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{:,:}),3);
            end
        end
    end
end
% Average across group and B/Z visits
for iB=1:2
    for iI=1:size(data.et.(options.runType{iB}).grpIdx,2)
        % Average heatmap for each of the 3 groups
        data.et.(options.runType{iB}).stats.fullAveHeatMap(iI,:,:) =...
            nanmean(cat(3,holderArray.(options.runType{iB}){...
            data.et.(options.runType{iB}).grpIdx{iI},:}),3);
        % Average corrected heatmap for each of the 3 groups
        data.et.(options.runType{iB}).stats.fullCorAveHeatMap(iI,:,:) =...
            nanmean(cat(3,holderArrayCor.(options.runType{iB}){...
            data.et.(options.runType{iB}).grpIdx{iI},:}),3);
    end
end
clear holderArray holderArrayCor

% Plot averaged heat maps
if options.displayETFigs
    aveHeatMapRealFig = figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [11.2683 6.439];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.75 ...
        (figSize.baseSize(3)*.75)*(figSize.aspectRatio(2)/figSize.aspectRatio(1))];   % Size/postion of fig
    set(gcf,'color','w')
    % Real
    for iI=1:length(data.et.A.grpIdx)
        % Uncorrected
        subplot(2,3,iI)
        imagesc(squeeze(data.et.A.stats.fullAveHeatMap(iI,heatMapIdx(1):heatMapIdx(3),...
            heatMapIdx(2):heatMapIdx(4)))');
        hold on
        th = 0:pi/50:2*pi;
        xDegCirc = (options.PPD/2) * cos(th) + (heatMapSize(3)/2);
        yDegCirc = (options.PPD/2) * sin(th) + (heatMapSize(4)/2);
        plot(xDegCirc, yDegCirc,'w','linewidth',2);
        plot([heatMapSize(3)/2 heatMapSize(3)/2],...
            [heatMapSize(4)/2 heatMapSize(4)/2],'.r');
        set(gca,'YDir','normal',...
            'xtick',heatMapXTicks,...
            'xticklabels',heatMapXLabels,...
            'ytick',heatMapXTicks,...
            'yticklabels',heatMapXLabels)
        if iI==2   % Plot fixation heat map on middle title
            title(sprintf('%s\n%s','Fixation heat map - Real Switch Task',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        else
            title(sprintf('%s\n%s','',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        end
        
        % Corrected
        subplot(2,3,iI+3)
        imagesc(squeeze(data.et.A.stats.fullCorAveHeatMap(iI,heatMapIdx(1):heatMapIdx(3),...
            heatMapIdx(2):heatMapIdx(4)))');
        hold on
        th = 0:pi/50:2*pi;
        xDegCirc = (options.PPD/2) * cos(th) + (heatMapSize(3)/2);
        yDegCirc = (options.PPD/2) * sin(th) + (heatMapSize(4)/2);
        plot(xDegCirc, yDegCirc,'w','linewidth',2);
        plot([heatMapSize(3)/2 heatMapSize(3)/2],...
            [heatMapSize(4)/2 heatMapSize(4)/2],'.r');
        set(gca,'YDir','normal',...
            'xtick',heatMapXTicks,...
            'xticklabels',heatMapXLabels,...
            'ytick',heatMapXTicks,...
            'yticklabels',heatMapXLabels)
        if iI==2   % Plot fixation heat map on middle title
            title(sprintf('%s\n%s','Fixation heat map (Corrected) - Real Switch Task',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        else
            title(sprintf('%s\n%s','',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        end
    end
    box off
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    % Add figure y axis title
    han=axes(aveHeatMapRealFig,'visible','off');
    han.YLabel.Visible = 'on';
    ylabel(han,sprintf('%s\n','Degrees of Visual Angle'),'FontSize',15)
    
    % Bi-stable
    aveHeatMapBistabFig = figure();
    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aspectRatio = [11.2683 6.439];   % Aspect ratio
    figSize.figSize = [0 0 ...
        figSize.baseSize(3)*.75 ...
        (figSize.baseSize(3)*.75)*(figSize.aspectRatio(2)/figSize.aspectRatio(1))];   % Size/postion of fig
    set(gcf,'color','w')
    for iI=1:length(data.et.A.grpIdx)
        % Uncorrected
        subplot(2,3,iI)
        imagesc(squeeze(data.et.B.stats.fullAveHeatMap(iI,heatMapIdx(1):heatMapIdx(3),...
            heatMapIdx(2):heatMapIdx(4)))');
        hold on
        th = 0:pi/50:2*pi;
        xDegCirc = (options.PPD/2) * cos(th) + (heatMapSize(3)/2);
        yDegCirc = (options.PPD/2) * sin(th) + (heatMapSize(4)/2);
        plot(xDegCirc, yDegCirc,'w','linewidth',2);
        plot([heatMapSize(3)/2 heatMapSize(3)/2],...
            [heatMapSize(4)/2 heatMapSize(4)/2],'.r');
        set(gca,'YDir','normal',...
            'xtick',heatMapXTicks,...
            'xticklabels',heatMapXLabels,...
            'ytick',heatMapXTicks,...
            'yticklabels',heatMapXLabels)
        if iI==2   % Plot fixation heat map on middle title
            title(sprintf('%s\n%s','Fixation heat map - Bi-stable Switch Task',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        else
            title(sprintf('%s\n%s','',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        end
        
        % Corrected
        subplot(2,3,3+iI)
        imagesc(squeeze(data.et.B.stats.fullCorAveHeatMap(iI,heatMapIdx(1):heatMapIdx(3),...
            heatMapIdx(2):heatMapIdx(4)))');
        hold on
        th = 0:pi/50:2*pi;
        xDegCirc = (options.PPD/2) * cos(th) + (heatMapSize(3)/2);
        yDegCirc = (options.PPD/2) * sin(th) + (heatMapSize(4)/2);
        plot(xDegCirc, yDegCirc,'w','linewidth',2);
        plot([heatMapSize(3)/2 heatMapSize(3)/2],...
            [heatMapSize(4)/2 heatMapSize(4)/2],'.r');
        set(gca,'YDir','normal',...
            'xtick',heatMapXTicks,...
            'xticklabels',heatMapXLabels,...
            'ytick',heatMapXTicks,...
            'yticklabels',heatMapXLabels)
        if iI==2   % Plot fixation heat map on middle title
            title(sprintf('%s\n%s','Fixation heat map (Corrected) - Bi-stable Switch Task',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        else
            title(sprintf('%s\n%s','',...
                options.(['group_labels_', options.runType{iB}]){iI}),'fontsize',15)
        end
    end
    box off
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Position', figSize.figSize)
    % Add figure y axis title
    han=axes(aveHeatMapBistabFig,'visible','off');
    han.YLabel.Visible = 'on';
    ylabel(han,sprintf('%s\n','Degrees of Visual Angle'),'FontSize',15)
end


%% Run inverse correlation analysis 
% For all behavioral events, grab X seconds of data before and after
options.downSampleBinWidth = 100;
options.downSampleBinSpacing = 50;
options.corrWindow = -10000:options.downSampleBinSpacing:2500;   % Num timepoints to search back from behav event for eye events
options.eyeType = {'left','right'};
[data,options] = SFM_ET_invCorrAnalysis(data,options);

%% Run pupil size analysis
% For all behavioral events, grab X seconds of data before and after
options.downSampleBinSpacing = 100;
options.downSampleStart = 6000;
options.downSampleEnd = 5000;
options.timeWindow = -options.downSampleStart:options.downSampleBinSpacing:options.downSampleEnd;   % Num timepoints to search back from behav event for eye events
[data,options] = SFM_ET_pupilSizeAnalysis(data,options);

end




