% Group analysis script to run multiple participants through structure from motion
% control analysis.

function [dataAve, datafull] = SFM_Behav_Control_Group_Analysis(options)

addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode'));

% clear all; close all;
if ~exist('options','var') % parameters
    options = [];
end
if ~isfield(options,'excludeRedcap')
    options.excludeRedcap = 1; % exclude subjects based on control task performance, 0 = no, 1 = yes
    addpath(genpath('/home/shaw-raid1/matlab_tools/COP_analysis.git')) % Add path for redcap excludion function (in COP_analysis.git)
end
if ~isfield(options,'dateCutoff')
    options.dateCutoff = 0;   % Only take participants after a given date
    options.dateCutoffVal = 20210815;   % Don't take ppt after this date
    options.dateCutoffValConverted = datenum(num2str(options.dateCutoffVal),'yyyymmdd');
end
if ~isfield(options,'displayControlFigs')
    options.displayControlFigs = 0; % 1 = on, 0 = off
end
if ~isfield(options,'buttonEventTesting')
    options.buttonEventTesting = 0;
end
if ~isfield(options,'subj_group_def')   % Define what groups you want to compare
    % 1=controls x probands x relatives
    % 2=controls x SZ x BP
    % 3=SZ x SCA x BP
    % 4=controls x SCZ+SCA x BP
    options.subj_group_def = 1;
end
if ~isfield(options,'includeSubjWithNoRealSwitch')
    options.includeSubjWithNoRealSwitch = 1;   % Include participants w/out A runs
end
if ~isfield(options,'excludeBenzoUsers')
   options.excludeBenzoUsers = 0; 
end
if ~isfield(options,'plotAccuracy')
    options.plotAccuracy = 0;
end
if ~isfield(options,'plotRT')
    options.plotRT = 0;
end
if ~isfield(options,'plotAvePerDur_TowardAway')
    options.plotAvePerDur_TowardAway = 0;
end
if ~isfield(options,'plotDistTotalResp')
    options.plotDistTotalResp = 0;
end
if ~isfield(options,'plotDistRespMade')
    options.plotDistRespMade = 0;
end

%% Pull in the data
% Make the subject list out of the file names
options.curDur = pwd;
cd /home/shaw-raid1/data/pHCP.git/subjectResponseFiles/
if options.buttonEventTesting == 0
    fileNames = dir('P*SFM*A*.txt');
    fileNamesB = dir('P*SFM*B*.txt');
elseif options.buttonEventTesting == 1
    buttonEventTesting_fileName = 'KWK_ButtonEvent_Testing*';
    fileNames = dir(buttonEventTesting_fileName);
end
cd(options.curDur)

dataAve.plotPartFigs = 0;

dataAve.responseNumCutoff = 7;   % Response # cuttoff
dataAve.reactionTimeCutoff = 4;   % Reaction time cutoff (s)

options.partExcludeList = {''}; % Particpant P1012343 was being included, but looks fine? Removed from list - KWK 20231206

% Clean up and grab all the relevant SFM participants from 'fileNames'
for iI = 1:size(fileNames)
    dataHolder(iI).fileName = fileNames(iI).name;
    dataHolder(iI).partName = dataHolder(iI).fileName(1:8);
    dataHolder(iI).partNum = str2double(dataHolder(iI).fileName(2:8));
    if dataHolder(iI).partNum<2000000
        dataHolder(iI).partGroup = 1;
    elseif dataHolder(iI).partNum>=2000000 && dataHolder(iI).partNum<6000000
        dataHolder(iI).partGroup = 2;
    elseif dataHolder(iI).partNum>=6000000
        dataHolder(iI).partGroup = 3;
    end
    dataHolder(iI).date = dataHolder(iI).fileName(10:17);
    dataHolder(iI).time = dataHolder(iI).fileName(19:22);
    dataHolder(iI).expName = dataHolder(iI).fileName(24:26);
    dataHolder(iI).runNum = dataHolder(iI).fileName(37:38);
    dataHolder(iI).expTypeLabel = dataHolder(iI).fileName(32);
    if strcmp(dataHolder(iI).expTypeLabel,'A')
        dataHolder(iI).expType = 1;
    elseif strcmp(dataHolder(iI).expTypeLabel,'B')
        dataHolder(iI).expType = 2;
    end
    dataHolder(iI).responseNumCutoff = dataAve.responseNumCutoff;
    dataHolder(iI).reactionTimeCutoff = dataAve.reactionTimeCutoff;
end

for iI = 1:size(fileNamesB)
    dataHolderB(iI).fileName = fileNamesB(iI).name;
    dataHolderB(iI).partName = dataHolderB(iI).fileName(1:8);
    dataHolderB(iI).partNum = str2double(dataHolderB(iI).fileName(2:8));
    if dataHolderB(iI).partNum<2000000
        dataHolderB(iI).partGroup = 1;
    elseif dataHolderB(iI).partNum>=2000000 && dataHolderB(iI).partNum<6000000
        dataHolderB(iI).partGroup = 2;
    elseif dataHolderB(iI).partNum>=6000000
        dataHolderB(iI).partGroup = 3;
    end
    dataHolderB(iI).date = dataHolderB(iI).fileName(10:17);
    dataHolderB(iI).time = dataHolderB(iI).fileName(19:22);
    dataHolderB(iI).expName = dataHolderB(iI).fileName(24:26);
    dataHolderB(iI).runNum = dataHolderB(iI).fileName(37:38);
    dataHolderB(iI).expTypeLabel = dataHolderB(iI).fileName(32);
    if strcmp(dataHolderB(iI).expTypeLabel,'A')
        dataHolderB(iI).expType = 1;
    elseif strcmp(dataHolderB(iI).expTypeLabel,'B')
        dataHolderB(iI).expType = 2;
    end
    dataHolderB(iI).responseNumCutoff = dataAve.responseNumCutoff;
    dataHolderB(iI).reactionTimeCutoff = dataAve.reactionTimeCutoff;
end



%% Exclusions
%% Exclude any participants based on redcap exclusions
% For A runs
if options.excludeRedcap==1
    exclusion_opts = [];
    exclusion_opts.subj_number = [dataHolder(:).partNum]'; % numeric, not a string
    exclusion_opts.date_number = datenum({dataHolder(:).date}','yyyymmdd'); % per the matlab function datenum, NOT a string!!
    exclusion_opts.overwrite_exclusion_csv = 0;
    
    options.redcap_exclusion_output = redcap_exclusion(exclusion_opts); % this function lives in the COP_analysis.git repo
    options.redcap_exclusion_output.excluded_binary_nonans = ...
        options.redcap_exclusion_output.excluded_binary; % Make list w/ 0s replacing NaNs
    options.redcap_exclusion_output.excluded_binary_nonans(...
        isnan(options.redcap_exclusion_output.excluded_binary)) = 0; 
elseif options.excludeRedcap==0
    options.redcap_exclusion_output.excluded_binary_nonans = zeros([length(dataHolder) 1]);
end
% For B runs
if options.excludeRedcap==1
    exclusion_opts = [];
    exclusion_opts.subj_number = [dataHolderB(:).partNum]'; % numeric, not a string
    exclusion_opts.date_number = datenum({dataHolderB(:).date}','yyyymmdd'); % per the matlab function datenum, NOT a string!!
    exclusion_opts.overwrite_exclusion_csv = 0;
    
    options.redcap_exclusion_outputB = redcap_exclusion(exclusion_opts); % this function lives in the COP_analysis.git repo
    options.redcap_exclusion_outputB.excluded_binary_nonans = ...
        options.redcap_exclusion_outputB.excluded_binary; % Make list w/ 0s replacing NaNs
    options.redcap_exclusion_outputB.excluded_binary_nonans(...
        isnan(options.redcap_exclusion_outputB.excluded_binary)) = 0; 
elseif options.excludeRedcap==0
    options.redcap_exclusion_outputB.excluded_binary_nonans = zeros([length(dataHolderB) 1]);
end

%% Exclude predefined subjects
% For A runs
options.partExcludeIdx = zeros([length(dataHolder),1]);
for iI=1:length(options.partExcludeList)
    options.partExcludeIdx(strcmp({dataHolder.partName}, options.partExcludeList(iI)))=1;
end
% For B runs
options.partExcludeIdxB = zeros([length(dataHolderB),1]);
for iI=1:length(options.partExcludeList)
    options.partExcludeIdxB(strcmp({dataHolderB.partName}, options.partExcludeList(iI)))=1;
end

%% Only include participants run before the cutoff date
% For A runs
if options.dateCutoff == 1
    options.dateCutoffIdx = zeros([length(dataHolder),1]);
    options.dateCutoffIdx(options.dateCutoffVal<=str2double({dataHolder.date})) = 1;
elseif options.dateCutoff == 0
    options.dateCutoffIdx = zeros([length(dataHolder),1]);
end
% For B runs
if options.dateCutoff == 1
    options.dateCutoffIdxB = zeros([length(dataHolderB),1]);
    options.dateCutoffIdxB(options.dateCutoffVal<=str2double({dataHolderB.date})) = 1;
elseif options.dateCutoff == 0
    options.dateCutoffIdxB = zeros([length(dataHolderB),1]);
end

%% Exclude datasets run in the same participant on the same
% day. Only include the second run from each day for any
% participant. (Do include multiple runs on different days)
% For A runs
[~,repeatIndex] = unique([dataHolder.partNum],'stable'); % Find participants with multiple runs
repeatIndex = setdiff(1:numel(dataHolder),repeatIndex);
dateDiff = str2double({dataHolder(repeatIndex(:)-1).date})...
    - str2double({dataHolder(repeatIndex(:)).date}); % Find dates that are the same for duplicate runs
options.repeatCutoffIdx = zeros([length(dataHolder),1]);
options.repeatCutoffIdx(repeatIndex(dateDiff==0)-1) = 1; % Find any differences = 0 (run on the same day)
% For A runs
[~,repeatIndexB] = unique([dataHolderB.partNum],'stable'); % Find participants with multiple runs
repeatIndexB = setdiff(1:numel(dataHolderB),repeatIndexB);
dateDiffB = str2double({dataHolderB(repeatIndexB(:)-1).date})...
    - str2double({dataHolderB(repeatIndexB(:)).date}); % Find dates that are the same for duplicate runs
options.repeatCutoffIdxB = zeros([length(dataHolderB),1]);
options.repeatCutoffIdxB(repeatIndexB(dateDiffB==0)-1) = 1; % Find any differences = 0 (run on the same day)

%% Make full exclusion index list
options.full_exclusion_ind = options.partExcludeIdx |...
    options.redcap_exclusion_output.excluded_binary_nonans |...
    options.dateCutoffIdx |...
    options.repeatCutoffIdx;
options.full_exclusion_indB = options.partExcludeIdxB |...
    options.redcap_exclusion_outputB.excluded_binary_nonans |...
    options.dateCutoffIdxB |...
    options.repeatCutoffIdxB;

% Remove participants based on all exclusion
dataHolder(options.full_exclusion_ind==1) = [];
dataHolderB(options.full_exclusion_indB==1) = [];

for iI=1:length(dataHolder)
    datafull(iI) = SFM_Behav_Control_Analysis(dataHolder(iI));
end

%% If we want to exclude participants based on benzo use - KWK 20231121
% For A runs
if options.excludeBenzoUsers == 1
    medListOpt.partNum = [dataHolder.partNum];
    medListOpt.partName = {dataHolder.partName};
    
    benzoExclusionList = medListCheck(medListOpt);
    benzoExclusionPartName = medListOpt.partNum(benzoExclusionList == 1);
end
% For B runs
if options.excludeBenzoUsers == 1
    medListOptB.partNum = [dataHolderB.partNum];
    medListOptB.partName = {dataHolderB.partName};
    
    benzoExclusionListB = medListCheck(medListOptB);
    benzoExclusionPartNameB = medListOptB.partNum(benzoExclusionListB == 1);
end

% Create list of participants using benzos
if options.excludeBenzoUsers == 1
    dataAve.benzoExclusionList = unique([benzoExclusionPartName benzoExclusionPartNameB]);
end

%% Make list of participants who have:
% 1) 0 A and 1/2 B visits
dataAve.part_w_no_A = setdiff([dataHolderB.partNum]', [dataHolder.partNum]');

% 2) 1 A and 2 B visits
% For each participant check who has 1 A and 2 B
partHolderA = [dataHolder.partNum]';
partHolderB = [dataHolderB.partNum]';
uniquePartHolderA = unique(partHolderA);
uniquePartHolderB = unique(partHolderB);
counter = 0;
for iI = 1:length(uniquePartHolderB)
    if sum(uniquePartHolderB(iI) == partHolderB)==2 & sum(uniquePartHolderB(iI) == partHolderA)==1
        counter = counter+1;
        dataAve.part_w_1_A_2_B(counter) = uniquePartHolderB(iI);
    end
end


%% Import diagnostic/demographic info
% DON'T NEED ANYMORE B/C USING MPS SCRIPT 'run_subj_group_def' - KWK 20220205
% % Must implement group level segmenting (probands, relatives, and
% % controls).
% dataAve.target_file = '/home/shaw-raid1/data/7T/demographics/PHCP7TfMRIDemo.csv';
% addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode/'))   % Add path for MPS function read_in_demog_data
% demogDataHolder = read_in_demog_data(dataAve);
% for iI=1:length(datafull)
%     subjIdx = find((ismember(demogDataHolder.Record_ID,...
%         datafull(iI).fileName(1:8))));
%     if isempty(subjIdx)
%         datafull(iI).diag = NaN;
%         datafull(iI).gender = NaN;
%         datafull(iI).age = NaN;
%         datafull(iI).iq = NaN;
%         datafull(iI).education = NaN;
%         datafull(iI).dxId = NaN;
%     else
%         datafull(iI).diag =...
%             table2array(demogDataHolder.P_HCP_Preliminary_Dx(subjIdx));
%         datafull(iI).gender =...
%             table2array(demogDataHolder.Gender(subjIdx));
%         datafull(iI).age =...
%             demogDataHolder.Age(subjIdx);
%         datafull(iI).iq =...
%             demogDataHolder.Estimated_IQ(subjIdx);
%         datafull(iI).education =...
%             demogDataHolder.Education(subjIdx);
%         datafull(iI).dxId =...
%             demogDataHolder.Dx_code(subjIdx);
%     end
% end

%% Define the groups
% 1 = controls, relatives, probands; 2 = controls, SZ, BP
    % 3 = SZ, schizoaffective (SCA), BP; 4 = healthy (con+rel),
    % SZ+SCA, bipolar,
group_def_opt = [];
group_def_opt.subj_group_def = options.subj_group_def;
group_def_opt.subj_number = [datafull(:).partNum]';

group_def_out = run_subj_group_def( group_def_opt ); % mps 20220127 changing how we use subj group def

% Set group colors/idxing from group_def_out to remain consistent w/ other
% analysis code group defs. 
options.col_list{1} = group_def_out.use_colors_RGB{1};
options.col_list{2} = group_def_out.use_colors_RGB{2};
options.col_list{3} = group_def_out.use_colors_RGB{3};
grpIdx{1} = group_def_out.g1_idx;
grpIdx{2} = group_def_out.g2_idx;
grpIdx{3} = group_def_out.g3_idx;
grpLabel{1} = group_def_out.g1_label;
grpLabelShort{1} = group_def_out.g1_short;
grpLabel{2} = group_def_out.g2_label;
grpLabelShort{2} = group_def_out.g2_short;
grpLabel{3} = group_def_out.g3_label;
grpLabelShort{3} = group_def_out.g3_short;
legendTitle = grpLabel; 



%% Group analysis


for n=1:length(grpIdx)   
    
    %% Number of responses (corrected and uncorrected)
    % Make a histogram of number of response
    dataAve.numResponses{n} = [datafull(grpIdx{n}).responseBinNum];
    dataAve.numResponsesCorrected{n} = [datafull(grpIdx{n}).responseBinNumCorrected];   % NumResponses where the response is made w/in 4 secs of stim onset
    
    %% Average response time (diff between switch onset and first response)
    % Average the response time across participants, only
    % including subjects w/ total responses >= responseNumCutoff
    holderIdx = [datafull(:).responseBinNum]>=dataAve.responseNumCutoff;
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseDiffFull{n} = nanmean(vertcat(datafull(holderIdx).responseDiff),2);
    dataAve.responseDiff{n} = nanmean(vertcat(datafull(holderIdx).responseDiff));
    dataAve.responseDiffSTE{n} = nanstd(vertcat(datafull(holderIdx).responseDiff))/...
        sqrt(length(vertcat(datafull(holderIdx).responseDiff))-1);
    clear holderIdx
    
    % Make a distribution of all response times
    dataAve.responseDiffDist{n} = vertcat(datafull(grpIdx{n}).responseDiff);
    
    % Also looks at RT including all participants and only responses made within 4s of switch onset
    holderIdx = ones([1, length([datafull(:).responseBinNumCorrected])]);
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseDiffNoCuttoffCorr{n} = nanmean(vertcat(datafull(holderIdx).responseDiff),2);
    clear holderIdx
    
    % Also looks at RT including all participants (no response time cutoff)
    holderIdx = ones([1, length([datafull(:).responseBinNum])]);
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseDiffNoCuttoff{n} = nanmean(vertcat(datafull(holderIdx).responseDiff),2);
    clear holderIdx
    
    %% Average accuracy
    % Average the accuracy for all participants with >= responseNumCutoff
    holderIdx = [datafull(:).responseBinNum]>=dataAve.responseNumCutoff;
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseAccFull{n} = vertcat(datafull(holderIdx).responseAcc);
    dataAve.responseAcc{n} = nanmean(vertcat(datafull(holderIdx).responseAcc));
    dataAve.responseAccSTE{n} = nanstd(vertcat(datafull(holderIdx).responseAcc))/sqrt(length(datafull)-1);
    clear holderIdx
    
    % Also looks at RT including all participants (no response time cutoff)
    holderIdx = ones([1, length([datafull(:).responseBinNum])]);
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseAccNoCuttoff{n} = vertcat(datafull(holderIdx).responseAcc);   
    clear holderIdx
    
    % Also looks at Acc including participants w/ >= 7 responses and only responses made within 4s of switch onset
    holderIdx = [datafull(:).responseBinNumCorrected]>=dataAve.responseNumCutoff;
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseAccCorrFull{n} = vertcat(datafull(holderIdx).responseAccCorrected);
    dataAve.responseAccCorr{n} = nanmean(vertcat(datafull(holderIdx).responseAccCorrected));
    dataAve.responseAccCorrSTE{n} = nanstd(vertcat(datafull(holderIdx).responseAccCorrected))/sqrt(length(datafull)-1);
    clear holderIdx
    
    % Also looks at Acc including all participants (no response time cutoff)
    holderIdx = ones([1, length([datafull(:).responseBinNumCorrected])]);
    holderIdx = holderIdx & grpIdx{n}';
    dataAve.responseAccNoCuttoffCorr{n} = vertcat(datafull(holderIdx).responseAccCorrected);   
    clear holderIdx
    
    %% Look at the type of response
    % Look at number of responses away/towards the control direction
    % (need average duration of the percept).
    respDirCounter = 0;
    % Create a holder index variable to chose the correct participants
    % out of the list for this particular group.
    holderIdx = find(grpIdx{n});
    for iI=1:length(holderIdx)
        
        % First thing to do is remove any repeat presses within on time bin.
        responseTypesHolder = datafull(holderIdx(iI)).responseTypesCell;
        responseTimesHolder = datafull(holderIdx(iI)).responseTimesCell;
        for j=1:length(responseTypesHolder)
            if length(responseTypesHolder{j})>1
                counter = 1;
                while 1
                    if ~(counter==length(responseTypesHolder{j}))
                        % Compare the kth response to the next response in the
                        % time bin.
                        if responseTypesHolder{j}(counter) == responseTypesHolder{j}(counter+1)
                            responseTypesHolder{j}(counter+1) = [];
                            responseTimesHolder{j}(counter+1) = [];
                        else
                            counter = counter+1;
                        end
                    else
                        break
                    end
                end
                clear counter
            end
        end
        datafull(holderIdx(iI)).responseTimesNoRepeats = responseTimesHolder;
        datafull(holderIdx(iI)).responseTypesNoRepeats = responseTypesHolder;
        
        % Then define each response as away/towards the physical direction of
        % motion of the stimuli.
        for j=1:length(responseTypesHolder)
            for k=1:length(responseTypesHolder{j})
                if responseTypesHolder{j}(k) == datafull(holderIdx(iI)).controlSwitchTypes(j)
                    datafull(holderIdx(iI)).responseDirection{j}(k) = 1;
                else
                    datafull(holderIdx(iI)).responseDirection{j}(k) = 0;
                end
            end
        end
        clear responseTypeHolder responseTimesHolder
        
        if datafull(holderIdx(iI)).responseBinNum >= dataAve.responseNumCutoff
            % Next, record how many times each participant responded either away or
            % toward the physical motion direction.
            respDirCounter = respDirCounter+1;
            dataAve.responseDirection{n}(respDirCounter,1) = 0;
            dataAve.responseDirection{n}(respDirCounter,2) = 0;
            for j=1:length(datafull(holderIdx(iI)).responseDirection)
                dataAve.responseDirection{n}(respDirCounter,1) = dataAve.responseDirection{n}(respDirCounter,1) +...
                    sum(datafull(holderIdx(iI)).responseDirection{j});
                dataAve.responseDirection{n}(respDirCounter,2) = dataAve.responseDirection{n}(respDirCounter,2) +...
                    (length(datafull(holderIdx(iI)).responseDirection{j}) - sum(datafull(holderIdx(iI)).responseDirection{j}));
            end
            
            % Make a holder matrix of the types and times for each
            % participant.
            counter=1;
            for j=1:length(datafull(holderIdx(iI)).responseDirection)
                for k=1:length(datafull(holderIdx(iI)).responseDirection{j})
                    dataAve.responseDirectionCombined{n}{respDirCounter}(counter,1) = datafull(holderIdx(iI)).responseDirection{j}(k);   % Types (1=towards, 2=away)
                    dataAve.responseDirectionCombined{n}{respDirCounter}(counter,2) = datafull(holderIdx(iI)).responseTimesNoRepeats{j}(k);   % Times
                    counter=counter+1;
                end
            end
            clear counter
            
            % Next, record the durations between the last response and the current
            % response for each type (away/towrads).
            counter=0;
            while 1
                counter=counter+1;
                if ~(counter == length(dataAve.responseDirectionCombined{n}{respDirCounter}))
                    if dataAve.responseDirectionCombined{n}{respDirCounter}(counter,1)==1   % Towards physical motion dir
                        dataAve.responseDirectionDuration{n}{respDirCounter,1}(counter) =...
                            dataAve.responseDirectionCombined{n}{respDirCounter}(counter+1,2)-...
                            dataAve.responseDirectionCombined{n}{respDirCounter}(counter,2);
                        dataAve.responseDirectionType{n}{respDirCounter}(counter) = 1;
                    elseif dataAve.responseDirectionCombined{n}{respDirCounter}(counter,1)==0   % Away from physical motion dir
                        dataAve.responseDirectionDuration{n}{respDirCounter,2}(counter) =...
                            dataAve.responseDirectionCombined{n}{respDirCounter}(counter+1,2)-...
                            dataAve.responseDirectionCombined{n}{respDirCounter}(counter,2);
                    end
                else   % If the last value break
                    break
                end
            end
            clear counter
            if ~isempty(dataAve.responseDirectionDuration{n}{respDirCounter,1})
                dataAve.responseDirectionDurationAve{n}(respDirCounter,1) =...
                    nanmean(dataAve.responseDirectionDuration{n}{respDirCounter,1});
            else
                dataAve.responseDirectionDurationAve{n}(respDirCounter,1) = NaN;
            end
            if ~isempty(dataAve.responseDirectionDuration{n}{respDirCounter,2})
                dataAve.responseDirectionDurationAve{n}(respDirCounter,2) =...
                    nanmean(dataAve.responseDirectionDuration{n}{respDirCounter,2});
            else
                dataAve.responseDirectionDurationAve{n}(respDirCounter,2) = NaN;
            end
        end
    end
end


%% Plot participant data
if dataAve.plotPartFigs
    for iI=1:length(data)
        figure()
        % Accuracy
        subplot(2,1,1)
        bar(data(iI).responseAcc,'b');
        ylim([0 100])
        ylabel('Percent of Switches Detected')
        title(sprintf('%s%s',subjList{iI}(1:4),' Accuracy'))
        
        % RT
        subplot(2,1,2)
        bar(1:length(data(iI).responseDiff),data(iI).responseDiff,'b');
        hold on
        bar(length(data(iI).responseDiff)+1,nanmean(data(iI).responseDiff),'r');
        errorbar(length(data(iI).responseDiff)+1,nanmean(data(iI).responseDiff),nanstd(data(iI).responseDiff)/sqrt(length(data(iI).responseDiff)-1),'.k')
        ylim([-1 10])
        xticks(1:length(data(iI).responseDiff)+1)
        xticklabels([num2cell(1:length(data(iI).responseDiff)) 'Average'])
        ylabel('Time (s)')
        title(sprintf('%s%s',subjList{iI}(1:4),' Difference From Actual Switch Onset'))
    end
end

%% Calculate unique n for each group
for iI=1:length(grpIdx)
    holder = datafull(grpIdx{iI});
    for j=1:length(holder)
        holderList(j) = holder(j).partNum;
    end
    dataAve.numUniquePart(iI) = length(unique(holderList));
    clear holder holderList
end

%% Run stats
% If we're not looking at combined data...
if options.subj_group_def ~= 1
    if options.displayControlFigs
        show_stats_fig = 'on';
    else
        show_stats_fig = 'off';
    end
    
    % Make grouping array
    grpIdxArray = zeros([1,length(grpIdx{1})]);
    grpIdxArray(grpIdx{1}) = 1;
    grpIdxArray(grpIdx{2}) = 2;
    grpIdxArray(grpIdx{3}) = 3;
    
    % Remove 0's from the grpIdxArray, if not using ALL participants 
    % E.g. only looking at subgroups of patients vs controls
    grpIdxArray(grpIdxArray(:)==0) = [];
    
    % Run KW for accuracy test between groups, as groups are not normally
    % distributed/equal variance.
    [kruskall_wallis_3_groups.p,...
        kruskall_wallis_3_groups.table , ...
        kruskall_wallis_3_groups.stats] = ...
        kruskalwallis([dataAve.responseAccNoCuttoffCorr{1};...
        dataAve.responseAccNoCuttoffCorr{2};...
        dataAve.responseAccNoCuttoffCorr{3}],...
        grpIdxArray,...
        show_stats_fig);
    
    % Run KW for reaction time test between groups, as groups are not normally
    % distributed/equal variance.
    [kruskall_wallis_3_groups.p,...
        kruskall_wallis_3_groups.table , ...
        kruskall_wallis_3_groups.stats] = ...
        kruskalwallis([dataAve.responseDiffNoCuttoffCorr{1};...
        dataAve.responseDiffNoCuttoffCorr{2};...
        dataAve. responseDiffNoCuttoffCorr{3}],...
        grpIdxArray,...
        show_stats_fig);
end

%% Plot acc/rt
if options.plotAccuracy
    %% Plot fig for paper (acc including all subj (except < 4s on average), w/ line indicating low response cuttoff)
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 10;
    statsFontSize = 10;
    % Set figure size
    figSize.aveAcc.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.aveAcc.aspectRatio = [6.5 5];   % Aspect ratio
    figSize.aveAcc.figSize = [0 0 ...
        figSize.aveAcc.aspectRatio];   % Size/postion of fig
    
    % Plot as boxplot w/ beeswarm underneath
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    % Boxplot
    hb = boxplot([dataAve.responseAccNoCuttoffCorr{1};dataAve.responseAccNoCuttoffCorr{2};dataAve.responseAccNoCuttoffCorr{3}],...
        [zeros([length(dataAve.responseAccNoCuttoffCorr{1}),1])+1;...
        zeros([length(dataAve.responseAccNoCuttoffCorr{2}),1])+2;...
        zeros([length(dataAve.responseAccNoCuttoffCorr{3}),1])+3]);
    hold on
    % Beeswarm
    x_val = [1 2 3];
    bee_bin_width = .1;
    bee_spread_width = .5;
    beePlot = plotSpread({dataAve.responseAccNoCuttoffCorr{1},...
        dataAve.responseAccNoCuttoffCorr{2},...
        dataAve.responseAccNoCuttoffCorr{3}},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.8 .8 .8]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    hold on
    for iI=1:length(legendTitle)
        legendTitleHolder{iI} = sprintf('%s%s%d%s',legendTitle{iI},' (n=',dataAve.numUniquePart(iI),')');
    end
    set(gca,'XTick',1:3,'XTickLabel',legendTitleHolder,'fontsize',axisLabelFontSize)
    set(hb,'linewidth',2)
    hb2 = findobj(gca,'type','line');
    for iHB = 1:size(hb,2)
        set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
        set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
    end
    box off
    ylim([0 100])
    ylabel(sprintf('%s\n%s','Percent of','Switches Detected'),'fontsize',axisTitleFontSize,'color','k')
    title(sprintf('\n%s%d%s\n','Average Accuracy (RT < ',dataAve.reactionTimeCutoff,'s)'),'fontsize',titleFontSize)
    set(gca,'XColor','k','YColor','k')
    
    % Plot horiz line at cuttoff
    plot([0 4],[(dataAve.responseNumCutoff/11)*100 (dataAve.responseNumCutoff/11)*100],'--k');
    for iI=1:length(legendTitle)
        legendTitleHolder{iI} = sprintf('%s%s%d%s',legendTitle{iI},' (n=',dataAve.numUniquePart(iI),')');
    end
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.aveAcc.figSize,'color','w')
    
%     % Export file
%     % Resize figure
%     set(gcf,'paperunits','centimeters')
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'papersize',figSize.aveAcc.aspectRatio)
%     set(gcf,'paperposition',[0,0,figSize.aveAcc.aspectRatio(1),figSize.aveAcc.aspectRatio(2)])
%     set(gcf, 'renderer', 'painters');
%     print('PaperFig_AccFULL.eps','-depsc');
end
   
if options.displayControlFigs
    %% Plot Accuracy and Histogram of resp made < 4s
    figure('DefaultAxesPosition', [0.1, 0.1, .9, .9])
    figSize.paperFig.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.paperFig.aspectRatio = [5.5408 5.5711];   % Aspect ratio
    figSize.paperFig.figSize = [0 0 ...
        figSize.paperFig.baseSize(4)*...
        (figSize.paperFig.aspectRatio(1)/figSize.paperFig.aspectRatio(2))...
        figSize.paperFig.baseSize(4)];   % Size/postion of fig
    subplot(2,1,1)
    set(gca,'LineWidth',1)
    % Plot as boxplot w/ beeswarm underneath
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    % Boxplot
    hb = boxplot([dataAve.responseAccCorrFull{1};dataAve.responseAccCorrFull{2};dataAve.responseAccCorrFull{3}],...
        [zeros([length(dataAve.responseAccCorrFull{1}),1])+1;...
        zeros([length(dataAve.responseAccCorrFull{2}),1])+2;...
        zeros([length(dataAve.responseAccCorrFull{3}),1])+3]);
    hold on
    % Beeswarm
    x_val = [1 2 3];
    bee_bin_width = .1;
    bee_spread_width = .5;
    beePlot = plotSpread({dataAve.responseAccCorrFull{1},...
        dataAve.responseAccCorrFull{2},...
        dataAve.responseAccCorrFull{3}},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.8 .8 .8]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    hold on
    for iI=1:length(legendTitle)
        legendTitleHolder{iI} = sprintf('%s%s%d%s',legendTitle{iI},' (n=',dataAve.numUniquePart(iI),')');
    end
    set(gca,'XTick',1:3,'XTickLabel',legendTitleHolder,'fontsize',15)
    set(hb,'linewidth',2)
    hb2 = findobj(gca,'type','line');
    for iHB = 1:size(hb,2)
        set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
        set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
    end
    set(gca,'fontsize',15)
    set(gcf,'Position', figSize.paperFig.figSize)
    box off
    ylim([40 100]);
    ylabel(sprintf('%s\n%s','Percent of','Switches Detected'),'fontsize',15,'color','k')
    title(sprintf('\n%s%d%s\n','Average Accuracy (RT < ',dataAve.reactionTimeCutoff,'s)'),'fontsize',18)
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    
    subplot(2,1,2)
    edgesHolder = [1.5 2.5 3.5 4.5 5.5 6.6 7.5 8.5 9.5 10.5 11.5];
    h1 = histcounts(dataAve.numResponsesCorrected{1},edgesHolder);
    h2 = histcounts(dataAve.numResponsesCorrected{2},edgesHolder);
    h3 = histcounts(dataAve.numResponsesCorrected{3},edgesHolder);
    hBar = bar(edgesHolder(1:end-1),[h1;h2;h3]');
    hBar(1).FaceColor = options.col_list{1};
    hBar(2).FaceColor = options.col_list{2};
    hBar(3).FaceColor = options.col_list{3};
    hBar(1).EdgeColor = [0 0 0];
    hBar(2).EdgeColor = [0 0 0];
    hBar(3).EdgeColor = [0 0 0];
    hBar(1).LineWidth = 1;
    hBar(2).LineWidth = 1;
    hBar(3).LineWidth = 1;
    hold on
    plot([dataAve.responseNumCutoff-1 dataAve.responseNumCutoff-1],[0 100],'--k');
    for iI=1:length(legendTitle)
        legendTitleHolder{iI} = sprintf('%s%s%d%s',legendTitle{iI},' (n=',dataAve.numUniquePart(iI),')');
    end
    set(gca,'fontsize',15)
    legend({legendTitleHolder{:}},'FontName','Arial','FontSize',12,'LineWidth',[1],'EdgeColor',[0 0 0],'Location','northwest')
    ylabel(sprintf('%s\n%s','Number of','Participants'),'FontName','Arial','FontSize',15)
    xlabel('Number of Responses Made','FontName','Arial','FontSize',15)
    ylim([0 30]);
    set(gca,'xtick',[.5 1.5 2.5 3.5 4.5 5.5 6.6 7.5 8.5 9.5 10.5],'xticklabels',1:11)
    title(sprintf('%s%d%s','Histogram of Participant Responses Made Before ',dataAve.reactionTimeCutoff,'s'),'FontName','Arial','FontSize',18)
    set(gcf,'Position', figSize.paperFig.figSize)
    box off
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    
    % Export file
    % Resize figure
    set(gcf,'paperunits','centimeters')
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'papersize',figSize.paperFig.aspectRatio)
    set(gcf,'paperposition',[0,0,figSize.paperFig.aspectRatio(1),figSize.paperFig.aspectRatio(2)])
    set(gcf, 'renderer', 'painters');
    print('PaperFig_AccHistOfResps.eps','-depsc');
end

if options.plotDistRespMade
    %% Historgram of participant responses
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 10;
    statsFontSize = 10;
    % Set figure size
    figSize.histPartResp.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.histPartResp.aspectRatio = [6.5 5];   % Aspect ratio
    figSize.histPartResp.figSize = [0 0 ...
        figSize.histPartResp.aspectRatio];   % Size/postion of fig
    
    subplot(2,1,1)
    set(gca,'LineWidth',1)
    edgesHolder = [1.5 2.5 3.5 4.5 5.5 6.6 7.5 8.5 9.5 10.5 11.5];
    h1 = histcounts(dataAve.numResponses{1},edgesHolder);
    h2 = histcounts(dataAve.numResponses{2},edgesHolder);
    h3 = histcounts(dataAve.numResponses{3},edgesHolder);
    hBar = bar(edgesHolder(1:end-1),[h1;h2;h3]');
    hBar(1).FaceColor = options.col_list{1};
    hBar(2).FaceColor = options.col_list{2};
    hBar(3).FaceColor = options.col_list{3};
    hBar(1).EdgeColor = [0 0 0];
    hBar(2).EdgeColor = [0 0 0];
    hBar(3).EdgeColor = [0 0 0];
    hBar(1).LineWidth = 1;
    hBar(2).LineWidth = 1;
    hBar(3).LineWidth = 1;
    for iI=1:length(legendTitle)
        legendTitleHolder{iI} = sprintf('%s%s%d%s',legendTitle{iI},' (n=',dataAve.numUniquePart(iI),')');
    end
    legend({legendTitleHolder{:}},'FontName','Arial','FontSize',axisLabelFontSize,'LineWidth',[1],'EdgeColor',[0 0 0],'Location','northwest')
    ylabel(sprintf('%s\n%s','Number of','Participants'),'FontName','Arial','FontSize',axisTitleFontSize)
%     xlabel('Number of Responses Made','FontName','Arial','FontSize',axisTitleFontSize)
    set(gca,'xtick',[.5 1.5 2.5 3.5 4.5 5.5 6.6 7.5 8.5 9.5 10.5],'xticklabels',1:11)
    title('Histogram of Participant Responses','FontName','Arial','FontSize',titleFontSize)
    box off
    ylim([0 30]);
    set(gca,'XColor','k','YColor','k')
    subplot(2,1,2)
    h1 = histcounts(dataAve.numResponsesCorrected{1},edgesHolder);
    h2 = histcounts(dataAve.numResponsesCorrected{2},edgesHolder);
    h3 = histcounts(dataAve.numResponsesCorrected{3},edgesHolder);
    hBar = bar(edgesHolder(1:end-1),[h1;h2;h3]');
    hBar(1).FaceColor = options.col_list{1};
    hBar(2).FaceColor = options.col_list{2};
    hBar(3).FaceColor = options.col_list{3};
    hBar(1).EdgeColor = [0 0 0];
    hBar(2).EdgeColor = [0 0 0];
    hBar(3).EdgeColor = [0 0 0];
    hBar(1).LineWidth = 1;
    hBar(2).LineWidth = 1;
    hBar(3).LineWidth = 1;
    hold on
    plot([dataAve.responseNumCutoff-1 dataAve.responseNumCutoff-1],[0 100],'--k');
    ylabel(sprintf('%s\n%s','Number of','Participants'),'FontName','Arial','FontSize',axisTitleFontSize)
    xlabel('Number of Responses Made','FontName','Arial','FontSize',axisTitleFontSize)
    ylim([0 30]);
    set(gca,'xtick',[.5 1.5 2.5 3.5 4.5 5.5 6.6 7.5 8.5 9.5 10.5],'xticklabels',1:11)
    title(sprintf('%s%d%s','Histogram of Participant Responses Made Before ',dataAve.reactionTimeCutoff,'s'),'FontName','Arial','FontSize',titleFontSize)
    box off
    set(gca,'XColor','k','YColor','k')
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.histPartResp.figSize,'color','w')
end

if options.displayControlFigs
    %% Accuracy
%     %         figure()
%     %         bar([dataAve.responseAcc{1:length(legendTitle)}; dataAve.responseAccCorr{1:length(legendTitle)}]);
%     %         nBars = size(dataAve.responseAcc,2);
%     %         nGroups = size([dataAve.responseAcc{1} dataAve.responseAccCorr{1}],2);
%     %         groupWidth = min(.8, nBars/(nBars+1.5));
%     %         for iI=1:nBars
%     %             x = (1:nGroups) - groupWidth/2 + (2*iI-1) * groupWidth / (2*nBars);
%     %             errorbar(x,[dataAve.responseAcc{iI} dataAve.responseAccCorr{iI}],...
%     %                 [dataAve.responseAccSTE{iI} dataAve.responseAccCorrSTE{iI}],'.k');
%     %         end
%     % Change to box plot instead of bar graph - KWK 20200507
%     figure('DefaultAxesPosition', [0.1, 0.1, .9, .85])
%     figSize.aveAcc.baseSize = get(0,'Screensize');   % Base size in pixels
%     figSize.aveAcc.aspectRatio = [5.4546 5.7077];   % Aspect ratio
%     figSize.aveAcc.figSize = [0 0 ...
%         figSize.aveAcc.baseSize(4)*...
%         (figSize.aveAcc.aspectRatio(1)/figSize.aveAcc.aspectRatio(2))...
%         figSize.aveAcc.baseSize(4)];   % Size/postion of fig
%     % Plot as boxplot w/ beeswarm underneath
%     addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
%     % Boxplot
%     hb = boxplot([dataAve.responseAccCorrFull{1};dataAve.responseAccCorrFull{2};dataAve.responseAccCorrFull{3}],...
%         [zeros([length(dataAve.responseAccCorrFull{1}),1])+1;...
%         zeros([length(dataAve.responseAccCorrFull{2}),1])+2;...
%         zeros([length(dataAve.responseAccCorrFull{3}),1])+3]);
%     hold on
%     % Beeswarm
%     x_val = [1 2 3];
%     bee_bin_width = .1;
%     bee_spread_width = .5;
%     beePlot = plotSpread({dataAve.responseAccCorrFull{1},...
%         dataAve.responseAccCorrFull{2},...
%         dataAve.responseAccCorrFull{3}},...
%         'binWidth', bee_bin_width,...
%         'distributionColors', {[.8 .8 .8]},...
%         'xValues', x_val,...
%         'spreadWidth', bee_spread_width);
%     set(beePlot{1},'MarkerSize',10)
%     hold on
%     for iI=1:length(legendTitle)
%         legendTitleHolder{iI} = sprintf('%s%s%d%s',legendTitle{iI},' (n=',dataAve.numUniquePart(iI),')');
%     end
%     set(gca,'XTick',1:3,'XTickLabel',legendTitleHolder,'fontsize',15)
%     set(hb,'linewidth',2)
%     hb2 = findobj(gca,'type','line');
%     for iHB = 1:size(hb,2)
%         set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
%         set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
%     end
%     set(gcf,'Position', figSize.aveAcc.figSize)
%     box off
%     ylim([0 100])
%     ylabel(sprintf('%s\n%s','Percent of','Switches Detected'),'fontsize',15,'color','k')
%     title(sprintf('\n%s%d%s\n','Average Accuracy (RT < ',dataAve.reactionTimeCutoff,'s)'),'fontsize',18)
%     set(gcf,'color','w')
%     set(gca,'XColor','k','YColor','k')
end

if options.plotRT
%% RT averages box plot
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 10;
    statsFontSize = 10;
    % Set figure size
    figSize.rtAve.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.rtAve.aspectRatio = [6.5 5];   % Aspect ratio
    figSize.rtAve.figSize = [0 0 ...
        figSize.rtAve.aspectRatio];   % Size/postion of fig
    
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    % Boxplot
    hb = boxplot([dataAve.responseDiffNoCuttoffCorr{1};dataAve.responseDiffNoCuttoffCorr{2};dataAve.responseDiffNoCuttoffCorr{3}],...
        [zeros([length(dataAve.responseDiffNoCuttoffCorr{1}),1])+1;...
        zeros([length(dataAve.responseDiffNoCuttoffCorr{2}),1])+2;...
        zeros([length(dataAve.responseDiffNoCuttoffCorr{3}),1])+3]);
    hold on
    % Beeswarm
    x_val = [1 2 3];
    bee_bin_width = .1;
    bee_spread_width = .5;
    beePlot = plotSpread({dataAve.responseDiffNoCuttoffCorr{1};...
        dataAve.responseDiffNoCuttoffCorr{2};...
        dataAve.responseDiffNoCuttoffCorr{3}},...
        'binWidth', bee_bin_width,...
        'distributionColors', {[.8 .8 .8]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',10)
    hold on
    set(gca,'XTick',1:3,'XTickLabel',legendTitleHolder,'fontsize',axisLabelFontSize)
    set(hb,'linewidth',2)
    hb2 = findobj(gca,'type','line');
    for iHB = 1:size(hb,2)
        set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
        set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
    end
    box off
    ylim([-1 10])
    ylabel('Time (s)','fontsize',axisTitleFontSize)
    title(sprintf('\n%s\n%s\n','Average Difference From Actual','Switch Onset To First Correct Response (RT)'),...
        'fontsize',titleFontSize)
    set(gca,'XColor','k','YColor','k')
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.rtAve.figSize,'color','w')
    
    
    %% RT Full set bar graph
    %         figure()
    %         bar(1:length(dataAve.responseDiff{1}),...
    %             [dataAve.responseDiff{1};dataAve.responseDiff{2};dataAve.responseDiff{3}]');
    %         hold on
    %         aveHolder = zeros([length(dataAve.responseDiff{1}) length(dataAve.responseDiff)]);
    %         aveHolder = [dataAve.responseDiff{1}' dataAve.responseDiff{2}' dataAve.responseDiff{3}'];
    %         aveHolder(12,:) = [nanmean(dataAve.responseDiff{1}) nanmean(dataAve.responseDiff{2}) nanmean(dataAve.responseDiff{3})];
    %         bar(1:length(aveHolder),aveHolder)
    %         steHolder = zeros([length(dataAve.responseDiffSTE{1}) length(dataAve.responseDiffSTE)]);
    %         steHolder = [dataAve.responseDiffSTE{1}' dataAve.responseDiffSTE{2}' dataAve.responseDiffSTE{3}'];
    %         steHolder(12,:) = [nanstd(dataAve.responseDiff{1})/sqrt(length(dataAve.responseDiff{1}-1))...
    %             nanstd(dataAve.responseDiff{2})/sqrt(length(dataAve.responseDiff{2}-1))...
    %             nanstd(dataAve.responseDiff{3})/sqrt(length(dataAve.responseDiff{3}-1))];
    %         nBars = size(aveHolder,2);
    %         nGroups = size(aveHolder,1);
    %         groupWidth = min(.8, nBars/(nBars+1.5));
    %         for iI=1:nBars
    %             x = (1:nGroups) - groupWidth/2 + (2*iI-1) * groupWidth / (2*nBars);
    %             errorbar(x,aveHolder(:,iI),steHolder(:,iI),'.k')'Probands';
    %         end
    %         ylim([-1 10])
    %         set(gca,'xtick',1:length(dataAve.responseDiff{1})+1,...
    %             'xticklabels',[num2cell(1:length(dataAve.responseDiff{1})) 'Average'])
    %         ylabel('Time (s)')
    %         title(sprintf('%s\n%s','Average Difference From Actual','Switch Onset To First Correct Response (RT)'))
end

if options.plotDistTotalResp
    %% Distribution of RTs
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 10;
    statsFontSize = 10;
    % Set figure size
    figSize.avePecpDur.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.avePecpDur.aspectRatio = [6.5 5];   % Aspect ratio
    figSize.avePecpDur.figSize = [0 0 ...
        figSize.rtAve.aspectRatio];   % Size/postion of fig
    
%     sTitle = suptitle(sprintf('%s\n','Histogram of Reaction Times by Time Bin'));
%     sTitle.FontSize = titleFontSize+1;
    subplot(3,1,1)
    hist([dataAve.responseDiffDist{1}(:,1),dataAve.responseDiffDist{1}(:,2),dataAve.responseDiffDist{1}(:,3),...
        dataAve.responseDiffDist{1}(:,4),dataAve.responseDiffDist{1}(:,5),dataAve.responseDiffDist{1}(:,6),...
        dataAve.responseDiffDist{1}(:,7),dataAve.responseDiffDist{1}(:,8),dataAve.responseDiffDist{1}(:,9),...
        dataAve.responseDiffDist{1}(:,10),dataAve.responseDiffDist{1}(:,11)],[0:11])
    hold on
    xlim([-1 12])
    % legend({'Time Bin 1','Time Bin 2','Time Bin 3','Time Bin 4','Time Bin 5','Time Bin 6',...
    %     'Time Bin 7','Time Bin 8','Time Bin 9','Time Bin 10','Time Bin 11'})
    title(legendTitleHolder{1},'color',options.col_list{1},'fontsize',titleFontSize)
    set(gca,'XColor','k','YColor','k','fontsize',axisLabelFontSize)
    box off
    subplot(3,1,2)
    hist([dataAve.responseDiffDist{2}(:,1),dataAve.responseDiffDist{2}(:,2),dataAve.responseDiffDist{2}(:,3),...
        dataAve.responseDiffDist{2}(:,4),dataAve.responseDiffDist{2}(:,5),dataAve.responseDiffDist{2}(:,6),...
        dataAve.responseDiffDist{2}(:,7),dataAve.responseDiffDist{2}(:,8),dataAve.responseDiffDist{2}(:,9),...
        dataAve.responseDiffDist{2}(:,10),dataAve.responseDiffDist{2}(:,11)],[0:11])
    hold on
    xlim([-1 12])
    % legend({'Time Bin 1','Time Bin 2','Time Bin 3','Time Bin 4','Time Bin 5','Time Bin 6',...
    %     'Time Bin 7','Time Bin 8','Time Bin 9','Time Bin 10','Time Bin 11'})
    title(legendTitleHolder{2},'color',options.col_list{2},'fontsize',titleFontSize)
    ylabel('Number of Time Bins','fontsize',axisTitleFontSize)
    set(gca,'XColor','k','YColor','k','fontsize',axisLabelFontSize)
    box off
    subplot(3,1,3)
    hist([dataAve.responseDiffDist{3}(:,1),dataAve.responseDiffDist{3}(:,2),dataAve.responseDiffDist{3}(:,3),...
        dataAve.responseDiffDist{3}(:,4),dataAve.responseDiffDist{3}(:,5),dataAve.responseDiffDist{3}(:,6),...
        dataAve.responseDiffDist{3}(:,7),dataAve.responseDiffDist{3}(:,8),dataAve.responseDiffDist{3}(:,9),...
        dataAve.responseDiffDist{3}(:,10),dataAve.responseDiffDist{3}(:,11)],[0:11])
    hold on
    xlim([-1 12])
    % legend({'Time Bin 1','Time Bin 2','Time Bin 3','Time Bin 4','Time Bin 5','Time Bin 6',...
    %     'Time Bin 7','Time Bin 8','Time Bin 9','Time Bin 10','Time Bin 11'})
    title(legendTitleHolder{3},'color',options.col_list{3},'fontsize',titleFontSize)
    xlabel('Time (s)','fontsize',axisLabelFontSize)
    set(gca,'XColor','k','YColor','k','fontsize',axisLabelFontSize)
    box off

    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.avePecpDur.figSize,'color','w')
end
   
if options.displayControlFigs
    %% Histogram of response types (towards vs away)
%     figure()
%     figSize.histRespType.baseSize = get(0,'Screensize');   % Base size in pixels
%     figSize.histRespType.aspectRatio = [11.2982 5.4653];   % Aspect ratio
%     figSize.histRespType.figSize = [0 0 ...
%         figSize.histRespType.baseSize(4)*...
%         (figSize.histRespType.aspectRatio(1)/figSize.histRespType.aspectRatio(2))...
%         figSize.histRespType.baseSize(4)];   % Size/postion of fig
%     suptitle(sprintf('%s\n%s\n','Histogram of Responses Made Towards or','Away From Actual Rotation Direction'))
%     subplot(3,1,1)
%     hist([dataAve.responseDirection{1}(:,1),dataAve.responseDirection{1}(:,2)],[0:20])
%     hb = findobj(gca,'Type','patch');
%     hb(1).FaceColor = [1 1 0];
%     hb(1).EdgeColor = [0 0 0];
%     hb(2).FaceColor = [0 1 1];
%     hb(2).EdgeColor = [0 0 0];
%     legend({'Towards','Away'})
%     title(legendTitleHolder{1},...
%         'color',options.col_list{1},'fontsize',18)
%     xlim([-1 21])
%     set(gcf,'color','w')
%     set(gca,'XColor','k','YColor','k')
%     box off
%     subplot(3,1,2)
%     hist([dataAve.responseDirection{2}(:,1),dataAve.responseDirection{2}(:,2)],[0:20])
%     hb = findobj(gca,'Type','patch');
%     hb(1).FaceColor = [1 1 0];
%     hb(1).EdgeColor = [0 0 0];
%     hb(2).FaceColor = [0 1 1];
%     hb(2).EdgeColor = [0 0 0];
%     title(legendTitleHolder{2},...
%         'color',options.col_list{2},'fontsize',18)
%     ylabel('Number of Participants','fontsize',15)
%     xlim([-1 21])
%     set(gcf,'color','w')
%     set(gca,'XColor','k','YColor','k')
%     box off
%     subplot(3,1,3)
%     hist([dataAve.responseDirection{3}(:,1),dataAve.responseDirection{3}(:,2)],[0:20])
%     hb = findobj(gca,'Type','patch');
%     hb(1).FaceColor = [1 1 0];
%     hb(1).EdgeColor = [0 0 0];
%     hb(2).FaceColor = [0 1 1];
%     hb(2).EdgeColor = [0 0 0];
%     title(legendTitleHolder{3},...
%         'color',options.col_list{3},'fontsize',18)
%     xlabel('Number of Responses','fontsize',15)
%     xlim([-1 21])
%     set(gcf,'color','w')
%     set(gca,'XColor','k','YColor','k')
%     box off
%     set(gcf,'Position', figSize.histRespType.figSize)
end

if options.plotAvePerDur_TowardAway    
    %% Average perception durations for both response types
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 10;
    statsFontSize = 10;
    % Set figure size
    figSize.avePecpDur.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.avePecpDur.aspectRatio = [6.5 5];   % Aspect ratio
    figSize.avePecpDur.figSize = [0 0 ...
        figSize.rtAve.aspectRatio];   % Size/postion of fig
    
%     sTitle = suptitle(sprintf('%s\n','Average Duration Away and Towards the Physical Direction'));
%     sTitle.FontSize = 15;
    subplot(3,1,1)
    bar([nanmean(dataAve.responseDirectionDurationAve{1}(:,1)) nanmean(dataAve.responseDirectionDurationAve{1}(:,2))],...
        'FaceColor',options.col_list{1})
    hold on
    errorbar([nanmean(dataAve.responseDirectionDurationAve{1}(:,1)) nanmean(dataAve.responseDirectionDurationAve{1}(:,2))],...
        [nanstd(dataAve.responseDirectionDurationAve{1}(:,1))/sqrt(sum(~isnan(dataAve.responseDirectionDurationAve{1}(:,1)))-1)...
        nanstd(dataAve.responseDirectionDurationAve{1}(:,2))/sqrt(sum(~isnan(dataAve.responseDirectionDurationAve{1}(:,2)))-1)],'.k')
    title(legendTitleHolder{1},'fontsize',titleFontSize)
%     ylabel('Time (s)','fontsize',axisTitleFontSize)
    set(gca,'xticklabels',...
        {sprintf('%s%d%s','Towards (n=',sum(~isnan(dataAve.responseDirectionDurationAve{1}(:,1))),')')...
        sprintf('%s%d%s','Away (n=',sum(~isnan(dataAve.responseDirectionDurationAve{1}(:,2))),')')},...
        'fontsize',axisLabelFontSize);
    set(gca,'XColor','k','YColor','k')
    box off
    subplot(3,1,2)
    bar([nanmean(dataAve.responseDirectionDurationAve{2}(:,1)) nanmean(dataAve.responseDirectionDurationAve{2}(:,2))],...
        'FaceColor',options.col_list{2})
    hold on
    errorbar([nanmean(dataAve.responseDirectionDurationAve{2}(:,1)) nanmean(dataAve.responseDirectionDurationAve{2}(:,2))],...
        [nanstd(dataAve.responseDirectionDurationAve{2}(:,1))/sqrt(sum(~isnan(dataAve.responseDirectionDurationAve{2}(:,1)))-1)...
        nanstd(dataAve.responseDirectionDurationAve{2}(:,2))/sqrt(sum(~isnan(dataAve.responseDirectionDurationAve{2}(:,2)))-1)],'.k')
    title(legendTitleHolder{2},'fontsize',titleFontSize)
    ylabel('Time (s)','fontsize',axisTitleFontSize)
    set(gca,'xticklabels',...
        {sprintf('%s%d%s','Towards (n=',sum(~isnan(dataAve.responseDirectionDurationAve{2}(:,1))),')')...
        sprintf('%s%d%s','Away (n=',sum(~isnan(dataAve.responseDirectionDurationAve{2}(:,2))),')')},...
        'fontsize',axisLabelFontSize);
    set(gcf,'color','w')
    set(gca,'XColor','k','YColor','k')
    box off
    subplot(3,1,3)
    bar([nanmean(dataAve.responseDirectionDurationAve{3}(:,1)) nanmean(dataAve.responseDirectionDurationAve{3}(:,2))],...
        'FaceColor',options.col_list{3})
    hold on
    errorbar([nanmean(dataAve.responseDirectionDurationAve{3}(:,1)) nanmean(dataAve.responseDirectionDurationAve{3}(:,2))],...
        [nanstd(dataAve.responseDirectionDurationAve{3}(:,1))/sqrt(sum(~isnan(dataAve.responseDirectionDurationAve{3}(:,1)))-1)...
        nanstd(dataAve.responseDirectionDurationAve{3}(:,2))/sqrt(sum(~isnan(dataAve.responseDirectionDurationAve{3}(:,2)))-1)],'.k')
    title(legendTitleHolder{3},'fontsize',axisLabelFontSize)
%     ylabel('Time (s)','fontsize',axisTitleFontSize)
    set(gca,'xticklabels',...'Probands'
        {sprintf('%s%d%s','Towards (n=',sum(~isnan(dataAve.responseDirectionDurationAve{3}(:,1))),')')...
        sprintf('%s%d%s','Away (n=',sum(~isnan(dataAve.responseDirectionDurationAve{3}(:,2))),')')},...
        'fontsize',axisLabelFontSize);
    set(gca,'XColor','k','YColor','k')
    box off
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.avePecpDur.figSize,'color','w')
end

%% Save/Clean up
% Create demographics table

% Create a subfield to use as a filter matrix for the illusory analysis.
for iI=1:length(datafull)
    dataAve.illusoryCutoff(iI,1) = str2double(datafull(iI).partName(2:end));
    dataAve.illusoryCutoff(iI,2) = str2double(datafull(iI).date);
    dataAve.illusoryCutoff(iI,3) = str2double(datafull(iI).time);
    if datafull(iI).responseBinNumCorrected >= dataAve.responseNumCutoff
        dataAve.illusoryCutoff(iI,4) = 1;
    else
        dataAve.illusoryCutoff(iI,4) = 0;
    end
end

dataAve.subj_number = dataAve.illusoryCutoff(:,1);
for iI=1:length(datafull)
    dataAve.date_number(iI,1) = datenum(num2str(dataAve.illusoryCutoff(iI,2)),...
        'yyyymmdd');
end

% Make illusoryCutoff that includes participants w/ no A runs but B runs
dataAve.illusoryCutoff_incAllIllTaskSubj = dataAve.illusoryCutoff;
for iI = 1:length(dataAve.part_w_no_A)
    counter = length(dataAve.illusoryCutoff_incAllIllTaskSubj) + 1;
    dataAve.illusoryCutoff_incAllIllTaskSubj(counter,1) = dataHolderB(partHolderB == dataAve.part_w_no_A(iI)).partNum;
    dataAve.illusoryCutoff_incAllIllTaskSubj(counter,2) = str2double(dataHolderB(partHolderB == dataAve.part_w_no_A(iI)).date);
    dataAve.illusoryCutoff_incAllIllTaskSubj(counter,3) = str2double(dataHolderB(partHolderB == dataAve.part_w_no_A(iI)).time);
    dataAve.illusoryCutoff_incAllIllTaskSubj(counter,4) = 1;
end

% Run stats to look for differences between number of excluded participants
n_cat = 2;
% Make new array to exclude repeats in counts
[~, illCutoff_uniqueIdx, ~] = unique(dataAve.illusoryCutoff(:,1));
illCutoff_noRepeats = dataAve.illusoryCutoff(illCutoff_uniqueIdx,:);
samples = [sum(illCutoff_noRepeats(:,1)<2000000)...
    sum(illCutoff_noRepeats(:,1)<2000000 & illCutoff_noRepeats(:,4)==0);...
    sum(illCutoff_noRepeats(:,1)>=2000000 & illCutoff_noRepeats(:,1)<6000000)...
    sum(illCutoff_noRepeats(:,1)>=2000000 & illCutoff_noRepeats(:,1)<6000000 & illCutoff_noRepeats(:,4)==0);...
    sum(illCutoff_noRepeats(:,1)>=6000000)...
    sum(illCutoff_noRepeats(:,1)>=6000000 & illCutoff_noRepeats(:,4)==0)];
yates = 0;
[stats,data] = mpsContingencyTable(n_cat, samples, yates);

% Make a table listing exclusion %s. 
% % responses made > 4s after onset
% % participants responded to > 6 switches
% % participants responded to > 6 switches (w/in 4s stim onset)
% Total part excluded 
% dataAve.exclusion_table.percRespExcluded = [...
%     ((sum(dataAve.numResponses{1})-sum(dataAve.numResponsesCorrected{1}))...
%     /sum(dataAve.numResponses{1}) ) *100;...
%     ((sum(dataAve.numResponses{2})-sum(dataAve.numResponsesCorrected{2}))...
%     /sum(dataAve.numResponses{2}) ) *100;...
%     ((sum(dataAve.numResponses{3})-sum(dataAve.numResponsesCorrected{3}))...
%     /sum(dataAve.numResponses{3}) ) *100;];
% dataAve.exclusion_table.percPartIncluded
% dataAve.exclusion_table.

% Set the date this analysis was run
dataAve.dateRun = datestr(now);

save('SFMGroupData','datafull','dataAve');

end




