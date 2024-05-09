% Analysis for PHCP SFM behavioral control data. 

function [data] = SFM_Behav_Control_Analysis(data)

%% Load and organize the data
% subjID = '0002_20191004_1401_nfy_SFM_typeA_run01';

SFM_git_dir = '/home/shaw-raid1/data/psychophysics/SFM.git/';

% Load in the control switch values and response data
data.controlSwitch = readtable(fullfile(SFM_git_dir,'SFMControlSwitches.xlsx'));
% mps 20220110 changing to full path
data.controlSwitch = table2cell(data.controlSwitch);

phcp_git_dir = '/home/shaw-raid1/data/pHCP.git/subjectResponseFiles/';
% Load in the rawdata for the participant
data.rawdata = readtable(fullfile(phcp_git_dir,data.fileName));
% mps 20220110 changing from using cd to full path
data.rawdata = table2cell(data.rawdata);

%% Organize the actual switch times
% Look only at block type 0 data
controlActualTrialIndex = [data.controlSwitch{:,3}]'==0;

% Pull out the control switch times
data.controlSwitchTimes = cell2mat(data.controlSwitch(:,2));
data.controlSwitchTimes(~controlActualTrialIndex) = [];   % Only look at control trials

% Pull out the control response types - code them as 1 (left) vs 2 (right)
for i=1:size(data.controlSwitch,1)
    if strcmp(data.controlSwitch(i,1),'left')
        data.controlSwitchTypes(i) = 1;
    elseif strcmp(data.controlSwitch(i,1),'right')
        data.controlSwitchTypes(i) = 2;
    end
end
data.controlSwitchTypes = data.controlSwitchTypes';
data.controlSwitchTypes(~controlActualTrialIndex) = [];   % Only look at control trials

% Create time bins that we can use to search for correct responses in
% Do we want a cutoff? What if they take 10 seconds to report? 6 seconds? - KWK
for i=1:length(data.controlSwitchTimes)
    if i<length(data.controlSwitchTimes)
        data.controlSwitchBins(i,1) = data.controlSwitchTimes(i);
        data.controlSwitchBins(i,2) = data.controlSwitchTimes(i+1);
    else
        data.controlSwitchBins(i,1) = data.controlSwitchTimes(i);
        data.controlSwitchBins(i,2) = 125;   % Add the max possible differnce? Is this time window available? - KWK
    end
end

%% Organize the response data
% Look only at block type 0 data
controlResponseTrialIndex = [data.rawdata{:,1}]'==0;

% Pull out the response times
data.responseTimes = cell2mat(data.rawdata(:,2));
data.responseTimes(~controlResponseTrialIndex) = [];   % Only look at control trials


% Pull out the response types - code them as 1 (left) vs 2 (right)
for i=1:size(data.rawdata,1)
    if strcmp(data.rawdata(i,3),'left')
        data.responseTypes(i) = 1;
    elseif strcmp(data.rawdata(i,3),'right')
        data.responseTypes(i) = 2;
    end
end
data.responseTypes = data.responseTypes';
data.responseTypes(~controlResponseTrialIndex) = [];   % Only look at control trials

% Organize into cell arrays with responses within correct time bins
counter=1;
idxCheck = find((data.responseTimes(1)>data.controlSwitchBins(:,1) &...
            data.responseTimes(1)<=data.controlSwitchBins(:,2)));
% Loop through the control switch times and find all the values in the
% response times that fill within the time windows specified in control
% switch times.
for i=1:length(data.controlSwitchTimes)
    
    % Check which response times fall in the control time bin
    idxCheck = find(data.responseTimes(:)>data.controlSwitchBins(i,1) &...
        data.responseTimes(:)<data.controlSwitchBins(i,2));
    
    data.responseTimesCell{i}(1:length(idxCheck)) = data.responseTimes(idxCheck);
    data.responseTypesCell{i}(1:length(idxCheck)) = data.responseTypes(idxCheck);
    
    clear idxCheck
end

%% Analysis

% Count the number of bins the participant made at least one response w/in
% 3 s of physical stimulus onset. 
responseBinCount = 0;
responseBinLessThanCount = 0;
for i=1:length(data.responseTimesCell)
    % Only count if the bin isn't empty and at least one of the responses
    % is towards the actual direction of motion.
    if ~isempty(data.responseTypesCell{i}) && sum((data.responseTypesCell{i}(:) == data.controlSwitchTypes(i))) >= 1
        responseBinCount = responseBinCount+1;
        
        % If response is made w/in the first 3 seconds after the physical
        % change happens count it otherwise don't. 
        if sum(((data.responseTimesCell{i}(data.responseTypesCell{i}(:) == data.controlSwitchTypes(i))) -...
                data.controlSwitchTimes(i)) <= data.reactionTimeCutoff) >= 1
            responseBinLessThanCount = responseBinLessThanCount+1;
        end
    end
end
data.responseBinNum = responseBinCount;
data.responseBinNumCorrected = responseBinLessThanCount;
clear responseBinCount
% Count the number of bins the participant made at least one response

% Make sure that data.responseTimesCell length is =
% data.controlSwitchTimes, otherwise throw a warning
if ~(data.responseBinNum==length(data.controlSwitchTimes))
   warning(sprintf('%s%d%s%d','Not enough response bins. Only ',data.responseBinNum,' out of ',...
       length(data.controlSwitchTimes))); 
end

% First double check that the responses are correct. Make an array where
% 1's are correct and 0's are incorrect. Check to see if any responses in
% the time bin are 'correct' or toward the physical direction. 
for i=1:length(data.responseTimesCell)
    if ~isempty(data.responseTypesCell{i})
        
        % Calculate acc as just a measure of whethere they made a response
        % during the time window. 
        counter = 0;
        accBreak = 0;
        while accBreak==0
            counter = counter+1;
            if data.responseTypesCell{i}(counter) == data.controlSwitchTypes(i)
                data.accList(i) = 1;
                data.accIdx(i) = counter;   % Record the idx of the correct response.
                accBreak = 1;
            else
                data.accList(i) = 0;
                data.accIdx(i) = 0;   % Record the idx as 0 if no corr response.
                if counter==length(data.responseTypesCell{i})
                    accBreak = 1;
                end
            end
        end
        
        % Also calculate acc using a cutoff time of 4 s (they must respond w/in 4 s)
        accBreak = 0;
        counter = 0;
        while accBreak==0
            counter = counter+1;
            if data.responseTypesCell{i}(counter) == data.controlSwitchTypes(i) &&...
                    (data.responseTimesCell{i}(counter) - data.controlSwitchTimes(i)) <= data.reactionTimeCutoff
                data.accListCorrected(i) = 1;
                accBreak = 1;
            else
                data.accListCorrected(i) = 0;
                if counter==length(data.responseTypesCell{i})
                    accBreak = 1;
                end
            end
        end
        
        
    else
        data.accList(i) = 0;
        data.accListCorrected(i) = 0;
        data.accIdx(i) = 0;   % Record the idx as 0 if no corr response.
    end
end

% Look at each actual switch and grab the responses that occured after that
% switch, but before the next switch. Only look at the first response (for
% now). Take the difference between that response and the actual switch. If
% there is no response store it as NaN.

% Response time
for i=1:length(data.accList)
    if data.accList(i) == 1
        data.responseDiff(i) = data.responseTimesCell{i}(data.accIdx(i))-data.controlSwitchTimes(i);
    elseif data.accList(i) == 0
        data.responseDiff(i) = NaN;
    end
end

% Accuracy
% Count up the number of non empty cells and divide by total number of
% cells
counter = 0;
for i=1:length(data.accList)
   if data.accList(i) == 1
       counter = counter+1;
   end
end
data.responseAcc = counter/length(data.responseTimesCell)*100;

counter = 0;
for i=1:length(data.accListCorrected)
   if data.accListCorrected(i) == 1
       counter = counter+1;
   end
end
data.responseAccCorrected = counter/length(data.responseTimesCell)*100;


save(sprintf('%s',fullfile(SFM_git_dir, 'SubjData', data.fileName)),'data');
% mps 20220110 changing from using cd to full path

end




