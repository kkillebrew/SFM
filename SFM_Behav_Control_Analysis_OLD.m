% Analysis for PHCP SFM behavioral control data. 

function [data] = SFM_Behav_Control_Analysis(subjID)

%% Load and organize the data
% subjID = '0001_20191003_0942_nfy_SFM_typeA_run01';

% Load in the control switch values\
cd ../
data.controlSwitch = readtable('SFMControlSwitches.xlsx');
cd ./SFMTestData
data.controlSwitch = table2cell(data.controlSwitch);

% Load in the rawdata for the participant
data.rawdata = readtable(subjID);
data.rawdata = table2cell(data.rawdata);

% Look only at block type 0 data
controlTrialIndex = [data.controlSwitch{:,3}]'==0;

% Pull out the response times
data.responseTimes = cell2mat(data.rawdata(:,2));
data.responseTimes(~controlTrialIndex) = [];   % Only look at control trials

% Pull out the response types - code them as 1 (left) vs 2 (right)
for i=1:size(data.rawdata,1)
    if strcmp(data.rawdata(i,3),'left')
        data.responseTypes(i) = 1;
    elseif strcmp(data.rawdata(i,3),'right')
        data.responseTypes(i) = 2;
    end
end
data.responseTypes = data.responseTypes';
data.responseTypes(~controlTrialIndex) = [];   % Only look at control trials

% Pull out the control switch times
data.controlSwitchTimes = cell2mat(data.controlSwitch(:,2));
data.controlSwitchTimes(~controlTrialIndex) = [];   % Only look at control trials

% Pull out the control response types - code them as 1 (left) vs 2 (right)
for i=1:size(data.controlSwitch,1)
    if strcmp(data.controlSwitch(i,1),'left')
        data.controlSwitchTypes(i) = 1;
    elseif strcmp(data.controlSwitch(i,1),'right')
        data.controlSwitchTypes(i) = 2;
    end
end
data.controlSwitchTypes = data.controlSwitchTypes';
data.controlSwitchTypes(~controlTrialIndex) = [];   % Only look at control trials



%% Analysis
% First compare the recorded response type with the actual switch type
% Compare lengths of the two arrays - if response is larger they responded 
% with too many switches. If response is smaller they responded with too few. 
data.responseTypeSizeCompare = length(data.responseTypes)-length(data.controlSwitchTypes);

if data.responseTypeSizeCompare == 0   % responses = actual switches
    
    % Check to see if the switches occur in the same order
    data.responseTypeCompare =  ~(data.controlSwitchTypes==data.responseTypes);
    data.responseTypeCompareSum = sum(data.responseTypeCompare);
    
    if data.responseTypeCompareSum == 0   % each response is = to its corresponding actual switch
        % Now that we know our response and actual switches are equal we
        % can simply compare the times for each response. 
        for i=1:length(data.responseTimes)
            data.responseDiff(i) = data.responseTimes(i)-data.controlSwitchTimes(i);
        end
            
    elseif data.responseTypeCompareSum >= 1   % somewhere in the list (or multiple places in the list) the participant made either an incorrect response, a double resonse, or an additional switch
    
    end
        
elseif data.responseTypeSizeCompare >= 1   % responses > actual switches


elseif data.responseTypeSizeCompare <= -1   % responses < actual switches
    
    
end

end




