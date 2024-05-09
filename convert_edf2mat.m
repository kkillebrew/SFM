function [output] = convert_edf2mat(options)
% usage: [output] = combine_sfm_and_mrs(options)
%
% mps 20190910
%%
addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode/'));
addpath(genpath('/home/shaw-raid1/data/psychophysics/SFM.git'));
addpath(genpath('/home/shaw-raid1/data/psychophysics/eyetracking/iTrack.git'));
%% opt
if ~exist('options','var')
    options = [];
end

% Load in participant file names
options.curDur = pwd;
cd /home/shaw-raid1/data/psychophysics/eyetracking/
output.fileNames = dir('P*');

% For each participant go into their respective folders, find any SFM
% files, convert the edf to matlab struct, and save the new mat file in
% their folder.
for iI = 1:length(output.fileNames)
   
    % cd into participant folder
    cd(['./' output.fileNames(iI).name])
    
    % Find SFM files
    output.fileNames(iI).ETName = dir('SFM*.edf');
    
    % Check to see if there are no files
    if isempty(output.fileNames(iI).ETName)
    else
        for iJ = 1:length(output.fileNames(iI).ETName)
            % Convert edf 2 mat
            output.etData(iI,iJ) = iTrack(output.fileNames(iI).ETName(iJ).name);
            
            % Save file in participant folder under same name
%             save()
        end
    end
    
    % cd back to ET folder
    cd ../
    
end

cd(options.curDur)

end





