% Funciton to load in and convert edf files to .mat files using edf_converter. Converts edf files into an edf_converter
% object and saves them back to their original locations on the server. Must have the edf_converter toolbox on local device
% and be connected to the server.

function [output] = convert_edf2mat(options)

%
%%
% addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode/'));
addpath(genpath('E:/GitRepos/SFM.git'));
addpath(genpath('E:/GitRepos/edf-converter.git'));

%% opt
if ~exist('options','var')
    options = [];
end
if ~isfield('options','buttonEventTesting')
    options.buttonEventTesting = 1;
end

% Load in participant file names
options.curDur = pwd;
options.etDataDur = 'Z:/eyetracking_data/psychophysics_eyetracking/';
cd(options.etDataDur)
if options.buttonEventTesting == 0
    output.fileNames = dir('P*');
elseif options.buttonEventTesting == 1
    buttonEventTesting_fileName = 'KWK_ButtonEvent_Testing*';
    output.fileNames = dir(buttonEventTesting_fileName);
end

% For each participant go into their respective folders, find any SFM
% files, convert the edf to matlab struct, and save the new mat file in
% their folder.
for iI = 1:length(output.fileNames)

        % cd into participant folder
        cd(['./' output.fileNames(iI).name])

        % Find SFM edf files
        output.fileNames(iI).ETNameEdf = dir('SFM*.edf');
        % Find SFM Edf2Mat files
        output.fileNames(iI).ETNameEDF2Mat = dir('SFM*EDF2MAT.mat');
        % Find SFM et useable files
        output.fileNames(iI).ETNameUseableMat = dir('SFM*Useable.mat');

        % Check to see if there are no files
        if isempty(output.fileNames(iI).ETNameEdf)
        else
            % Create idx for edf files that don't have corresponding edf2mat files
            clear ETEdfConvIdxHolder
            if isempty(output.fileNames(iI).ETNameEDF2Mat)
                output.fileNames(iI).ETEdfConvIdx = ones([1 size(output.fileNames(iI).ETNameEdf,1)]);
            else
                output.fileNames(iI).ETEdfConvIdx = ones([1 size(output.fileNames(iI).ETNameEdf,1)]);
                for iK = 1:size(output.fileNames(iI).ETNameEDF2Mat,1)
                    ETEdfConvIdxHolder(iK) = find(strcmp(extractBefore(output.fileNames(iI).ETNameEDF2Mat(iK).name,'_EDF2MAT.'),...
                        extractBefore({output.fileNames(iI).ETNameEdf.name},'.')));
                end
                output.fileNames(iI).ETEdfConvIdx(ETEdfConvIdxHolder) = 0;
            end
            
            % Create idx for edf files that don't have corresponding useable mat files
            clear ETEdfUseableIdxHolder
            if isempty(output.fileNames(iI).ETNameUseableMat)
                output.fileNames(iI).ETEdfUseableIdx = ones([1 size(output.fileNames(iI).ETNameEdf,1)]);
            else
                output.fileNames(iI).ETEdfUseableIdx = ones([1 size(output.fileNames(iI).ETNameEdf,1)]);
                for iK = 1:size(output.fileNames(iI).ETNameUseableMat,1)
                    ETEdfUseableIdxHolder(iK) = find(strcmp(extractBefore(output.fileNames(iI).ETNameUseableMat(iK).name,'_Useable.'),...
                        extractBefore({output.fileNames(iI).ETNameEdf.name},'.')));
                end
                output.fileNames(iI).ETEdfUseableIdx(ETEdfUseableIdxHolder) = 0;
            end


            for iJ = 1:length(output.fileNames(iI).ETNameEdf)
                if output.fileNames(iI).ETEdfConvIdx(iJ) || output.fileNames(iI).ETEdfUseableIdx(iJ)
                    % Check write permissions
                    tempFileDir = [options.etDataDur output.fileNames(iI).name '/tempFile.txt'];
                    [fid,errmsg] = fopen(tempFileDir, 'w');
                    if ~isempty(errmsg)&&strcmp(errmsg,'Permission denied')
                        % Do not have write permission!
                        warning(sprintf('You do not have write permission to the folder (%s).',output.fileNames(iI).name));
                    else
                        % Convert edf 2 mat
                        if output.fileNames(iI).ETEdfConvIdx(iJ) % IF the file doesn't exist create it
                            etData = Edf2Mat(output.fileNames(iI).ETNameEdf(iJ).name);
                        else   % Otherwise load it in to create useable file
                            load(output.fileNames(iI).ETNameEDF2Mat(iJ).name);
                        end
                        
                        % Create a non 'Edf2Mat' object .mat file with all relevant data, that can be imported on the server
                        % (without needing the edf-converter toolbox)
                        etDataUseable.Events.eventCode = etData.Events;   % Event info (messages contains events labels and time stamps)
                        etDataUseable.Samples = etData.Samples;   % All data, includes position and pupil size, and many other unknonw vals - KWK 20211006
                        etDataUseable.RawEdf = etData.RawEdf;   % Raw data (converted from asc using edfConvert)
                        etDataUseable.PUPIL = etData.PUPIL;   % Pupil index vals
                        etDataUseable.EYES = etData.EYES;   % Eye index vals
                        etDataUseable.EVENT_TYPES = etData.EVENT_TYPES;   % List of event types and corresponding values
                        etDataUseable.origFileName = etData.filename;   % Filename of original edf file

                        % Save file in participant folder under same name
                        if output.fileNames(iI).ETEdfConvIdx(iJ)
                            save([extractBefore(output.fileNames(iI).ETNameEdf(iJ).name,'.') '_EDF2MAT'],'etData')
                        end
                        if output.fileNames(iI).ETEdfUseableIdx(iJ)
                            save([extractBefore(output.fileNames(iI).ETNameEdf(iJ).name,'.') '_Useable'],'etDataUseable')
                        end

                        % Clear holders
                        disp(['Converted edf for ' extractBefore(output.fileNames(iI).ETNameEdf(iJ).name,'.')]);
                        % Close the temp file created to check for permissions
                        fclose(fid);
                        delete(tempFileDir);
                        clear etData etDataUseable filePerm fid tempFileDir
                    end
                else
                    disp(['Mat file found for ' extractBefore(output.fileNames(iI).ETNameEdf(iJ).name,'.')]);
                end
            end
        end
    
    % cd back to ET folder
    cd ../
    
end

cd(options.curDur)

end


