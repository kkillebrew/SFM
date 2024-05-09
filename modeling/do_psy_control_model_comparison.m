%% check matlab version
if datenum(version('-date')) < datenum('December 14, 2021')
    error(['This code is intended to run in MATLAB 2021b (or newer)...'])
end

%% Start parameters
    %Vanessa Morgan, 2/28/23

runControl = 1; %If you actually want to run the control or not
outputToCSV = 1; %If you want the script to append to the specified file
overwriteCSV = 1; % do you want to delete the existing .csv (if present)?
outFileName = 'ProbModelResults.csv'; %Name of .csv file to append to
matOutput = 'ProbModelBootstrapped_'; %Name of .mat file to save to
matOutDir = '/home/shaw-raid1/data/psychophysics/SFM_modeling_data';

outputToMat = 1; %If you want to output the model struct to a .mat file
    %NOTE: Only use this if a single variable (not including nBoot) is changed from default
    %If you've changed multiple variables, export the output manually

do_n_boot = 10000; % # times to boostrap

%% Empirical Data Values

              %Mean Switch Rate   Median Switch Count  Median Duration  Median CV
%Controls     .0902               11                   12               0.69
%Relatives    .15                 15                   7.2              0.67
%Patients     .166                17                   7.45             0.63

%% Create a struct of default options
defaultOptions =[];
defaultOptions.displayFigs = 1;%1 = on, 0 = off
defaultOptions.nLayers = 3;%1st = near depth, 2 = far depth, 3 = depth invariant 
defaultOptions.normGain = 1;% scale normalization factor
defaultOptions.tau = 1200;% set time constant - changed from 3250 KWK 20230912
defaultOptions.noiseAmp = 0.2;% arb. units
defaultOptions.noisefilter_t = 800;% set noise filtering in time
defaultOptions.sigma = 0.5;% semisaturation constant
defaultOptions.sigma_opp = 0.9;% semisaturation constant for opponency cells
defaultOptions.attnGain = 1.25;% scale attention factor - changed from 1 KWK 20230912

%% first do control model

%If you want to run model with default parameters
if(runControl)
    options = defaultOptions;
    options.nBoot = do_n_boot;  %Number of times to bootstrap
    options.bar_color = [0.5 1 0.5]; % color for histogram

    ctl_model_output = bootstrap_prob_model(options); %#ok<*UNRCH>
end

%% do psychosis model bootstrap
dnmG = defaultOptions.normGain;
dnsG = defaultOptions.noiseAmp;
datG = defaultOptions.attnGain;
dtau = defaultOptions.tau;
dsig = defaultOptions.sigma;
                            
all_params = [0.8   dnsG  datG  dtau  dsig;
              dnmG  0.24  datG  dtau  dsig; % Changed from .16 to .24 - KWK 20240405
              dnmG  dnsG  1.05  dtau  dsig; % Changed from .6 - KWK 20230912
              dnmG  dnsG  datG  800   dsig; % x / 3250 = 0.0902 / .166 = ~1770 (previously used 210) - KWK changed from 1400 20230912 (mess with number to get correct values)
              dnmG  dnsG  datG  dtau  0.4]; % n model versions x m parameters to vary - KWK changed from 0.6 20230912
% norm Gain
% noise Gain
% Attn Gain
% tau
% sigma

for iP = 1:size(all_params,1)
    
    options.nBoot = do_n_boot; % # times to boostrap
    options.normGain = all_params(iP,1); %Normalization gain
    options.noiseAmp = all_params(iP,2); % fixed bootstrap_prob_model to use this mps 2023.07.13
    options.attnGain = all_params(iP,3); % scale attention factor, .8 = PwP group
    options.tau = all_params(iP,4); % set time constant -- 3250 matches controls, 2100 matches SZ
    options.sigma = all_params(iP,5); %Semi-saturation constant
    options.sigma_opp = 0.9; %Semi-saturation constant (opponency cells), * does not get used *
    options.bar_color = [1 0.5 0.5]; % color for histogram
    options.nLayers = 3; %Number of layers in model
    %only specified in run_prob_model.m

    psy_model_output = bootstrap_prob_model(options);

    %% Table Preparation
    %Vanessa Morgan, 2/28/23
    
    %Copy the names of the fields of the 'changedOptions' struct
    optionNames = fieldnames(psy_model_output.changedOptions);
    
    %Removes the parameters 'bar_color' and 'displayFigs' from list of options
    %(These options aren't necessary and only complicate things)
    optionNames(strcmp(optionNames,'bar_color')) = [];
    optionNames(strcmp(optionNames,'displayFigs')) = [];
    
    %Output statistics for the table
    outVarNames = {'medianReversals', 'meanSwitchRate', ...
        'medianDuration', 'medianCV', ...
        'timeRun', 'changedOptions'};
    
    %Adds the names of each option to table
    outVarNames = vertcat(outVarNames(:), string(optionNames(:,1)));
    
    %Creates a string that holds names of changed variables
    changedFields = '';
    
    %Creates a 0 long table for data
    outputTable = cell2table(cell(0,numel(outVarNames)),'VariableNames', outVarNames(:));
    
    %Loops through each field in the options of the output
    for changed = 1:numel(optionNames)
        %If it has been changed from default
        if psy_model_output.changedOptions.(optionNames{changed}) ~= 0
            %Appends any changed variables
            changedFields = append(changedFields, [optionNames{changed} ' ']);
        end
    end
    
    %Trims trailing whitespace from the string
    changedFields = strtrim(changedFields);
    
    %Changes 'timeRun' and 'changedOptions' to String to allow for correct input
    outputTable = convertvars(outputTable, {'timeRun', 'changedOptions'}, 'string');
    
    %Assign variables to table
    outputTable.medianReversals(1) = psy_model_output.reversals.median; %This throws a warning you should ignore
    outputTable.meanSwitchRate(1) = psy_model_output.Hz.mean;
    outputTable.medianDuration(1) = psy_model_output.duration.median;
    outputTable.medianCV(1) = psy_model_output.CV.median;
    outputTable.timeRun(1) = string(psy_model_output.date_run);
    outputTable.changedOptions(1) = changedFields;
    
    %For each option
    for opt = 1:length(optionNames)
        %Assign value of option to the matching column in output
        outputTable.(optionNames{opt})(1) = psy_model_output.options.(optionNames{opt});
    end
    
    %If you're running a control model
    if(runControl) && (iP == 1)
        %Assign control variables to table
        outputTable.medianReversals(2) = ctl_model_output.reversals.median;
        outputTable.meanSwitchRate(2) = ctl_model_output.Hz.mean;
        outputTable.medianDuration(2) = ctl_model_output.duration.median;
        outputTable.medianCV(2) = ctl_model_output.CV.median;
        outputTable.timeRun(2) = string(ctl_model_output.date_run);
        outputTable.changedOptions(2) = '';
        %For each option
        for opt = 1:length(optionNames)
            outputTable.(optionNames{opt})(2) = ctl_model_output.options.(optionNames{opt});
        end
    end
    
    %% Writing table to CSV
    %Vanessa Morgan, 2/12/23
    
    %If you want to output model results and parameters to a .csv file
    if outputToCSV
        if overwriteCSV && (iP == 1)
            cmd = ['rm ' outFileName];
            system(cmd);
        end

        %Appends an additional row to the table of the .csv file already submitted
        writetable(outputTable, outFileName, 'WriteMode', 'Append');
    end
    
    %% Append to .mat file
    %Vanessa Morgan, 3/3/23
    %Currently needs modification as resulting file is 7.9 GB at 10k bootstraps
    
    %If you want to output to a .mat file
    if outputToMat
        %Remove 'nBoot ' from changedFields (otherwise, script would not work
        changedField = erase(changedFields, 'nBoot');
        changedField = strtrim(changedField); %Remove whitespace
        
        %If no value (besides nBoot) is changed
        if strcmp(changedField , '')
            matOutputFile = [matOutput  'Default'];
            %If changed value is greater than default value
        elseif psy_model_output.changedOptions.(changedField) > 0
            %Append matOut with the variable changed and it being 25% higher
            matOutputFile = [matOutput changedField '_20Higher'];
        else
            matOutputFile = [matOutput changedField '_20Lower'];
        end
        matOutputFile = fullfile(matOutDir,[matOutputFile '.mat']);
        save(matOutputFile, 'psy_model_output', '-v7.3') %Save file
    end
end