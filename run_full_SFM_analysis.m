function output = run_full_SFM_analysis(options)
% usage: output = run_full_SFM_analysis(options)
%
% Run all associated SFM analysis, including behavioral (real switch and
% bistable task analysis), symptom, (both in summarize_SFM_results), MRS 
% (combine_sfm_and_mrs)
%
% KWK - 20240115

%% opts
if ~exist('options','var') % parameters
    options = [];
end
if ~isfield(options,'displayControlFigs')   % Display figures in the control task analysis
    options.displayControlFigs = 0; % 1 = on, 0 = off
end
if ~isfield(options,'displayFigs')   % plot analysis figures
    options.displayFigs = 0; % 1 = on, 0 = off
end
if ~isfield(options,'displayFigs_stats')   % plot the auto generated stats figures for anova / K-W / etc
    options.displayFigs_stats = 0; % 1 = on, 0 = off
end
if ~isfield(options,'plot_stats')   % plot stats on top of the graphs
    options.plot_stats = 0; % 1 = on, 0 = off
end
if ~isfield(options,'excludeTypeA')
    options.excludeTypeA = 1; % exclude subjects based on control task performance, 0 = no, 1 = yes
end
if ~isfield(options,'excludeRedcap')
    options.excludeRedcap = 1; % exclude subjects based on control task performance, 0 = no, 1 = yes
    addpath(genpath('/home/shaw-raid1/matlab_tools/COP_analysis.git')) % Add path for redcap excludion function (in COP_analysis.git)
end
if ~isfield(options,'normalize')
    options.normalize = 1;   % Normalize the data being analyzed (1=log; 2=sqrt; 0=non normalized)
end
if ~isfield(options,'normalize_plot')
    options.normalize_plot = 0;   % Normalize the data being plotted (1=log; 2=sqrt; 0=non normalized)
end
if ~isfield(options,'combinePatsRels')
    options.combinePatsRels = 0;   % Combine across patients and relatives
end
if ~isfield(options,'dateCutoff')
    options.dateCutoff = 0;   % Only take participants after a given date
    options.dateCutoffVal = 20210815;   % Don't take ppt after this date
    options.dateCutoffValConverted = datenum(num2str(options.dateCutoffVal),'yyyymmdd');
end
if ~isfield(options,'diagnosis_corr_data_type')
    options.diagnosis_corr_data_type = 1;   % Correlate using hz=0; log=1
end
if ~isfield(options,'diagnosis_corr_data_plot')
    options.diagnosis_corr_data_plot = 0;   % Plot correlations using hz=0; log=1
end
if ~isfield(options,'diagnosis_corr_include_run2')
    options.diagnosis_corr_include_run2 = 1;   % Include=1; exclude=0
end
if ~isfield(options,'subj_group_def')   % Define what groups you want to compare
    % 1=controls x probands x relatives
    % 2=controls x SZ x BP
    % 3=SZ x SCA x BP
    % 4=controls x SCZ+SCA x BP
    options.subj_group_def = 1;
%     options.subj_group_def = 2;
end
if ~isfield(options,'includeSubjWithNoRealSwitch')
    options.includeSubjWithNoRealSwitch = 1;   % Include participants w/out A runs
end
if ~isfield(options,'excludeBenzoUsers')
   options.excludeBenzoUsers = 0; 
end
if ~isfield(options,'makeDemoTable_switchRateSubjsOnly')
    options.makeDemoTable_switchRateSubjsOnly = 0; 
end
if ~isfield(options,'run_corr_w_acc_sr')
    options.run_corr_w_acc_sr = 0;
end

addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode'))
options.curDur = pwd;

%% First pull in SFM data
% Set all switches to plot supplemental figures
options.plotAccuracy = 1;
options.plotRT = 1;
options.plotAvePerDur_TowardAway = 1;
options.plotDistTotalResp = 1;
options.plotDistRespMade = 1;

if options.excludeTypeA
    [dataAve] = SFM_Behav_Control_Group_Analysis(options);
    if options.includeSubjWithNoRealSwitch == 0
        illusoryCutoff = dataAve.illusoryCutoff;
        [~,illusoryCutoffIdx] = unique(illusoryCutoff(:,1));
        illusoryCutoff = illusoryCutoff(illusoryCutoffIdx,:);
        options.illusoryCutoff = illusoryCutoff;
    elseif options.includeSubjWithNoRealSwitch == 1
        illusoryCutoff = dataAve.illusoryCutoff_incAllIllTaskSubj;
        [~,illusoryCutoffIdx] = unique(illusoryCutoff(:,1));
        illusoryCutoff = illusoryCutoff(illusoryCutoffIdx,:);
        options.illusoryCutoff = illusoryCutoff;
    end
    SFMdata = analyze_SFM_data(options);
else
    SFMdata = analyze_SFM_data(options);
end

options.SFMdata = SFMdata;
options.dataAve = dataAve;

%% Next, run switch rate / percept duration analysis
options.pull_in_sfm_data = 0;
options.plotAccuracy = 0;
options.plotRT = 0;
options.plotAvePerDur_TowardAway = 0;
options.plotDistTotalResp = 0;
options.plotDistRespMade = 0;
options.plotTestRetest = 1;
options.plotMedianRangeTestRetest = 1;
options.plotSwitchRateByBlock = 1;
options.plotHistOfPerDurs = 1; 
options.plotAveSwitchRate = 1;
options.plotAvePerDur = 1;
options.plotAveCoV = 1;
options.plotPsychSympCorrs = 1;

output = summarize_SFM_results(options);

%% Next, plot MRS correlations
% options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20220823_phcp_OCC_192subj_H2O_scaled.csv';
% In MPS script, set conditional plot_retest to 0 to not plot the retest
% values

%% Lastly, plot switch rates for controls, BP, scz
clear SFMdata dataAve
options.pull_in_sfm_data = 0;
options.plotAccuracy = 0;
options.plotRT = 0;
options.plotAvePerDur_TowardAway = 0;
options.plotDistTotalResp = 0;
options.plotDistRespMade = 0;

options.subj_group_def = 1;
if options.excludeTypeA
    [dataAve] = SFM_Behav_Control_Group_Analysis(options);
    if options.includeSubjWithNoRealSwitch == 0
        illusoryCutoff = dataAve.illusoryCutoff;
        [~,illusoryCutoffIdx] = unique(illusoryCutoff(:,1));
        illusoryCutoff = illusoryCutoff(illusoryCutoffIdx,:);
        options.illusoryCutoff = illusoryCutoff;
    elseif options.includeSubjWithNoRealSwitch == 1
        illusoryCutoff = dataAve.illusoryCutoff_incAllIllTaskSubj;
        [~,illusoryCutoffIdx] = unique(illusoryCutoff(:,1));
        illusoryCutoff = illusoryCutoff(illusoryCutoffIdx,:);
        options.illusoryCutoff = illusoryCutoff;
    end
    SFMdata = analyze_SFM_data(options);
else
    SFMdata = analyze_SFM_data(options);
end

options.SFMdata = SFMdata;
options.dataAve = dataAve;
options.pull_in_sfm_data = 0;
options.plotTestRetest = 0;
options.plotMedianRangeTestRetest = 0;
options.plotSwitchRateByBlock = 0;
options.plotHistOfPerDurs = 0; 
options.plotAveSwitchRate = 1;
options.plotAvePerDur = 0;
options.plotAveCoV = 0;
options.plotPsychSympCorrs = 0;

output = summarize_SFM_results(options);


end