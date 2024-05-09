function output = summarize_SFM_results(options)
% usage: output = summarize_SFM_results(options)
%
% mps c. 2019#	modified:   analyze_SFM_data.m

%% opts
if ~exist('options','var') % parameters
    options = [];
end
if ~isfield(options,'displayControlFigs')   % Display figures in the control task analysis
    options.displayControlFigs = 0; % 1 = on, 0 = off
end
if ~isfield(options,'displayFigs')   % plot analysis figures
    options.displayFigs = 1; % 1 = on, 0 = off
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
if ~isfield(options,'pull_in_sfm_data')
    options.pull_in_sfm_data = 1;
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
if ~isfield(options,'plotTestRetest')
    options.plotTestRetest = 1;
end
if ~isfield(options,'plotMedianRangeTestRetest')
    options.plotMedianRangeTestRetest = 1;
end
if ~isfield(options,'plotSwitchRateByBlock')
    options.plotSwitchRateByBlock = 1;
end
if ~isfield(options,'plotHistOfPerDurs')
    options.plotHistOfPerDurs = 1;
end
if ~isfield(options,'plotAveSwitchRate')
    options.plotAveSwitchRate = 1;
end
if ~isfield(options,'plotAvePerDur')
    options.plotAvePerDur = 1;
end
if ~isfield(options,'plotAveCoV')
    options.plotAveCoV = 1;
end
if ~isfield(options,'plotPsychSympCorrs')
    options.plotPsychSympCorrs = 1;
end

addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode'))
options.curDur = pwd;

%% pull in SFM data
if options.pull_in_sfm_data == 1
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
else
    % Set SFMdata and dataAve using values already in options
    SFMdata = options.SFMdata;
    dataAve = options.dataAve;
    options = rmfield(options,'SFMdata');
    options = rmfield(options,'dataAve');
end

%% plot all data, all subj.
n_subj = numel(SFMdata.subj_number);
n_blocks = [1 5]; % A & B respectively
type_labels = {'SFM Real Switch Task','SFM Bistable Task'};
colors = {'c','y'};

if options.displayFigs == 1
    
    figure;
    % Set figure size
    figSize.allData.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.allData.aspectRatio = [11.1585 4.9261];   % Aspect ratio
    figSize.allData.figSize = [0 0 ...
        figSize.allData.baseSize(4)*...
        (figSize.allData.aspectRatio(1)/figSize.allData.aspectRatio(2))...
        figSize.allData.baseSize(4)];   % Size/postion of fig
    for iSession = 1:size(SFMdata.Hz_flips,2)
        for iType = 1:size(SFMdata.Hz_flips,3)
            subplot(3,1,iType); hold on
            errorbar([1:n_subj]',squeeze(nanmean(SFMdata.Hz_flips(:,iSession,iType,:),4)),...
                squeeze(nanstd(SFMdata.Hz_flips(:,iSession,iType,:),0,4))/sqrt(n_blocks(iType)),...
                ['o' colors{iSession}])
            title(type_labels{iType},'fontsize',18)
            %set(gca,'XTick',1:n_subj,'XTickLabel',SFMdata.subj_number)
            ylabel(sprintf('%s\n%s','Rotation switch','rate (Hz)'),'fontsize',15)
            box off 
            set(gcf,'color','w')
            set(gca,'XColor','k','YColor','k')
        end
        
        set(gcf,'color','w')
        
        switch_diff = SFMdata.Hz_flips(:,:,1,:) - SFMdata.Hz_flips(:,:,2,:);
        subplot(3,1,3); hold on
        plot([0 n_subj+1],[0 0],'k--')
        errorbar([1:n_subj]',squeeze(nanmean(switch_diff(:,iSession,1,:),4)),...
            squeeze(nanstd(switch_diff(:,iSession,1,:),0,4))/sqrt(n_blocks(iType)),...
            ['o' colors{iSession}])
        title('Control - Illusory','fontsize',18)
        %set(gca,'XTick',1:n_subj,'XTickLabel',SFMdata.subj_number)
        ylabel(sprintf('%s\n%s','Switch rate','difference (Hz)'),'fontsize',15)
        xlabel('Subject','fontsize',15)
    end
    subplot(3,1,1)
    legend('session 1','session 2','fontsize',15)
    set(gcf,'Position', figSize.allData.figSize)
end


%% Define the groups
% 1 = controls, relatives, probands; 2 = controls, SZ, BP
    % 3 = SZ, schizoaffective (SCA), BP; 4 = healthy (con+rel),
    % SZ+SCA, bipolar,
group_def_opt = [];
group_def_opt.subj_group_def = options.subj_group_def;
group_def_opt.subj_number = SFMdata.subj_number;

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

options.group_def = group_def_out;


%% Look at coorelations between patients and their relatives for SFM switch rates. 

%%%%% UPDATE THIS W/ SWITCH TO ONLY RUN DURING SUBJGROUPDEF = 1

typeIdx = [1 2];   % Illusory condition

% Find group
% c_idx = find(SFMdata.subj_number < 2000000);   % Control indices
p_idx = find(SFMdata.subj_number >= 6000000);   % Patient indices
r_idx = find(SFMdata.subj_number >= 2000000 & SFMdata.subj_number < 6000000);   % Relative indices

% Find family IDs (grab last 5 digits 
p_fid = num2str(SFMdata.subj_number(p_idx));
p_fid = str2num(p_fid(:,3:7));
r_fid = num2str(SFMdata.subj_number(r_idx));
r_fid = str2num(r_fid(:,3:7));

% Find only patients who have corresponding relatives
counter = 0;
for i=1:length(r_fid)
    if ~isempty(find(p_fid==r_fid(i)))
        counter = counter+1;
        r_intersect(counter) = r_idx(i);
        p_intersect(counter) = p_idx(find(p_fid==r_fid(i)));
    end
end

if options.displayFigs
    figure(2); hold on
    suptitle(sprintf('%s\n\n','Correlation of Switch Rates Between Patients and Their Relatives'));
    % Set figure size
    figSize.testRe.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.testRe.aspectRatio = [5.5933 4.9806];   % Aspect ratio
    figSize.testRe.figSize = [0 0 ...
        figSize.testRe.baseSize(4)*...
        (figSize.testRe.aspectRatio(1)/figSize.testRe.aspectRatio(2))...
        figSize.testRe.baseSize(4)];   % Size/postion of fig
end

% Average across runs
for i=typeIdx
    r_ave = nanmean(SFMdata.Hz_flips(r_intersect,1,i,:),4);
    p_ave = nanmean(SFMdata.Hz_flips(p_intersect,1,i,:),4);
    
    % Run correlation
    output.relativePatientCorr.ICC(i) = ICC(1, 'k', [r_ave,p_ave]);
    
    % Plot the scatter plot
    if options.displayFigs
        plotTitle = {'Control Task','Illusory Task'};
        
        [poly_fit] = polyfit(r_ave, p_ave, 1);
        
        fit_x = [min(r_ave) max(r_ave)];
        fit_y = poly_fit(1).*fit_x + poly_fit(2);
        y_range = [min(p_ave) max(p_ave)];
        
        subplot(1,2,i)
        plot(fit_x,fit_y,'k-','linewidth',2)
        
        hold on
        
        plot(r_ave, p_ave, ...
            'ko','MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',8)
        
        if options.plot_stats == 1
            text(max(r_ave)*0.1,max(p_ave)*1.05,...
                ['n = ' num2str(numel(r_ave))],'fontsize',15)
            text(max(r_ave)*0.1,max(p_ave)*1,...
                ['ICC = ' num2str(round(100.*output.relativePatientCorr.ICC(i))/100)],'fontsize',15)
            % text(max(data_t1)*0.9,max(data_t2)*0.6,...
            %     ['r = ' num2str(round(100.*corr_retest_r)/100)],'fontsize',18)
            % text(max(data_t1)*0.9,max(data_t2)*0.5,...
            %     ['p = ' num2str(round(100.*corr_retest_p)/100)],'fontsize',18)
        end
        
        set(gca,'xlim',[0 max(r_ave)*1.1])   % Fix axis limits - KWK 20200507
        set(gca,'ylim',[0 max(p_ave)*1.1])
        set(gca,'fontsize',18,'xcolor','k','ycolor','k')
        set(gcf,'color','w')
        xlabel('Relatives switch rate (Hz)','color','k','fontsize',15)
        ylabel('Probands switch rate (Hz)','color','k','fontsize',15)
        title(plotTitle{i},'fontsize',18)
        set(gcf,'Position',figSize.testRe.figSize)
        
    end
end

%% look at test-retest
ill_idx = 2; % illusory condition
retest_idx = ~isnan(SFMdata.date_num(:,1)) & ...
    ~isnan(SFMdata.date_num(:,2));
retest_subjid = SFMdata.subj_number(retest_idx);
retest_group_idx(retest_subjid<2000000) = 2;   % Controls
retest_group_idx(retest_subjid>=2000000 & retest_subjid<6000000) = 2;   % Relatives
retest_group_idx(retest_subjid>=6000000) = 3;   % Patients

data_t1 = nanmean(SFMdata.Hz_flips(retest_idx,1,ill_idx,:),4);
data_t2 = nanmean(SFMdata.Hz_flips(retest_idx,2,ill_idx,:),4);
% Get ride of nans - KWK 20200507
new_data_t1 = data_t1(~isnan(data_t1) & ~isnan(data_t2));
new_data_t2 = data_t2(~isnan(data_t1) & ~isnan(data_t2));
new_retest_group_idx = retest_group_idx(~isnan(data_t1) & ~isnan(data_t2));
data_t1 = new_data_t1;
data_t2 = new_data_t2;
retest_group_idx = new_retest_group_idx;
output.illusory_task.retest.corr.type = 'Spearman';   % KWK chaged from pearson 20220609
[output.illusory_task.retest.corr.r, output.illusory_task.retest.corr.p] = ...
    corr(data_t1, data_t2, 'type', output.illusory_task.retest.corr.type);
output.illusory_task.retest.ICC = ICC(3, 'k', [data_t1 data_t2]);

if options.plotTestRetest
    [poly_fit] = polyfit(data_t1, data_t2, 1);

    fit_x = [min(data_t1) max(data_t1)];
    fit_y = poly_fit(1).*fit_x + poly_fit(2);
    y_range = [min(data_t2) max(data_t2)];

    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 12;
    statsFontSize = 12;
    % Set figure size
    figSize.testRe.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.testRe.aspectRatio = [3 3];   % Aspect ratio
    figSize.testRe.figSize = [0 0 ...
        figSize.testRe.aspectRatio];   % Size/postion of fig
    
    plot(fit_x,fit_y,'k-','linewidth',2)
    hold on
    plot(data_t1(retest_group_idx==2), data_t2(retest_group_idx==2), ...
        'bo','MarkerFaceColor','w','linewidth',2,...
        'MarkerSize',6)
    plot(data_t1(retest_group_idx==3), data_t2(retest_group_idx==3), ...
        'ro','MarkerFaceColor','w','linewidth',2,...
        'MarkerSize',6)

    %     if options.plot_stats == 1
    text(max(data_t1)*0.1,max(data_t2)*1.15,...
        ['ICC = ' num2str(round(100.*output.illusory_task.retest.ICC)/100)],'fontsize',statsFontSize)
    text(max(data_t1)*0.1,max(data_t2)*1.05,...
        ['n = ' num2str(numel(data_t1))],'fontsize',statsFontSize)    
    % text(max(data_t1)*0.9,max(data_t2)*0.6,...
    %     ['r = ' num2str(round(100.*corr_retest_r)/100)],'fontsize',18)
    % text(max(data_t1)*0.9,max(data_t2)*0.5,...
    %     ['p = ' num2str(round(100.*corr_retest_p)/100)],'fontsize',18)
    %     end
    
    set(gca,'xlim',[0 .5])   % Fix axis limits - KWK 20200507
    set(gca,'ylim',[0 .5])
    set(gca,'xcolor','k','ycolor','k')
    xlabel('Session 1 switch rate (Hz)','color','k','fontsize',axisTitleFontSize)
    ylabel('Session 2 switch rate (Hz)','color','k','fontsize',axisTitleFontSize)
    title(sprintf('%s\n%s','Test-Retest Reliability','of Switch Rates'),'fontsize',titleFontSize)
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.testRe.figSize,'color','w')
end


%% Average difference between test-retest
output.illusory_task.retest.date_diff = SFMdata.date_num(:,2) - SFMdata.date_num(:,1);
output.illusory_task.retest.mean = mean(output.illusory_task.retest.date_diff(retest_idx));
output.illusory_task.retest.median = median(output.illusory_task.retest.date_diff(retest_idx));
output.illusory_task.retest.range = [min(output.illusory_task.retest.date_diff(retest_idx)) ...
    max(output.illusory_task.retest.date_diff(retest_idx))];

% Look at mean between groups
output.illusory_task.retest.meanGroup(1) =...
    mean(output.illusory_task.retest.date_diff(grpIdx{1} & retest_idx));
output.illusory_task.retest.meanGroup(2) =...
    mean(output.illusory_task.retest.date_diff(grpIdx{2} & retest_idx));
output.illusory_task.retest.meanGroup(3) =...
    mean(output.illusory_task.retest.date_diff(grpIdx{3} & retest_idx));

% Look at std between groups
output.illusory_task.retest.stdGroup(1) =...
    std(output.illusory_task.retest.date_diff(grpIdx{1} & retest_idx));
output.illusory_task.retest.stdGroup(2) =...
    std(output.illusory_task.retest.date_diff(grpIdx{2} & retest_idx));
output.illusory_task.retest.stdGroup(3) =...
    std(output.illusory_task.retest.date_diff(grpIdx{3} & retest_idx));

% Look at median between groups
output.illusory_task.retest.medianGroup(1) =...
    median(output.illusory_task.retest.date_diff(grpIdx{1} & retest_idx));
output.illusory_task.retest.medianGroup(2) =...
    median(output.illusory_task.retest.date_diff(grpIdx{2} & retest_idx));
output.illusory_task.retest.medianGroup(3) =...
    median(output.illusory_task.retest.date_diff(grpIdx{3} & retest_idx));

% Look at range between groups
output.illusory_task.retest.rangeGroup(1,:) = [min(output.illusory_task.retest.date_diff(grpIdx{1} & retest_idx)) ...
    max(output.illusory_task.retest.date_diff(grpIdx{1} & retest_idx))];
if ~isnan(output.illusory_task.retest.medianGroup(2))
    output.illusory_task.retest.rangeGroup(2,:) = [min(output.illusory_task.retest.date_diff(grpIdx{2} & retest_idx)) ...
        max(output.illusory_task.retest.date_diff(grpIdx{2} & retest_idx))];
else
    output.illusory_task.retest.rangeGroup(2,1:2) = NaN;
end
output.illusory_task.retest.rangeGroup(3,:) = [min(output.illusory_task.retest.date_diff(grpIdx{3} & retest_idx)) ...
    max(output.illusory_task.retest.date_diff(grpIdx{3} & retest_idx))];

% Plot the median and range of date diff  - KWK 20200507
if options.plotMedianRangeTestRetest == 1
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 12;
    statsFontSize = 12;
    % Set figure size
    figSize.testRe.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.testRe.aspectRatio = [3 3];   % Aspect ratio
    figSize.testRe.figSize = [0 0 ...
        figSize.testRe.aspectRatio];   % Size/postion of fig
    
    % Plot as boxplot w/ beeswarm underneath
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    % Box plot
    boxPlot = boxplot(output.illusory_task.retest.date_diff);
    pause(0.5)
    set(boxPlot,'linewidth',2)
    % Beeswarm
    x_val = 1;
    bee_bin_width = .1;
    bee_spread_width = .5;
    beePlot = plotSpread(output.illusory_task.retest.date_diff, 'binWidth', bee_bin_width,...
        'distributionColors', {[.8 .8 .8]},...
        'xValues', x_val,...
        'spreadWidth', bee_spread_width);
    set(beePlot{1},'MarkerSize',15)
    % Set y axis to log scale
    set(gca,'YScale','log')
    set(gca,'ylim',[1 1000],'ytick',[10 100 250 500 1000])
    ylabel('Difference (days)','fontsize',axisTitleFontSize)
    xlabel('Participants','fontsize',axisTitleFontSize)
    set(gca,'xticklabel',{[]})
    title(sprintf('%s\n%s','Median and Range Between','Test and Re-Test'),'fontsize',titleFontSize)
    box off
    set(gca,'XColor','k','YColor','k')
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.testRe.figSize,'color','w')
end

output.illusory_task.retest.n_grp1 = numel(SFMdata.subj_number(retest_idx & grpIdx{1}));
output.illusory_task.retest.n_grp2 = numel(SFMdata.subj_number(retest_idx & grpIdx{2}));
output.illusory_task.retest.n_grp3 = numel(SFMdata.subj_number(retest_idx & grpIdx{3}));

%% Look at correlations between diagnostic test and swtich rate
if options.plotPsychSympCorrs
    % Grab diagnostic info
    % symptom_list = {'Mars contrast','BACS Composite Z Scores',...
    %     'BPRS Total Score','SGI',...
    %     'spq total','PID-5 psychoticism','SAPS Total Scores',...
    %     'SANS Total Scores'};
    % symptom_short = {'Mars','BACS','BPRS','SGI','SPQ','PID5','SAPS','SANS'};
    
    
    % Can try looking at PID5 and SGI first (possibly SANS and SAPS and MARS)
    symptom_list = {'BACS Composite Z Scores','BPRS Total Score','spq total'};
    symptom_short = {'BACS','BPRS','SPQ'};
    
%     symptom_list = {'BPRS Disorganization'};
%     symptom_short = {'BPRS_D'};
    
    if options.diagnosis_corr_include_run2==0
        sym_opt.subj_number = SFMdata.subj_number;
        sym_opt.date_number = SFMdata.date_num(:,1);
    elseif options.diagnosis_corr_include_run2==1
        sym_opt.subj_number = [SFMdata.subj_number, SFMdata.subj_number];
        sym_opt.date_number = SFMdata.date_num(:,:);
    end
    
    for iS = 1:numel(symptom_list)
        sym_opt.symptom_measure = symptom_list{iS};
        % Call MPS get_phcp_symptoms function
        cd /home/shaw-raid1/matlab_tools/mpsCode/
        sym_out = get_phcp_symptoms(sym_opt);
        cd(options.curDur);
        symptoms.(symptom_short{iS}).data = sym_out.psy_list;
    end
    
    % Grab/average Hz data
    if options.diagnosis_corr_include_run2==0   % Exclude session 2
        holderHz = squeeze(nanmean(SFMdata.Hz_flips(:,1,:,:),4));
        holderHzLog = log10(squeeze(nanmean(SFMdata.Hz_flips(:,1,:,:),4)));
        holderSubjNum = [SFMdata.subj_number];
        holderGrpIdx{1} = [grpIdx{1}];
        holderGrpIdx{2} = [grpIdx{2}];
        holderGrpIdx{3} = [grpIdx{3}];
        
        % Make a run list, to exclude all run 2 subjects during correlation
        runIdx = ones([length(holderHz),1]);
    elseif options.diagnosis_corr_include_run2==1   % Include session 2
        holderHz = [squeeze(nanmean(SFMdata.Hz_flips(:,1,:,:),4));...
            squeeze(nanmean(SFMdata.Hz_flips(:,2,:,:),4))];
        holderHzLog = [log10(squeeze(nanmean(SFMdata.Hz_flips(:,1,:,:),4)));...
            log10(squeeze(nanmean(SFMdata.Hz_flips(:,2,:,:),4)))];
        holderSubjNum = [SFMdata.subj_number; SFMdata.subj_number];
        holderGrpIdx{1} = [grpIdx{1}; grpIdx{1}];
        holderGrpIdx{2} = [grpIdx{2}; grpIdx{2}];
        holderGrpIdx{3} = [grpIdx{3}; grpIdx{3}];
        
        % Make a run list, to exclude all run 2 subjects during correlation
        runIdx = [ones([length(holderHz)/2,1]); zeros([length(holderHz)/2,1])];
    end
    
    groupColorArray = options.group_def.use_colors_RGB;
    
    % Correlate dagnosis values with Hz in illusory task
    for iS = 1:numel(symptom_list)
        
        % Concatinate a second list of diagnostic values to the end of
        % original, to match the concat'd second list of Hz values. All 'extra'
        % subjects (that don't have a second rund) will be excluded b/c of
        % nanIdx (where all values in second concat'd list will be nans, where
        % no second session is run).
        clear holderSymptoms
        holderSymptoms = symptoms.(symptom_short{iS}).data;
        %     if options.diagnosis_corr_include_run2==1   % Include session 2
        %         holderSymptoms = [holderSymptoms;...
        %             holderSymptoms];
        %     end
        
        % Remove nans to do correlations
        clear nanIdx
        nanIdx = isnan(holderHz);
        nanIdx(:,1) = isnan(holderSymptoms) | nanIdx(:,1);
        nanIdx(:,2) = isnan(holderSymptoms) | nanIdx(:,2);
        
        % Create group index arrays to plot points in group specific colors
        group_idx = zeros([1,size(holderSubjNum(nanIdx(:,2)~=1),1)])';
        grp1_idx = find(holderGrpIdx{1}(nanIdx(:,2)~=1));
        grp2_idx = find(holderGrpIdx{2}(nanIdx(:,2)~=1));
        grp3_idx = find(holderGrpIdx{3}(nanIdx(:,2)~=1));
        group_idx(grp1_idx) = 1;
        group_idx(grp2_idx) = 2;
        group_idx(grp3_idx) = 3;
        
        % Choose what data to plot (hz vs log)
        if options.diagnosis_corr_data_plot==1
            data_t1_plot = holderHzLog(nanIdx(:,2)~=1,2);
        elseif options.diagnosis_corr_data_plot==0
            data_t1_plot = holderHz(nanIdx(:,2)~=1,2);
        end
        if options.diagnosis_corr_data_type==1
            data_t1_corr = holderHzLog([nanIdx(:,2)~=1 & runIdx],2);
            data_t1_lme = holderHzLog(nanIdx(:,2)~=1,2);
        elseif options.diagnosis_corr_data_type==0
            data_t1_corr = holderHz([nanIdx(:,2)~=1 & runIdx],2);
            data_t1_lme = holderHz(nanIdx(:,2)~=1,2);
        end
        data_t2 = holderSymptoms(nanIdx(:,2)~=1);
        data_t2_corr = holderSymptoms(nanIdx(:,2)~=1 & runIdx);
        
        % Run correlation
        output.illusory_task.diagnosis.(symptom_short{iS}).corr.type = 'Spearman';
        [output.illusory_task.diagnosis.(symptom_short{iS}).corr.r, output.illusory_task.diagnosis.(symptom_short{iS}).corr.p] = ...
            corr(data_t1_corr, data_t2_corr, 'type', output.illusory_task.diagnosis.(symptom_short{iS}).corr.type);
        
        % Run linear mixed effects model
        output.illusory_task.diagnosis.(symptom_short{iS}).lme.lmeTable = table(...
            data_t1_lme, holderSubjNum(nanIdx(:,2)~=1),group_idx,data_t2,...
            'VariableNames',{'Switch','Subj','Group',symptom_short{iS}});
        % Turn the categorical variables into 'categorical' type
        output.illusory_task.diagnosis.(symptom_short{iS}).lme.lmeTable.Subj = categorical(...
            output.illusory_task.diagnosis.(symptom_short{iS}).lme.lmeTable.Subj);
        output.illusory_task.diagnosis.(symptom_short{iS}).lme.lmeTable.Group = categorical(...
            output.illusory_task.diagnosis.(symptom_short{iS}).lme.lmeTable.Group);
        
        %     output.illusory_task.diagnosis.(symptom_short{iS}).lme.formula =...
        %         ['Switch ~ Group + ' symptom_short{iS} ' + (1 | Subj)'];
        output.illusory_task.diagnosis.(symptom_short{iS}).lme.formula =...
            ['Switch ~ ' symptom_short{iS} ' + (1 | Subj)'];
        %     output.illusory_task.diagnosis.(symptom_short{iS}).lme.formula =...
        %         ['Switch ~ Group * ' symptom_short{iS} ' + (1 | Subj)'];
        output.illusory_task.diagnosis.(symptom_short{iS}).lme.output =...
            fitlme(output.illusory_task.diagnosis.(symptom_short{iS}).lme.lmeTable,...
            output.illusory_task.diagnosis.(symptom_short{iS}).lme.formula);
        
        [poly_fit] = polyfit(data_t2, data_t1_plot, 1);
        
        fit_x = [min(data_t2) max(data_t2)];
        fit_y = poly_fit(1).*fit_x + poly_fit(2);
        x_range = [min(data_t2) max(data_t2)];
        if options.diagnosis_corr_data_plot==1
            y_range = [-2 0];
        elseif options.diagnosis_corr_data_plot==0
            y_range = [0 0.5];
        end
        
        figure; hold on
        % Set font sizes
        titleFontSize = 12;
        axisTitleFontSize = 12;
        axisLabelFontSize = 12;
        statsFontSize = 10;
        % Set figure size
        figSize.diag.baseSize = get(0,'Screensize');   % Base size in pixels
        figSize.diag.aspectRatio = [3 3];   % Aspect ratio
        figSize.diag.figSize = [0 0 ...
            figSize.diag.aspectRatio];   % Size/postion of fig
        
        plot(fit_x,fit_y,'k-','linewidth',2)
        
        plot(data_t2(grp1_idx), data_t1_plot(grp1_idx), ...
            'o','Color',groupColorArray{1},'MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        plot(data_t2(grp2_idx), data_t1_plot(grp2_idx), ...
            'o','Color',groupColorArray{2},'MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        plot(data_t2(grp3_idx), data_t1_plot(grp3_idx), ...
            'o','Color',groupColorArray{3},'MarkerFaceColor','w','linewidth',2,...
            'MarkerSize',6)
        % Set y axis to log scale
        set(gca,'YScale','log')
        
        if options.plot_stats == 1
            text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
                .02,...
                ['n = ' num2str(sum(group_idx~=0))],'fontsize',statsFontSize)
        end
        text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
            .02,...
            ['r = ' sprintf('%1.3f',output.illusory_task.diagnosis.(symptom_short{iS}).corr.r)],'fontsize',statsFontSize)
        text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
            .015,...
            ['p = ' sprintf('%1.3f',output.illusory_task.diagnosis.(symptom_short{iS}).corr.p)],'fontsize',statsFontSize)
        if options.plot_stats == 1
            % Find the correct index for the LME symptom output
            lmeOutputIdx = strcmp(output.illusory_task.diagnosis.(symptom_short{iS}).lme.output.Coefficients.Name,symptom_short(iS));
            text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
                .005,...
                ['t(', sprintf('%d',output.illusory_task.diagnosis.(symptom_short{iS}).lme.output.Coefficients.DF(lmeOutputIdx)), ') = ' ...
                sprintf('%1.3f',output.illusory_task.diagnosis.(symptom_short{iS}).lme.output.Coefficients.tStat(lmeOutputIdx)) ...
                ', p=' sprintf('%1.3f',output.illusory_task.diagnosis.(symptom_short{iS}).lme.output.Coefficients.pValue(lmeOutputIdx))],'fontsize',statsFontSize)
        end
        % text(max(data_t1)*0.9,max(data_t2)*0.6,...
        %     ['r = ' num2str(round(100.*corr_retest_r)/100)],'fontsize',18)
        % text(max(data_t1)*0.9,max(data_t2)*0.5,...
        %     ['p = ' num2str(round(100.*corr_retest_p)/100)],'fontsize',18)
        %         end
        
        set(gca,'xlim',x_range)   % Fix axis limits - KWK 20200507
        set(gca,'ylim',[0.007 0.5],'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
        set(gca,'fontsize',axisLabelFontSize,'xcolor','k','ycolor','k')
        set(gcf,'color','w')
        if options.diagnosis_corr_data_plot==1
            ylabel('Switch rate (Log10(Hz))','color','k','fontsize',axisTitleFontSize)
        elseif options.diagnosis_corr_data_plot==0
            ylabel('Switch rate (Hz)','color','k','fontsize',axisTitleFontSize)
        end
        xlabel(symptom_short{iS},'color','k','fontsize',axisTitleFontSize)
        %         title(sprintf('%s','Relationship between diagnostic test (',symptom_short{iS},') and switch rate (Hz)'),...
        %             'fontsize',18)
        
        set(gcf,'Units','inches')
        set(gcf,'Position',figSize.diag.figSize,'color','w')
    end
end

%% Look at the data by block
blockHolder = squeeze(nanmean(SFMdata.Hz_flips(:,:,2,:),2));   % Average across retests

% Create grouping arrays
clear grp1_idx grp2_idx grp3_idx
grp1_idx = grpIdx{1};
grp2_idx = grpIdx{2};
grp3_idx = grpIdx{3};

blockGroupAve.controlsFull(:,:) = blockHolder(grp1_idx,:);
blockGroupAve.relativesFull(:,:) = blockHolder(grp2_idx,:);
blockGroupAve.patientsFull(:,:) = blockHolder(grp3_idx,:);

blockGroupAve.controlsAve(:) = nanmean(blockGroupAve.controlsFull,1);
blockGroupAve.controlsSte(:) = (nanstd(blockGroupAve.controlsFull,1)) ./ sqrt(length(blockGroupAve.controlsFull));
blockGroupAve.relativesAve(:) = nanmean(blockGroupAve.relativesFull,1);
blockGroupAve.relativesSte(:) = (nanstd(blockGroupAve.relativesFull,1)) ./ sqrt(length(blockGroupAve.relativesFull));
blockGroupAve.patientsAve(:) = nanmean(blockGroupAve.patientsFull,1);
blockGroupAve.patientsSte(:) = (nanstd(blockGroupAve.patientsFull,1)) ./ sqrt(length(blockGroupAve.patientsFull));

% Plot 
if options.plotSwitchRateByBlock
    figure; hold on
    % Set font sizes
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 12;
    statsFontSize = 10;
    % Set figure size
    figSize.blockAve.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.blockAve.aspectRatio = [6.5 5];   % Aspect ratio
    figSize.blockAve.figSize = [0 0 ...
        figSize.blockAve.aspectRatio];   % Size/postion of fig
    figSize.blockAve.lineOffset = 0.3;
    
    addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
    
    xAxis = [1:2:length(blockGroupAve.controlsAve')*2]';
    
    % Beeswarm
    for iI=1:length(blockGroupAve.controlsAve)
        bee_bin_width = .1;
        bee_spread_width = .5;
        beePlot = plotSpread({blockGroupAve.controlsFull(:,iI),blockGroupAve.relativesFull(:,iI),blockGroupAve.patientsFull(:,iI)},...
            'binWidth', bee_bin_width,...
            'distributionColors', {options.group_def.corr_colors{1},...
            options.group_def.corr_colors{2},options.group_def.corr_colors{3}},...
            'xValues', [xAxis(iI)-figSize.blockAve.lineOffset, xAxis(iI), xAxis(iI)+figSize.blockAve.lineOffset],...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
        hold on
    end
    
    % Lineplots
    hb1 = plot(xAxis-figSize.blockAve.lineOffset, blockGroupAve.controlsAve', 'o--', 'Color', options.col_list{1}, 'MarkerEdgeColor',options.col_list{1});
    hold on
    hb2 = plot(xAxis, blockGroupAve.relativesAve', 'o--', 'Color', options.col_list{2}, 'MarkerEdgeColor',options.col_list{2});
    hb3 = plot(xAxis+figSize.blockAve.lineOffset, blockGroupAve.patientsAve', 'o--', 'Color', options.col_list{3}, 'MarkerEdgeColor',options.col_list{3});
    
    errorbar(xAxis-figSize.blockAve.lineOffset,blockGroupAve.controlsAve',blockGroupAve.controlsSte,'.','Color', options.col_list{1});
    errorbar(xAxis,blockGroupAve.relativesAve',blockGroupAve.relativesSte,'.','Color', options.col_list{2});
    errorbar(xAxis+figSize.blockAve.lineOffset,blockGroupAve.patientsAve',blockGroupAve.patientsSte,'.','Color', options.col_list{3});
        
    title('SFM Bistable Task Switch Rate by Block','fontsize',18)
    set(gca,'xticklabels',{'1','2','3','4','5'},'xtick',xAxis,'xLim',[0 10]);
    set(gca,'yscale','log')
    set(gca,'ylim',[0.007 0.5],'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
    box off
    ylabel('Switch rate (Hz)','fontsize',15)
    xlabel('Block number','fontsize',15)
    set(gca,'fontsize',15,'xcolor','k','ycolor','k')
    set(gca,'XColor','k','YColor','k')
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.blockAve.figSize,'color','w')
end

%% Plot historgram of all percept durations across all blocks for each switch across all subjs
% Create grouping arrays
clear grp1_idx grp2_idx grp3_idx
grp1_idx = grpIdx{1};
grp2_idx = grpIdx{2};
grp3_idx = grpIdx{3};
grp1_idx_val = find(grp1_idx==1);
grp2_idx_val = find(grp2_idx==1);
grp3_idx_val = find(grp3_idx==1);

% Concatinate all the perc durations and put them into one array for each subject
clear percDurHolder
percDurHolder = cell([size(SFMdata.percept_dur,1) 2]);
for iS=1:size(SFMdata.percept_dur,1)   % For all subects
    for iR=1:size(SFMdata.percept_dur,2)   % For repeated runs
        for iB=1:size(SFMdata.percept_dur,4)
            percDurHolder{iS,1} = [percDurHolder{iS,1}(:); SFMdata.percept_dur{iS,iR,1,iB}];
            percDurHolder{iS,2} = [percDurHolder{iS,2}(:); SFMdata.percept_dur{iS,iR,2,iB}];
        end
    end
end

% Combine across subjects into groups
clear percDurHolderFull
percDurHolderFull = cell([3,2]);
for iS=1:length(grp1_idx_val)
    percDurHolderFull{1,1} = [percDurHolderFull{1,1}(:); percDurHolder{grp1_idx_val(iS),1}];
    percDurHolderFull{1,2} = [percDurHolderFull{1,2}(:); percDurHolder{grp1_idx_val(iS),2}];
end
for iS=1:length(grp2_idx_val)
    percDurHolderFull{2,1} = [percDurHolderFull{2,1}(:); percDurHolder{grp2_idx_val(iS),1}];
    percDurHolderFull{2,2} = [percDurHolderFull{2,2}(:); percDurHolder{grp2_idx_val(iS),2}];
end
for iS=1:length(grp3_idx_val)
    percDurHolderFull{3,1} = [percDurHolderFull{3,1}(:); percDurHolder{grp3_idx_val(iS),1}];
    percDurHolderFull{3,2} = [percDurHolderFull{3,2}(:); percDurHolder{grp3_idx_val(iS),2}];
end
% 
% % Plot histograms
% figure()
% subplot(1,3,1)% Set figure size
% histogram(percDurHolderFull{1,1},'BinWidth',1,'BinLimits',[0,60]);
% subplot(1,3,2)
% histogram(percDurHolderFull{2,1},'BinWidth',1,'BinLimits',[0,60]);
% subplot(1,3,3)
% histogram(percDurHolderFull{3,1},'BinWidth',1,'BinLimits',[0,60]);
% 

if options.plotHistOfPerDurs
    titleFontSize = 12;
    axisTitleFontSize = 12;
    axisLabelFontSize = 12;
    plotKuipSwitch = 0;
    
    %% Plot illusory task
    figure()
    figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
    figSize.switchRate.aspectRatio = [6.5 2];   % Aspect ratio
    figSize.switchRate.figSize = [0 0 ...
        figSize.switchRate.baseSize(3)*.75 ...
        (figSize.switchRate.baseSize(3)*.75)*...
        (figSize.switchRate.aspectRatio(2)/figSize.switchRate.aspectRatio(1))];   % Size/postion of fig
    
    subplot(1,2,1)
    histRealG1 = histogram(percDurHolderFull{1,1},'BinWidth',.5,'BinLimits',[0,60],'Normalization','probability','DisplayStyle','stairs');
    hold on
    set(histRealG1,'EdgeColor',options.group_def.use_colors{1})
    histRealG2 = histogram(percDurHolderFull{2,1},'BinWidth',.5,'BinLimits',[0,60],'Normalization','probability','DisplayStyle','stairs');
    set(histRealG2,'EdgeColor',options.group_def.use_colors{2})
    histRealG3 = histogram(percDurHolderFull{3,1},'BinWidth',.5,'BinLimits',[0,60],'Normalization','probability','DisplayStyle','stairs');
    axis([0 20 0 0.15]);
    set(histRealG3,'EdgeColor',options.group_def.use_colors{3})
    
    % Kolmogorov-Smirnov test
    [output.control_task.Hz.KSTest.h(1) output.control_task.Hz.KSTest.p(1) ...
        output.control_task.Hz.KSTest.stats(1)] = kstest2(percDurHolderFull{1,1},percDurHolderFull{2,1});
    [output.control_task.Hz.KSTest.h(2) output.control_task.Hz.KSTest.p(2) ...
        output.control_task.Hz.KSTest.stats(2)] = kstest2(percDurHolderFull{1,1},percDurHolderFull{3,1});
    [output.control_task.Hz.KSTest.h(3) output.control_task.Hz.KSTest.p(3) ...
        output.control_task.Hz.KSTest.stats(3)] = kstest2(percDurHolderFull{2,1},percDurHolderFull{3,1});
    
    % Kuipers's test
    [output.control_task.Hz.KupTest.stats(1) output.control_task.Hz.KupTest.stats(1)] = ...
        kuipertest2(percDurHolderFull{1,1},...
        percDurHolderFull{2,1},numel(percDurHolderFull{1,1}),plotKuipSwitch);
    [output.control_task.Hz.KupTest.stats(2) output.control_task.Hz.KupTest.stats(2)] = ... 
        kuipertest2(percDurHolderFull{1,1},...
        percDurHolderFull{3,1},numel(percDurHolderFull{1,1}),plotKuipSwitch);
    [output.control_task.Hz.KupTest.stats(3) output.control_task.Hz.KupTest.stats(3)] = ...
        kuipertest2(percDurHolderFull{2,1},...
        percDurHolderFull{3,1},numel(percDurHolderFull{1,1}),plotKuipSwitch);
    
    % Plot stats
    if options.plot_stats == 1
        text(9.5,0.13,...
            sprintf('%s%s%s%s%d%s%d%s%.3f%s%.3f',grpLabelShort{1},'v',grpLabelShort{2},...
            ': D(',size(grp1_idx_val,1),',',size(grp2_idx_val,1),')=',output.control_task.Hz.KSTest.stats(1),...
            ', p=',output.control_task.Hz.KSTest.p(1)),'FontSize',axisLabelFontSize);
        text(9.5,0.12,...
            sprintf('%s%s%s%s%d%s%d%s%.3f%s%.3f',grpLabelShort{1},'v',grpLabelShort{3},...
            ': D(',size(grp1_idx_val,1),',',size(grp3_idx_val,1),')=',output.control_task.Hz.KSTest.stats(2),...
            ', p=',output.control_task.Hz.KSTest.p(2)),'FontSize',axisLabelFontSize);
        text(9.5,0.11,...
            sprintf('%s%s%s%s%d%s%d%s%.3f%s%.3f',grpLabelShort{2},'v',grpLabelShort{3},...
            ': D(',size(grp2_idx_val,1),',',size(grp3_idx_val,1),')=',output.control_task.Hz.KSTest.stats(3),...
            ', p=',output.control_task.Hz.KSTest.p(3)),'FontSize',axisLabelFontSize);
    end
    
    set(gca,'YLim',[0 .15],'fontsize',axisLabelFontSize)
    title(sprintf('%s\n%s','Histogram of Percept Durations',...
        '(Real Switch Task)'),...
        'fontsize',titleFontSize)
    box off
    ylabel('Number of Percepts Reported','fontsize',axisTitleFontSize)
    xlabel('Time (s)','fontsize',axisTitleFontSize)
    set(gca,'XColor','k','YColor','k')
    
    
    %% Plot illusory task
    %     figure()
    %     figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
    %     figSize.switchRate.aspectRatio = [10.9849 9.2814];   % Aspect ratio
    %     figSize.switchRate.figSize = [0 0 ...
    %         3.54 3];   % Size/postion of fig
    %     set(gcf,'color','w')
    %     set(gca,'XColor','k','YColor','k')
    %     set(gcf,'units','inches')
    %     set(gcf,'Position', figSize.switchRate.figSize)

    subplot(1,2,2)
    histIllG1 = histogram(percDurHolderFull{1,2},'BinWidth',.5,'BinLimits',[0,60],'Normalization','probability','DisplayStyle','stairs');
    hold on
    set(histIllG1,'EdgeColor',options.group_def.use_colors{1})
    histIllG2 = histogram(percDurHolderFull{2,2},'BinWidth',.5,'BinLimits',[0,60],'Normalization','probability','DisplayStyle','stairs');
    set(histIllG2,'EdgeColor',options.group_def.use_colors{2})
    histIllG3 = histogram(percDurHolderFull{3,2},'BinWidth',.5,'BinLimits',[0,60],'Normalization','probability','DisplayStyle','stairs');
    axis([0 20 0 0.15]);
    set(histIllG3,'EdgeColor',options.group_def.use_colors{3})
    
    % Kolmogorov-Smirnov test
    [output.illusory_task.Hz.KSTest.h(1) output.illusory_task.Hz.KSTest.p(1) ...
        output.illusory_task.Hz.KSTest.stats(1)] = kstest2(percDurHolderFull{1,2},percDurHolderFull{2,2});
    [output.illusory_task.Hz.KSTest.h(2) output.illusory_task.Hz.KSTest.p(2) ...
        output.illusory_task.Hz.KSTest.stats(2)] = kstest2(percDurHolderFull{1,2},percDurHolderFull{3,2});
    [output.illusory_task.Hz.KSTest.h(3) output.illusory_task.Hz.KSTest.p(3) ...
        output.illusory_task.Hz.KSTest.stats(3)] = kstest2(percDurHolderFull{2,2},percDurHolderFull{3,2});
    
    % Kuipers's test
    [output.illusory_task.Hz.KupTest.stats(1) output.illusory_task.Hz.KupTest.stats(1)] = ...
        kuipertest2(percDurHolderFull{1,2},...
        percDurHolderFull{2,2},numel(percDurHolderFull{1,2}),plotKuipSwitch);
    [output.illusory_task.Hz.KupTest.stats(2) output.illusory_task.Hz.KupTest.stats(2)] = ... 
        kuipertest2(percDurHolderFull{1,2},...
        percDurHolderFull{3,2},numel(percDurHolderFull{1,2}),plotKuipSwitch);
    [output.illusory_task.Hz.KupTest.stats(3) output.illusory_task.Hz.KupTest.stats(3)] = ...
        kuipertest2(percDurHolderFull{2,2},...
        percDurHolderFull{3,2},numel(percDurHolderFull{1,2}),plotKuipSwitch);
    
    % Plot stats
    if options.plot_stats == 1
        text(9.5,0.13,...
            sprintf('%s%s%s%s%d%s%d%s%.3f%s%.3f',grpLabelShort{1},'v',grpLabelShort{2},...
            ': D(',size(grp1_idx_val,1),',',size(grp2_idx_val,1),')=',output.illusory_task.Hz.KSTest.stats(1),...
            ', p=',output.illusory_task.Hz.KSTest.p(1)),'FontSize',axisLabelFontSize);
        text(9.5,0.12,...
            sprintf('%s%s%s%s%d%s%d%s%.3f%s%.3f',grpLabelShort{1},'v',grpLabelShort{3},...
            ': D(',size(grp1_idx_val,1),',',size(grp3_idx_val,1),')=',output.illusory_task.Hz.KSTest.stats(2),...
            ', p=',output.illusory_task.Hz.KSTest.p(2)),'FontSize',axisLabelFontSize);
        text(9.5,0.11,...
            sprintf('%s%s%s%s%d%s%d%s%.3f%s%.3f',grpLabelShort{2},'v',grpLabelShort{3},...
            ': D(',size(grp2_idx_val,1),',',size(grp3_idx_val,1),')=',output.illusory_task.Hz.KSTest.stats(3),...
            ', p=',output.illusory_task.Hz.KSTest.p(3)),'FontSize',axisLabelFontSize);
    end
    
    set(gca,'YLim',[0 .15],'fontsize',axisLabelFontSize)
    title(sprintf('%s\n%s','Histogram of Percept Durations',...
        '(Bistable Task)'),...
        'fontsize',titleFontSize)
    box off
    ylabel('Number of Percepts Reported','fontsize',axisTitleFontSize)
    xlabel('Time (s)','fontsize',axisTitleFontSize)
    set(gca,'XColor','k','YColor','k')
    
    set(gcf,'Units','inches')
    set(gcf,'Position',figSize.blockAve.figSize,'color','w')
    
end

%% Create demographics table using only subjects included in switch rate analysis
if options.makeDemoTable_switchRateSubjsOnly == 1
    demogOpts.subj_number = [SFMdata.subj_number(grpIdx{1}); ...
        SFMdata.subj_number(grpIdx{2}); ...
        SFMdata.subj_number(grpIdx{3})];
    demogOpts.date_number = [SFMdata.date_num(grpIdx{1},1); ...
        SFMdata.date_num(grpIdx{2},1); ...
        SFMdata.date_num(grpIdx{3},1)];
    
    output.demoTalbe_switchRateSubjOnly = make_phcp_methods_table(demogOpts);
end

%% break up by group

% Make a list of subjects that are taking benzo's for exclusion - KWK 20231215
if options.excludeBenzoUsers == 1
    SFMdata.excludeBenzoIdx = zeros([length(SFMdata.subj_number) 1]);
    for iI = 1:length(SFMdata.subj_number)
        if sum(SFMdata.subj_number(iI) == dataAve.benzoExclusionList') >= 1
            SFMdata.excludeBenzoIdx(iI) = 1; 
        end
    end
else
    SFMdata.excludeBenzoIdx = zeros([length(SFMdata.subj_number) 1]);
end

clear grp1_idx grp3_idx grp2_idx
grp1_idx = find(grpIdx{1});
grp2_idx = find(grpIdx{2});
grp3_idx = find(grpIdx{3});
% Make index to identify what data to use if not all subjects are included
total_grp_idx_val = [grp1_idx; grp2_idx; grp3_idx];
total_grp_idx = [ones([numel(grp1_idx) 1]); ones([numel(grp2_idx) 1])+1; ...
    ones([numel(grp3_idx) 1])+2];

% Exclude benzo users - KWK 20231215
if options.excludeBenzoUsers == 1
    total_grp_idx(SFMdata.excludeBenzoIdx==1) = [];
    total_grp_idx_val(SFMdata.excludeBenzoIdx==1) = [];
end

% Update the grouping variables to index from the correct participant list
% (E.g. don't index from the full participant list if you are using a subset)
clear grp1_idx grp3_idx grp2_idx
grp1_idx = find(total_grp_idx==1);
grp2_idx = find(total_grp_idx==2);
grp3_idx = find(total_grp_idx==3);

if options.plotAveSwitchRate
    h = figure;
end
if options.plotAvePerDur
    j = figure;
end
if options.plotAveCoV
    o = figure;
end

% Add in subject numbers and date numbers for each group to the output struct
output.subj_number{1} = SFMdata.subj_number(grpIdx{1});
output.subj_number{2} = SFMdata.subj_number(grpIdx{2});
output.subj_number{3} = SFMdata.subj_number(grpIdx{3});

output.date_number{1} = SFMdata.date_num(grpIdx{1},:);
output.date_number{2} = SFMdata.date_num(grpIdx{2},:);
output.date_number{3} = SFMdata.date_num(grpIdx{3},:);

task_names = {'control_task','illusory_task'};

for iTask = 1:size(SFMdata.Hz_flips,3)
    
    %% First look at Hz
    % first look at an ANOVA, which allows us to take test-retest into
    % account, as well as using data from each block
    % -- log transform does NOT make data more normal -- does for
    % patients but makes controls and relatives look more skewed
    % Add in switch to plot normalized data or non normalized data
    % if options.normalize=1 (log) options.normalize=2 (sq root)
    % options.normalize=0 (nothing)
    if options.normalize == 0
        all_data = squeeze(SFMdata.Hz_flips(:,:,iTask,:));
    elseif options.normalize == 1
        all_data = squeeze(log10(SFMdata.Hz_flips(:,:,iTask,:)));
        all_data(all_data == -Inf | all_data == Inf) = NaN;
    elseif options.normalize == 2
        all_data = squeeze(sqrt(SFMdata.Hz_flips(:,:,iTask,:)));
    end
    % Only grab subjects to be included
    all_data = all_data(total_grp_idx_val,:,:);
    output.(task_names{iTask}).Hz.all_data = all_data;
    
    if options.normalize_plot == 0
        all_data_plot = squeeze(SFMdata.Hz_flips(:,:,iTask,:));
    elseif options.normalize_plot == 1
        all_data_plot = squeeze(log10(SFMdata.Hz_flips(:,:,iTask,:)));
        all_data_plot(all_data_plot == -Inf | all_data_plot == Inf) = NaN;
    elseif options.normalize_plot == 2
        all_data_plot = squeeze(sqrt(SFMdata.Hz_flips(:,:,iTask,:)));
    end
    % Only grab subjects to be included
    all_data_plot = all_data_plot(total_grp_idx_val,:,:);
    
    all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2) size(all_data,3)]);
    if options.combinePatsRels == 0
        all_group = [ones(numel(grp1_idx), size(all_data,2), size(all_data,3)) ; ...
            2*ones(numel(grp2_idx), size(all_data,2), size(all_data,3)) ; ...
            3*ones(numel(grp3_idx), size(all_data,2), size(all_data,3)) ];
    elseif options.combinePatsRels == 1
        all_group = [ones(numel(grp1_idx), size(all_data,2), size(all_data,3)) ; ...
            2*ones(numel(grp2_idx), size(all_data,2), size(all_data,3)) ; ...
            2*ones(numel(grp3_idx), size(all_data,2), size(all_data,3)) ];
    end
    all_retest = repmat([1 2],[size(all_data,1) 1 size(all_data,3)]);
    all_block = [];
    for iB = 1:size(all_data,3)
        all_block = cat(3, all_block, iB*ones(size(all_data,1), size(all_data,2)));
    end
    output.(task_names{iTask}).Hz.all_group = all_group;
    
    x_labels = {[grpLabelShort{1} ', n=' num2str(sum(~isnan(nanmean(nanmean(all_data(grp1_idx,:,:),3),2))))],...
        [grpLabelShort{2} ', n=' num2str(sum(~isnan(nanmean(nanmean(all_data(grp2_idx,:,:),3),2))))],...
        [grpLabelShort{3} ', n=' num2str(sum(~isnan(nanmean(nanmean(all_data(grp3_idx,:,:),3),2))))]};
    
    nest = zeros(3,3);
    nest(1,2) = 1;
    
    if options.displayFigs_stats
        show_stats_fig = 'on';
    else
        show_stats_fig = 'off';
    end
        
    if iTask == 2 % only do ANOVA for the illusory task, which has different blocks
        [panova_Hz, output.(task_names{iTask}).Hz.anova_table] = anovan(all_data(:),{all_subj(:),...
            all_group(:),all_block(:)},'random',1,'continuous',[3],...
            'nested',nest,'model','full','varnames',{'subj','group','block'},...
            'display',show_stats_fig);
    end
    
    % Set averaged all_data to plot and do further 
    Hz_data = nanmean(nanmean(all_data,3),2);
    Hz_data_plot = nanmean(nanmean(all_data_plot,3),2);
    
    % do a k-w test, because data may not be normally distributed & equal
    % variance across groups... (but can't use test-retest)
    comp_group_idx = [1 2 ; 1 3 ; 2 3]; % c vs. r, c vs. p, r vs. p
    grouping = [ones(1,numel(grp1_idx)) 2*ones(1,numel(grp2_idx)) 3*ones(1,numel(grp3_idx))];
%     grouping = [zeros(1,numel(grp1_idx)) zeros(1,numel(grp2_idx)) zeros(1,numel(grp3_idx))];
%     grouping(grp1_idx) = 1;
%     grouping(grp2_idx) = 2;
%     grouping(grp3_idx) = 3;

    % Run Levene's test for equal variance
    [output.(task_names{iTask}).Hz.levene.p, output.(task_names{iTask}).Hz.levene.stats] =...
        vartestn(Hz_data,grouping,'display','off','testtype','LeveneAbsolute');

    % KW and ttest2 between individual groups
    for iComp = 1:size(comp_group_idx,1)
        % KW
        [p_group_Hz_kw(iTask, iComp) , ...
            output.(task_names{iTask}).Hz.kruskall_wallis_2_groups{iComp}.table , ...
            output.(task_names{iTask}).Hz.kruskall_wallis_2_groups{iComp}.stats ] = ...
            kruskalwallis([ Hz_data( grouping == ...
            comp_group_idx(iComp,1) ) ; Hz_data( grouping == ...
            comp_group_idx(iComp,2) ) ] , [ grouping( grouping == ...
            comp_group_idx(iComp,1))' ; grouping( grouping == ...
            comp_group_idx(iComp,2))' ], 'off');
        output.(task_names{iTask}).Hz.kruskall_wallis_2_groups{iComp}.p = ...
            p_group_Hz_kw(iTask, iComp);
        
        % Ttest2
        [~, p_group_Hz_ttest(iTask, iComp) , ...
            output.(task_names{iTask}).Hz.t_test_2{iComp}.CI , ...
            output.(task_names{iTask}).Hz.t_test_2{iComp}.stats ] = ...
            ttest2( all_data( grouping == ...
            comp_group_idx(iComp,1) ) , all_data( grouping == ...
            comp_group_idx(iComp,2) ) );
        output.(task_names{iTask}).Hz.t_test_2{iComp}.p = p_group_Hz_ttest(iTask, iComp);
    end
    
    % Reset all_data to average for k-w  between 3 groups
    % then do k-w for Hz
    clear all_data
    if options.normalize==0
        all_data = [nanmean(nanmean(SFMdata.Hz_flips(grp1_idx,:,iTask,:),4),2) ; ...
            nanmean(nanmean(SFMdata.Hz_flips(grp2_idx,:,iTask,:),4),2) ; ...
            nanmean(nanmean(SFMdata.Hz_flips(grp3_idx,:,iTask,:),4),2)];
    elseif options.normalize==1
        all_data = [nanmean(nanmean(log10(SFMdata.Hz_flips(grp1_idx,:,iTask,:)),4),2) ; ...
            nanmean(nanmean(log10(SFMdata.Hz_flips(grp2_idx,:,iTask,:)),4),2) ; ...
            nanmean(nanmean(log10(SFMdata.Hz_flips(grp3_idx,:,iTask,:)),4),2)];
    end
    
    [p_full_Hz_kw(iTask), output.(task_names{iTask}).Hz.kruskall_wallis_3_groups.table , ...
        output.(task_names{iTask}).Hz.kruskall_wallis_3_groups.stats] = ...
        kruskalwallis(all_data, grouping, show_stats_fig);
    output.(task_names{iTask}).Hz.kruskall_wallis_3_groups.p = p_full_Hz_kw(iTask);
    
    % Create array of switches over time for control subjects to plot for
    % 7T Methods paper.
    if iTask == 2
        xAxis = 0:.001:120;
        flipState = 1;
        counter = 1;
        for iI = 1:length(grp1_idx)
            if ~isempty(SFMdata.flip_time{grp1_idx(iI),1,2,5})
                flipTimes_data_idx{counter} = SFMdata.flip_time{grp1_idx(iI),1,2,5};
                for iJ = 1:length(xAxis)
                    if sum(xAxis(iJ) == round(flipTimes_data_idx{counter},3)) > 0
                        flipState = 3 - flipState;
                    end
                    flipTime_data(counter,iJ) = flipState;
                end
                counter = counter+1;
            end
        end
    end
    
    % Look at effect size between controls and PwPP - 20231211
    cohenDOpt.mean(1) = nanmean(Hz_data(grouping == 1));
    cohenDOpt.mean(2) = nanmean(Hz_data(grouping == 3));
    cohenDOpt.sds(1) = nanstd(Hz_data(grouping == 1));
    cohenDOpt.sds(2) = nanstd(Hz_data(grouping == 3));
    cohenDOpt.n(1) = numel(Hz_data(grouping == 1));
    cohenDOpt.n(2) = numel(Hz_data(grouping == 3));
    output.(task_names{iTask}).Hz.effectSize = mpsCohensD(cohenDOpt.mean, cohenDOpt.sds, cohenDOpt.n);
    
    % Run correlation between accuracy and swtich rates - KWK 20231212
    if options.run_corr_w_acc_sr
        if iTask == 1
            % Find accuracy values
            accPartNumHolder = dataAve.subj_number;
            accValHolder = [dataAve.responseAccNoCuttoffCorr{1}; ...
                dataAve.responseAccNoCuttoffCorr{2}; ...
                dataAve.responseAccNoCuttoffCorr{3}];
            
            % First make index to grab the correct participants that are
            % included in both the switch rate analysis and have acc data
            clear partAccHolderIdx
            partAccHolderIdx = zeros([1 length(accPartNumHolder)]);
            for iS = 1:length(SFMdata.subj_number)
                if sum(SFMdata.subj_number(iS)==accPartNumHolder) >= 1
                    firstPartHolder = find(SFMdata.subj_number(iS)==accPartNumHolder);
                    partAccHolderIdx(firstPartHolder(1)) = 1;
                end
            end
            clear partSRHolderIdx
            partSRHolderIdx = zeros([1 length(SFMdata.subj_number)]);
            for iS = 1:length(SFMdata.subj_number)
                if sum(SFMdata.subj_number(iS)==accPartNumHolder(partAccHolderIdx==1)) >= 1
                    partSRHolderIdx(iS) = 1;
                end
            end
            
            % Make list of acc and SR values
            subjList = [SFMdata.subj_number(partSRHolderIdx==1) ...
                dataAve.subj_number(partAccHolderIdx==1)];
            groupingList(subjList(:,1) < 2000000) = 1;
            groupingList(subjList(:,1) >= 2000000 & subjList(:,1) < 6000000) = 2;
            groupingList(subjList(:,1) >= 6000000) = 3;
            srList = Hz_data(partSRHolderIdx==1);
            srList_Plot = Hz_data_plot(partSRHolderIdx==1);
            accList = accValHolder(partAccHolderIdx==1);
            
            %% Run correlations
            % Run across groups
            [output.control_task.corr_w_accuracy.all.r, output.control_task.corr_w_accuracy.all.p] = ...
                corr(srList, accList, 'type', output.illusory_task.retest.corr.type);
            
            % Now run correlations between groups
            for iG = 1:3
                [output.control_task.corr_w_accuracy.group(iG).r, output.control_task.corr_w_accuracy.group(iG).p] = ...
                    corr(srList(groupingList==iG), accList(groupingList==iG), 'type', output.illusory_task.retest.corr.type);
            end
            
            % Plot the correlations
            figure(); hold on
            % Set font sizes
            titleFontSize = 12;
            axisTitleFontSize = 12;
            axisLabelFontSize = 10;
            statsFontSize = 10;
            % Set figure size
            figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
            figSize.switchRate.aspectRatio = [6.5 5];   % Aspect ratio
            figSize.switchRate.figSize = [0 0 ...
                figSize.switchRate.aspectRatio];   % Size/postion of fig
            
            % Fit all data
            [poly_fit_all] = polyfit(accList, srList_Plot, 1);
            
            fit_x_all = [min(accList) max(accList)];
            fit_y_all = poly_fit_all(1).*fit_x_all + poly_fit_all(2);
            x_range_all = [min(accList) max(accList)];
            %             if options.diagnosis_corr_data_plot==1
            %                 y_range_all = [-2 0];
            %             elseif options.diagnosis_corr_data_plot==0
            %                 y_range_all = [0 0.5];
            %             end
            
            % Fit each group correlation
            for iG=1:3
                [poly_fit_group{iG}] = polyfit(accList(groupingList==iG), srList_Plot(groupingList==iG), 1);
                
                fit_x_group{iG} = [min(accList(groupingList==iG)) max(accList(groupingList==iG))];
                fit_y_group{iG} = poly_fit_group{iG}(1).*fit_x_group{iG} + poly_fit_group{iG}(2);
                x_range_group{iG} = [min(accList(groupingList==iG)) max(accList(groupingList==iG))];
                %                 if options.diagnosis_corr_data_plot==1
                %                     y_range = [-2 0];
                %                 elseif options.diagnosis_corr_data_plot==0
                %                     y_range = [0 0.5];
                %                 end
            end
            
            plot(accList(groupingList==1), srList_Plot(groupingList==1), ...
                'o','Color',groupColorArray{1},'MarkerFaceColor','w','linewidth',2,...
                'MarkerSize',6)
            hold on
            plot(accList(groupingList==2), srList_Plot(groupingList==2), ...
                'o','Color',groupColorArray{2},'MarkerFaceColor','w','linewidth',2,...
                'MarkerSize',6)
            plot(accList(groupingList==3), srList_Plot(groupingList==3), ...
                'o','Color',groupColorArray{3},'MarkerFaceColor','w','linewidth',2,...
                'MarkerSize',6)
            
            plot(fit_x_all,fit_y_all,'k-','linewidth',2)
            
            groupColorText = {'g-','b-','r-'};
            for iG=1:3
               plot(fit_x_group{iG},fit_y_group{iG},groupColorText{iG},'linewidth',2) 
            end
            
            % Set y axis to log scale
            set(gca,'YScale','log')
        
            %         if options.plot_stats == 1
            %             text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
            %                 .02,...
            %                 ['n = ' num2str(sum(group_idx~=0))],'fontsize',statsFontSize)
            %         end
            %         text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
            %             .02,...
            %             ['r = ' sprintf('%1.3f',output.illusory_task.diagnosis.(symptom_short{iS}).corr.r)],'fontsize',statsFontSize)
            %         text(max(data_t2)-((max(data_t2)-min(data_t2))*0.3),...
            %             .015,...
            %             ['p = ' sprintf('%1.3f',output.illusory_task.diagnosis.(symptom_short{iS}).corr.p)],'fontsize',statsFontSize)
        
            %             set(gca,'xlim',x_range)   % Fix axis limits - KWK 20200507
            set(gca,'ylim',[0.007 0.5],'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
            set(gca,'fontsize',axisLabelFontSize,'xcolor','k','ycolor','k')
            set(gcf,'color','w')
            ylabel('Switch rate (Hz)','color','k','fontsize',axisTitleFontSize)
            xlabel('Accuracy (% Correct)','color','k','fontsize',axisTitleFontSize)
            %         title(sprintf('%s','Relationship between diagnostic test (',symptom_short{iS},') and switch rate (Hz)'),...
            %             'fontsize',18)
            
            set(gcf,'Units','inches')
            set(gcf,'Position',figSize.diag.figSize,'color','w')
            
        end
    end
    
    %% Plot switch rate (Hz)
    
    if options.plotAveSwitchRate
        figure(h); hold on
        % Set font sizes
        titleFontSize = 12;
        axisTitleFontSize = 12;
        axisLabelFontSize = 10;
        statsFontSize = 10;
        % Set figure size
        figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
        figSize.switchRate.aspectRatio = [6.5 5];   % Aspect ratio
        figSize.switchRate.figSize = [0 0 ...
            figSize.switchRate.aspectRatio];   % Size/postion of fig        
        
        subplot(1,2,iTask)
        addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
        
        % Boxplots
        hb{iTask} = boxplot(Hz_data_plot,grouping);
        pause(0.5)
        set(gca,'XTick',1:3,'XTickLabel',x_labels,'fontsize',axisLabelFontSize)
        set(hb{iTask},'linewidth',2)
        hb2 = findobj(gca,'type','line');
        hb3 = findobj(gca,'type','Outliers');
        for iHB = 1:size(hb{iTask},2)
            set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
            set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
            set(hb2((iHB)+3:3:end),'MarkerFaceColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'MarkerFaceColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'Color',options.col_list{4-iHB})
        end
        hold on
        
        % Beeswarm
        x_val = [1 2 3];        set(gca,'XColor','k','YColor','k')

        bee_bin_width = .1;
        bee_spread_width = .5;
        beePlot = plotSpread({Hz_data_plot(grouping==1),Hz_data_plot(grouping==2),Hz_data_plot(grouping==3)},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.8 .8 .8]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
        hold on
        hbCurr = findobj(gca,'type','line');
        for iHB = 1:size(hb{iTask},2)
            set(hbCurr((iHB)+3:3:end),'color',options.col_list{4-iHB})
            set(hbCurr((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
        end
        
        % Plot line at physical switch
        if iTask == 1
            if options.normalize_plot == 0
                plot([0 4],[0.09 0.09],'--k');
            elseif options.normalize_plot == 1
                plot([0 4],log10([0.09 0.09]),'--k');
            end
        end
        
        % Plot significance
%         max_Hz = max([Hz_data_plot(grouping==1); Hz_data_plot(grouping==2); Hz_data_plot(grouping==3)]);
        max_Hz = .5;
        if options.plot_stats == 1
            if iTask == 1
                % Plot 3-K-W
                text(1,max_Hz*.75,...
                    ['X2(' sprintf('%d',output.(task_names{iTask}).Hz.kruskall_wallis_3_groups.table{2,3})  ') = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).Hz.kruskall_wallis_3_groups.table{2,5}) ', p = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).Hz.kruskall_wallis_3_groups.table{2,6})],...
                    'fontsize',statsFontSize);
            elseif iTask == 2
                % Plot anova results
                text(1,max_Hz*.85,...
                    ['F(' sprintf('%d',output.(task_names{iTask}).Hz.anova_table{2,3}) ',' ...
                    sprintf('%d',output.(task_names{iTask}).Hz.anova_table{3,3}) ') = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).Hz.anova_table{3,6}) ', p = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).Hz.anova_table{3,7})],...
                    'fontsize',statsFontSize);
                % Plot post hoc 2-K-W for controls vs patients
                text(1,max_Hz*.75,...
                    ['X2(' sprintf('%d',output.(task_names{iTask}).Hz.kruskall_wallis_2_groups{2}.table{2,3})  ') = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).Hz.kruskall_wallis_2_groups{2}.table{2,5}) ', p = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).Hz.kruskall_wallis_2_groups{2}.table{2,6})],...
                    'fontsize',statsFontSize);
            end
        end
                 
        title(type_labels{iTask},'fontsize',titleFontSize)
        box off
        if iTask == 1
            ylabel('Switch rate (Hz)','fontsize',axisTitleFontSize)
        end
        set(gca,'YScale','log')
        set(gca,'ylim',[0.007 0.5],'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
        set(gca,'XColor','k','YColor','k')
        
        set(gcf,'Units','inches')
        set(gcf,'Position',figSize.blockAve.figSize,'color','w')
    end
    clear all_data
    
    %% Then do percept duration
    % first look at an ANOVA, which allows us to take test-retest into
    % account, as well as using data from each block
    % -- log transform does NOT make data more normal -- does for
    % patients but makes controls and relatives look more skewed
    % Add in switch to plot normalized data or non normalized data
    % if options.normalize=1 (log) options.normalize=2 (sq root)
    % options.normalize=0 (nothing)
    if options.normalize == 0
        all_data = squeeze(SFMdata.mean_percept_dur(:,:,iTask,:));
    elseif options.normalize == 1
        all_data = squeeze(log10(SFMdata.mean_percept_dur(:,:,iTask,:)));
        all_data(all_data == -Inf | all_data == Inf) = NaN;
    elseif options.normalize == 2
        all_data = squeeze(sqrt(SFMdata.mean_percept_dur(:,:,iTask,:)));
    end
    % Only grab subjects to be included
    all_data = all_data(total_grp_idx_val,:,:);
    output.(task_names{iTask}).duration.all_data = all_data;
    
    if options.normalize_plot == 0
        all_data_plot = squeeze(SFMdata.mean_percept_dur(:,:,iTask,:));
    elseif options.normalize_plot == 1
        all_data_plot = squeeze(log10(SFMdata.mean_percept_dur(:,:,iTask,:)));
        all_data_plot(all_data_plot == -Inf | all_data_plot == Inf) = NaN;
    elseif options.normalize_plot == 2
        all_data_plot = squeeze(sqrt(SFMdata.mean_percept_dur(:,:,iTask,:)));
    end
    % Only grab subjects to be included
    all_data_plot = all_data_plot(total_grp_idx_val,:,:);
    
    all_subj = repmat([1:size(all_data,1)]',[1 size(all_data,2) size(all_data,3)]);
    if options.combinePatsRels == 0
        all_group = [ones(numel(grp1_idx), size(all_data,2), size(all_data,3)) ; ...
            2*ones(numel(grp2_idx), size(all_data,2), size(all_data,3)) ; ...
            3*ones(numel(grp3_idx), size(all_data,2), size(all_data,3)) ];
    elseif options.combinePatsRels == 1
        all_group = [ones(numel(grp1_idx), size(all_data,2), size(all_data,3)) ; ...
            2*ones(numel(grp2_idx), size(all_data,2), size(all_data,3)) ; ...
            2*ones(numel(grp3_idx), size(all_data,2), size(all_data,3)) ];
    end
    all_retest = repmat([1 2],[size(all_data,1) 1 size(all_data,3)]);
    all_block = [];
    for iB = 1:size(all_data,3)
        all_block = cat(3, all_block, iB*ones(size(all_data,1), size(all_data,2)));
    end
    output.(task_names{iTask}).duration.all_group = all_group;
    
    x_labels = {[grpLabelShort{1} ', n=' num2str(sum(~isnan(nanmean(nanmean(all_data(grp1_idx,:,:),3),2))))],...
        [grpLabelShort{2} ', n=' num2str(sum(~isnan(nanmean(nanmean(all_data(grp2_idx,:,:),3),2))))],...
        [grpLabelShort{3} ', n=' num2str(sum(~isnan(nanmean(nanmean(all_data(grp3_idx,:,:),3),2))))]};
    
    nest = zeros(3,3);
    nest(1,2) = 1;
    
    if options.displayFigs_stats
        show_stats_fig = 'on';
    else
        show_stats_fig = 'off';
    end
        
    if iTask == 2 % only do ANOVA for the illusory task, which has different blocks
        [panova, output.(task_names{iTask}).duration.anova_table] = anovan(all_data(:),{all_subj(:),...
            all_group(:),all_block(:)},'random',1,'continuous',[3],...
            'nested',nest,'model','full','varnames',{'subj','group','block'},...
            'display',show_stats_fig);
    end
    
    % Set averaged all_data to plot and do further 
    dur_data = nanmean(nanmean(all_data,3),2);
    dur_data_plot = nanmean(nanmean(all_data_plot,3),2);
    
    % do a k-w test, because data may not be normally distributed & equal
    % variance across groups... (but can't use test-retest)
    comp_group_idx = [1 2 ; 1 3 ; 2 3]; % c vs. r, c vs. p, r vs. p
    grouping = [ones(1,numel(grp1_idx)) 2*ones(1,numel(grp2_idx)) 3*ones(1,numel(grp3_idx))];
    %     grouping = [zeros(1,numel(grp1_idx)) zeros(1,numel(grp2_idx)) zeros(1,numel(grp3_idx))];
    %     grouping(grp1_idx) = 1;
    %     grouping(grp2_idx) = 2;
    %     grouping(grp3_idx) = 3;
    
    % Run Levene's test for equal variance
    [output.(task_names{iTask}).duration.levene.p, output.(task_names{iTask}).duration.levene.stats] =...
        vartestn(dur_data,grouping,'display','off','testtype','LeveneAbsolute');
    
    % KW and ttest2 between individual groups
    for iComp = 1:size(comp_group_idx,1)
        % KW
        [p_group_duration_kw(iTask, iComp) , ...
            output.(task_names{iTask}).duration.kruskall_wallis_2_groups{iComp}.table , ...
            output.(task_names{iTask}).duration.kruskall_wallis_2_groups{iComp}.stats ] = ...
            kruskalwallis([ dur_data( grouping == ...
            comp_group_idx(iComp,1) ) ; dur_data( grouping == ...
            comp_group_idx(iComp,2) ) ] , [ grouping( grouping == ...
            comp_group_idx(iComp,1))' ; grouping( grouping == ...
            comp_group_idx(iComp,2))' ], 'off');
        output.(task_names{iTask}).duration.kruskall_wallis_2_groups{iComp}.p = ...
            p_group_duration_kw(iTask, iComp);
        
        % Ttest2
        [~, p_group_duration_ttest(iTask, iComp) , ...
            output.(task_names{iTask}).duration.t_test_2{iComp}.CI , ...
            output.(task_names{iTask}).duration.t_test_2{iComp}.stats ] = ...
            ttest2( all_data( grouping == ...
            comp_group_idx(iComp,1) ) , all_data( grouping == ...
            comp_group_idx(iComp,2) ) );
        output.(task_names{iTask}).duration.t_test_2{iComp}.p = p_group_duration_ttest(iTask, iComp);
    end
    
    % do k-w for mean percept duration
    clear all_data
    if options.normalize==0
        all_data = [nanmean(nanmean(SFMdata.mean_percept_dur(grp1_idx,:,iTask,:),4),2) ; ...
            nanmean(nanmean(SFMdata.mean_percept_dur(grp2_idx,:,iTask,:),4),2) ; ...
            nanmean(nanmean(SFMdata.mean_percept_dur(grp3_idx,:,iTask,:),4),2)];
    elseif options.normalize==1
        all_data = [nanmean(nanmean(log10(SFMdata.mean_percept_dur(grp1_idx,:,iTask,:)),4),2) ; ...
            nanmean(nanmean(log10(SFMdata.mean_percept_dur(grp2_idx,:,iTask,:)),4),2) ; ...
            nanmean(nanmean(log10(SFMdata.mean_percept_dur(grp3_idx,:,iTask,:)),4),2)];
    end
        
    [p_full_dur_kw(iTask), output.(task_names{iTask}).duration.kruskall_wallis_3_groups.table , ...
        output.(task_names{iTask}).duration.kruskall_wallis_3_groups.stats] = ...
        kruskalwallis(all_data, grouping, show_stats_fig);
    output.(task_names{iTask}).duration.kruskall_wallis_3_groups.p = p_full_dur_kw(iTask);
     
    
    %% Plot percept duration
    if options.plotAvePerDur
        figure(j); hold on
        % Set font sizes
        titleFontSize = 12;
        axisTitleFontSize = 12;
        axisLabelFontSize = 10;
        statsFontSize = 10;
        % Set figure size
        figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
        figSize.switchRate.aspectRatio = [6.5 5];   % Aspect ratio
        figSize.switchRate.figSize = [0 0 ...
            figSize.switchRate.aspectRatio];   % Size/postion of fig
        
        subplot(1,2,iTask)
        addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
        
        % Boxplots
        hb{iTask} = boxplot(dur_data_plot,grouping);
        pause(0.5)
        set(gca,'XTick',1:3,'XTickLabel',x_labels,'fontsize',axisLabelFontSize)
        set(hb{iTask},'linewidth',2)
        hb2 = findobj(gca,'type','line');
        hb3 = findobj(gca,'type','Outliers');
        for iHB = 1:size(hb{iTask},2)
            set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
            set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
            set(hb2((iHB)+3:3:end),'MarkerFaceColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'MarkerFaceColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'Color',options.col_list{4-iHB})
        end
        hold on
        
        % Beeswarm
        x_val = [1 2 3];
        bee_bin_width = .1;
        bee_spread_width = .5;
        beePlot = plotSpread({dur_data_plot(grouping==1),dur_data_plot(grouping==2),dur_data_plot(grouping==3)},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.8 .8 .8]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
        hold on
        hbCurr = findobj(gca,'type','line');
        for iHB = 1:size(hb{iTask},2)
            set(hbCurr((iHB)+3:3:end),'color',options.col_list{4-iHB})
            set(hbCurr((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
        end
        
        if iTask == 1
            if options.normalize_plot == 0
                plot([0 4],[11 11],'--k');
            elseif options.normalize_plot == 1
                plot([0 4],log10([11 11]),'--k');
            end
        end
        
        % Plot significance
%         max_Dur = max([dur_data_plot(grouping==1); dur_data_plot(grouping==2); dur_data_plot(grouping==3)]);
        max_Dur = 50;
        if options.plot_stats == 1
            if iTask == 1
                % Plot 3-K-W
                text(1,max_Dur*.75,...
                    ['X2(' sprintf('%d',output.(task_names{iTask}).duration.kruskall_wallis_3_groups.table{2,3})  ') = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).duration.kruskall_wallis_3_groups.table{2,5}) ', p = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).duration.kruskall_wallis_3_groups.table{2,6})],...
                    'fontsize',statsFontSize);
            elseif iTask == 2
                % Plot anova results
                text(1,max_Dur*.85,...
                    ['F(' sprintf('%d',output.(task_names{iTask}).duration.anova_table{2,3}) ',' ...
                    sprintf('%d',output.(task_names{iTask}).duration.anova_table{3,3}) ') = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).duration.anova_table{3,6}) ', p = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).duration.anova_table{3,7})],...
                    'fontsize',statsFontSize);
                % Plot post hoc 2-K-W for controls vs patients
                text(1,max_Dur*.75,...
                    ['X2(' sprintf('%d',output.(task_names{iTask}).duration.kruskall_wallis_2_groups{2}.table{2,3})  ') = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).duration.kruskall_wallis_2_groups{2}.table{2,5}) ', p = ' ...
                    sprintf('%1.3f',output.(task_names{iTask}).duration.kruskall_wallis_2_groups{2}.table{2,6})],...
                    'fontsize',statsFontSize);
            end
        end
       
        title(type_labels{iTask},'fontsize',titleFontSize)
        box off
        if iTask == 1
            ylabel('Average Percept Duration (sec)','fontsize',axisTitleFontSize)
        end
        set(gca,'YScale','log')
        set(gca,'ylim',[1 60],'ytick',[1 2.5 5 10 25 50])
        set(gca,'XColor','k','YColor','k')
        
        set(gcf,'Units','inches')
        set(gcf,'Position',figSize.blockAve.figSize,'color','w')
    end
    
    
    %% Plot time course of responses and switch rate for 7T Methods
    titleFontSize = 12;
    axisTitleFontSize = 10;
    axisLabelFontSize = 8;
    if iTask == 2   % Only for illusory task
        if options.displayFigs
            figure()
            figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels
            figSize.switchRate.aspectRatio = [10.9849 9.2814];   % Aspect ratio
            figSize.switchRate.figSize = [0 0 ...
                3.54 6];   % Size/postion of fig
            set(gcf,'color','w')
            set(gca,'XColor','k','YColor','k')
            set(gcf,'units','inches')
            set(gcf,'Position', figSize.switchRate.figSize)
            subplot(1,2,iTask)
            addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
            
            %% Time course
            % Determine x axis
            xAxis = 0:.001:120;
            xAxisTicks = 0:15:120;
            %             directionSwitch = 1;
            %             for iI=1:length(xAxis)
            %                 if sum(round(SFMdata.flip_time{3,1,2,3},3) == xAxis(iI))
            %                     directionSwitch = 3-directionSwitch;
            %                 end
            %                 timeCourseData(iI) = directionSwitch;
            %             end
            %             yAxisLabels = {'CW','CCW'};
            
            % Plot
            subplot(2,1,1)
            imagesc(flipTime_data)
            %             plot(xAxis,timeCourseData,'linewidth',2)
            hold on
            set(gca,'fontsize',axisLabelFontSize)
            set(gca,'XTickLabel',xAxisTicks,...
                'XTick',xAxisTicks*1000,...
                'fontsize',axisLabelFontSize,...
                'YTick',[])
            title('Example Timecourses','fontsize',titleFontSize)
            box off
            %             ylabel('Behavioral Response','fontsize',axisTitleFontSize)
            ylabel('Participant','fontsize',axisTitleFontSize)
            xlabel('Time (s)','fontsize',axisTitleFontSize)
            set(gca,'XColor','k','YColor','k')
            cb1 = colorbar;
            cb1.Ticks = [1 2];
            cb1.TickLabels = {'CW' 'CCW'};
            
            %% Boxplots
            subplot(2,1,2)
            hb{iTask} = boxplot(Hz_data_plot(grouping==1),grouping(grouping==1),'widths',3.5);
            pause(0.5);
            set(gca,'XTick',1:3,'XTickLabel',[x_labels{1}(1) 'ontrols' x_labels{1}(2:end)],...
                'fontsize',axisLabelFontSize)
            set(hb{iTask},'linewidth',2)
            hb2 = findobj(gca,'type','line');
            hb3 = findobj(gca,'type','Outliers');
            for iHB = 1:size(hb{iTask},2)
                set(hb2((iHB):end),'color',options.col_list{iHB})
                set(hb2((iHB):end),'MarkerEdgeColor',options.col_list{iHB})
                set(hb3((iHB):end),'MarkerEdgeColor',options.col_list{iHB})
                set(hb3((iHB):end),'Color',options.col_list{iHB})
            end
            set(hb2(2),'LineWidth',3)
            hold on
            
            % Plot for use in legend
            hold on;
            plot([0 0],[0 0],'-','color',[0 1 0],'linewidth',3)
            plot([0 0],[0 0],'-','color',[0 1 0],'linewidth',2)
            plot([0 0],[0 0],'--','color',[0 1 0],'linewidth',2)
            
            % Beeswarm
            x_val = [1];
            bee_bin_width = .075;
            bee_spread_width = 5;
            beePlot = plotSpread({Hz_data_plot(grouping==1)},...
                'binWidth', bee_bin_width,...
                'distributionColors', {[.8 .8 .8]},...
                'xValues', x_val,...
                'spreadWidth', bee_spread_width);
            set(beePlot{1},'MarkerSize',10)
            hold on
            
            % Plot line at physical switch
            if iTask == 1
                if options.normalize == 0
                    plot([0 4],[0.09 0.09],'--k');
                elseif options.normalize == 1
                    plot([0 4],log10([0.09 0.09]),'--k');
                end
            end
            
            % Plot legend
            legend('Median','25%-75%','Range','Data','Location','NorthEast')
            
            % Plot significance
            
            title(type_labels{iTask},'fontsize',titleFontSize)
            box off
            ylabel('Switch rate (Hz)','fontsize',axisTitleFontSize)
            set(gca,'XColor','k','YColor','k')
            set(gca,'YLim',[0 .33])
            
            % Plot figure labels
            annotation('textbox',[0, .99, 0, 0],'string','A)','fontsize',titleFontSize,'fontweight','bold')
            annotation('textbox',[0, .5, 0, 0],'string','B)','fontsize',titleFontSize,'fontweight','bold')
            
            %% Export file
            % Resize figure
%             set(gcf,'paperunits','inches')
%             set(gcf, 'PaperPositionMode', 'manual');
%             set(gcf,'papersize',[3.25 5])
%             set(gcf,'paperposition',[0,0,3,6])
%             set(gcf, 'renderer', 'painters');
%             cd ./Figures/'Paper Figs'/  
%             print('PaperFig_SFM_Example.png','-dpng');
%             cd ../../
        end
    end
    
    
    %% Coefficient of variance
    % quantifies relationship between percept duration and variaibility
    % theory: more variable = neural noise, more regular = neural adaptation
    % see Robertson (2013), Shpiro (2009)
%     all_data_kw = [[nanmean(SFMdata.coeff_var(grp1_idx,1,iTask,:),4); ...
%         nanmean(SFMdata.coeff_var(grp1_idx,2,iTask,:),4)] ; ...
%         [nanmean(SFMdata.coeff_var(grp2_idx,1,iTask,:),4); ...
%         nanmean(SFMdata.coeff_var(grp2_idx,2,iTask,:),4)] ; ...
%         [nanmean(SFMdata.coeff_var(grp3_idx,1,iTask,:),4); ...
%         nanmean(SFMdata.coeff_var(grp3_idx,2,iTask,:),4)]];
%     grouping_kw = [ones(1,numel(grp1_idx)*2) 2*ones(1,numel(grp2_idx)*2) 3*ones(1,numel(grp3_idx)*2)];
    all_data_kw = [nanmean(SFMdata.coeff_var(grp1_idx,:,iTask,:),4); ...
        nanmean(SFMdata.coeff_var(grp2_idx,:,iTask,:),4); ...
        nanmean(SFMdata.coeff_var(grp3_idx,:,iTask,:),4)];
    grouping_kw = [ones(numel(grp1_idx), size(all_data_kw,2)) ; ...
        2*ones(numel(grp2_idx), size(all_data_kw,2)) ; ...
        3*ones(numel(grp3_idx), size(all_data_kw,2)) ];
    output.(task_names{iTask}).CV.all_data = all_data_kw;
    output.(task_names{iTask}).CV.all_group = grouping_kw;
    
    % Average data to plot
    CV_data = nanmean(all_data_kw,2);
    CV_grouping = nanmean(grouping_kw,2);

    % Kruskall-Wallis
    [p_full_CV_kw(iTask), output.(task_names{iTask}).CV.kruskall_wallis_3_groups.table , ...
        output.(task_names{iTask}).CV.kruskall_wallis_3_groups.stats ] = ...
        kruskalwallis(CV_data, CV_grouping, show_stats_fig);
    output.(task_names{iTask}).CV.kruskall_wallis_3_groups.p = p_full_CV_kw(iTask);

%     % ANOVA
%     all_data_anova = [nanmean(SFMdata.coeff_var(grp1_idx,:,iTask,:),4); ...
%         nanmean(SFMdata.coeff_var(grp2_idx,:,iTask,:),4); ...
%         nanmean(SFMdata.coeff_var(grp3_idx,:,iTask,:),4)];
%     
%     all_subj_anova = repmat([1:size(all_data_anova,1)]',[1 size(all_data_anova,2) size(all_data_anova,3)]);
%     all_group_anova = [ones(numel(grp1_idx), size(all_data_anova,2)) ; ...
%         2*ones(numel(grp2_idx), size(all_data_anova,2)) ; ...
%         3*ones(numel(grp3_idx), size(all_data_anova,2)) ];
%     
%     %     % No nested variables (?)
%     %     nest = zeros(2,2);
%     %     nest(1,2) = 1;
%     
%     [panova_CV, output.(task_names{iTask}).CV.anova_table] = anovan(all_data_anova(:),{all_subj_anova(:),...
%         all_group_anova(:)},'random',1,'continuous',[2],...
%         'model','full','varnames',{'subj','group'},...
%         'display',show_stats_fig);

    %% Plot CV data
    if options.plotAveCoV
        figure(o)
        figSize.switchRate.baseSize = get(0,'Screensize');   % Base size in pixels 
        figSize.switchRate.aspectRatio = [5.5 5.5];   % Aspect ratio
        figSize.switchRate.figSize = [0 0 ...
            figSize.switchRate.baseSize(3)*.75 ...
            (figSize.switchRate.baseSize(3)*.75)*...
            (figSize.switchRate.aspectRatio(2)/figSize.switchRate.aspectRatio(1))];   % Size/postion of fig
        subplot(1,2,iTask)
        addpath(genpath('/home/shaw-raid/matlab_tools/mpsCode/plotSpread'))
        
        % Boxplots
        hb{iTask} = boxplot(CV_data,grouping);
        pause(0.5)
        set(gca,'XTick',1:3,'XTickLabel',x_labels,'fontsize',15)
        set(hb{iTask},'linewidth',2)
        hb2 = findobj(gca,'type','line');
        hb3 = findobj(gca,'type','Outliers');
        for iHB = 1:size(hb{iTask},2)
            set(hb2((iHB)+3:3:end),'color',options.col_list{4-iHB})
            set(hb2((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
            set(hb2((iHB)+3:3:end),'MarkerFaceColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'MarkerFaceColor',options.col_list{4-iHB})
            set(hb3((iHB)+3:3:end),'Color',options.col_list{4-iHB})
        end
        hold on
        
        % Beeswarm
        x_val = [1 2 3];
        bee_bin_width = .1;
        bee_spread_width = .5;
        beePlot = plotSpread({CV_data(grouping==1),CV_data(grouping==2),CV_data(grouping==3)},...
            'binWidth', bee_bin_width,...
            'distributionColors', {[.8 .8 .8]},...
            'xValues', x_val,...
            'spreadWidth', bee_spread_width);
        set(beePlot{1},'MarkerSize',10)
        hold on
        hbCurr = findobj(gca,'type','line');
        for iHB = 1:size(hb{iTask},2)
            set(hbCurr((iHB)+3:3:end),'color',options.col_list{4-iHB})
            set(hbCurr((iHB)+3:3:end),'MarkerEdgeColor',options.col_list{4-iHB})
        end
%         if iTask == 1
%             if options.normalize_plot == 0
%                 plot([0 4],[11 11],'--k');
%             elseif options.normalize_plot == 1
%                 plot([0 4],log10([11 11]),'--k');
%             end
%         end
    
       
        title(type_labels{iTask},'fontsize',18)
        box off
        if iTask == 1
            ylabel('Average Cumulative Variance','fontsize',15)
        end
%         set(gca,'YScale','log')
%         set(gca,'ylim',[1 60],'ytick',[1 2.5 5 10 25 50])
        set(gcf,'color','w')
        set(gca,'XColor','k','YColor','k')
        set(gcf,'Position', figSize.switchRate.figSize)
    end
    
end

%% output
output.options = options;
output.date_run = datestr(now);

end