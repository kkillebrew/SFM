function [output] = combine_sfm_and_mrs(options)
% usage: [output] = combine_sfm_and_mrs(options)
%
% mps 20190910
%%
addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode/'));
addpath(genpath('/home/shaw-raid1/data/psychophysics/SFM.git'));
%% opt
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 1; % 1 = yes, 0 = no
end
if ~isfield(options,'displayControlFigs')   % Display figures in the control task analysis
    options.displayControlFigs = 0; % 1 = on, 0 = off
end
if ~isfield(options,'normalize')
    options.normalize = 1;   % 0=no norm; 1=log10
end
if ~isfield(options,'mrs_file')
    error('No options.mrs_file provided!')
    % e.g., options.mrs_file = '/home/shaw-raid1/data/MRS/processed_data/20220823_phcp_OCC_192subj_H2O_scaled.csv';
end
if ~isfield(options,'mrs_n_col')
    options.mrs_n_col = 504; % LCM default, if using Gosia's notes = 446
end
if ~isfield(options,'mrs_header_lines')
    options.mrs_header_lines = 6; % LCM default, if using Gosia's notes = 6
end
mrs_opt.target_file = options.mrs_file;
mrs_opt.n_col = options.mrs_n_col;
mrs_opt.header_lines = options.mrs_header_lines;

if ~isfield(options,'mrs_struct')
    options.mrs_struct = read_in_LCModel_results(mrs_opt);
end

if ~isfield(options,'sfm_struct')
    warning('Running analyze_SFM_data.m to generate options.sfm_struct...');
    options.excludeTypeA = 1;
    [dataAve] = SFM_Behav_Control_Group_Analysis;
    illusoryCutoff = dataAve.illusoryCutoff;
    [~,illusoryCutoffIdx] = unique(illusoryCutoff(:,1));
    illusoryCutoff = illusoryCutoff(illusoryCutoffIdx,:);
    options.illusoryCutoff = illusoryCutoff;
    options.sfm_struct = analyze_SFM_data(options);
end

if ~isfield(options,'which_corr')
    options.which_corr = 'Spearman'; % thresholds look like the are not equal variance between groups, so let's go with this
end
if ~isfield(options,'toss_subj_num_date')
    options.toss_subj_num_date = []; % these subjects look like the have bad MRS data on 20190331
%     options.toss_subj_num_date = [6004213 datenum('20190116','yyyymmdd')]; % these subjects look like the have bad MRS data on 20190331
end
if ~isfield(options,'which_group')
    options.which_group = 'all'; % options are all, controls, relatives, patients
end
ctrl_names = {'controls','control','ctrl','ctrls','c','ctr','ctl'};
rel_names = {'relatives','relative','rel','r','rels'};
pat_names = {'patients','patient','pat','p','pats','pts','ot'};
all_names{1} = ctrl_names;
all_names{2} = rel_names;
all_names{3} = pat_names;
if strcmp(options.which_group,'all')
    which_group_idx = 0;
elseif sum(strcmp(options.which_group,ctrl_names))
    which_group_idx = 1;
elseif sum(strcmp(options.which_group,rel_names))
    which_group_idx = 2;
elseif sum(strcmp(options.which_group,pat_names))
    which_group_idx = 3;
else
    error(['options.which_group not recognized: ' options.which_group]);
end
if ~isfield(options,'which_metab')
        options.which_metab = {'Glu','GABA','Gln'};
%         options.which_metab = {'Glu','GABA','Gln','GSH','NAA','NAAG'};
    %     options.which_metab = {'Glu','GABA','Gln','GSH','NAA','NAAG','Asp','Asc','MacY','Glc'};
%     options.which_metab = {'Glu','GABA','Gln','MacY','NAA','NAAG','Glc'};
    warning('options.which_metab not specified, assuming you want to look at only Glu and GABA...')
end
if ~isfield(options,'which_sfm_metric')
    options.which_sfm_metric = 'Hz_flips';
%     options.which_sfm_metric = 'percept_dur';
end

% if ~isfield(options,'toss_sfm_outliers')
%     options.toss_sfm_outliers = 0;
% end

%% get mrs file info
for iFile = 1:numel(options.mrs_struct.row_name)-4 % skip last 2, mean and sd/mean
    name_idx = regexp(options.mrs_struct.row_name{iFile},'P\d\d\d\d\d\d\d');
    options.mrs_struct.subj_number(iFile,1) = str2num(options.mrs_struct.row_name{iFile}...
        (name_idx+1:name_idx+7));
    date_idx = regexp(options.mrs_struct.row_name{iFile},'\d\d\d\d\d\d\d\d');
    options.mrs_struct.date_number(iFile,1) = datenum(options.mrs_struct.row_name{iFile}...
        (date_idx:date_idx+7),'yyyymmdd');
end
%% find subjects & dates that overlap
clear corr_sfm_hz

% start with MRS, because all the subjects should have psychophysics...?
MRS_subj_date = [options.mrs_struct.subj_number options.mrs_struct.date_number];

% Grab subjects based on which_group_idx (take all participants or a subset)
toss_idx = []; wrong_grp_idx = [];
if ~isempty(options.toss_subj_num_date) || which_group_idx
    for iSubj = 1:size(MRS_subj_date,1)

        if which_group_idx == 1
            if MRS_subj_date(iSubj,1) >= 2000000
                wrong_grp_idx = [wrong_grp_idx ; iSubj];
            end
        elseif which_group_idx == 2
            if MRS_subj_date(iSubj,1) < 2000000 || ...
                    MRS_subj_date(iSubj,1) >= 6000000 
                wrong_grp_idx = [wrong_grp_idx ; iSubj];
            end
        elseif which_group_idx == 3
            if MRS_subj_date(iSubj,1) < 6000000;
                wrong_grp_idx = [wrong_grp_idx ; iSubj];
            end
        end
        
        if ~isempty(options.toss_subj_num_date)
            if which_group_idx && ~isempty(wrong_grp_idx)
                if wrong_grp_idx(end) ~= iSubj % if we haven't already tossed this one...
                    for iToss = 1:size(options.toss_subj_num_date,1)
                        if sum(MRS_subj_date(iSubj,:) == options.toss_subj_num_date(iToss,:)) == 2
                            % both subject num and date are the same, so toss
                            toss_idx = [toss_idx ; iSubj];
                        end
                    end
                end
            else
                for iToss = 1:size(options.toss_subj_num_date,1)
                    if sum(MRS_subj_date(iSubj,:) == options.toss_subj_num_date(iToss,:)) == 2
                        % both subject num and date are the same, so toss
                        toss_idx = [toss_idx ; iSubj];
                    end
                end
            end
        end
    end
end

start_date = datenum('20170714','yyyymmdd')-1; % this was the day we ran the first scan, minus 1, 
% use it to cut down matrix size for indexing
max_date = max([reshape(options.sfm_struct.date_num(:,1), ...
    numel(options.sfm_struct.date_num(:,1)),1) ; options.mrs_struct.date_number]);
all_sfm = zeros(7000000, max_date - start_date); % looks like using zeros takes up a lot less memory then NaNs... 16GB YIKES
%%%% this approach is not perfect, because it treats all scans equally -
%%%% doesn't take repeated scans within subjects into account!!

rep_sfm_subj = repmat(options.sfm_struct.subj_number,[1 2]);
sfm_idx = ~isnan(options.sfm_struct.date_num);
all_sfm_subj = rep_sfm_subj(sfm_idx);
all_sfm_date = options.sfm_struct.date_num(sfm_idx);
illusory_idx = 2; % index for illusory stimulus blocks (i.e., block B)
sfm_data = squeeze(nanmean(options.sfm_struct.(options.which_sfm_metric)...
    (:,:,illusory_idx,:),4));
sfm_data = sfm_data(sfm_idx);

thresh_idx = sub2ind(size(all_sfm),all_sfm_subj,...
    all_sfm_date - start_date);
all_sfm(thresh_idx) = sfm_data;

mrs_idx = sub2ind(size(all_sfm),options.mrs_struct.subj_number,...
    options.mrs_struct.date_number - start_date); % need to use the same matrix size, but this time with mrs indices
corr_sfm_hz = all_sfm(mrs_idx);
clear all_sfm


% if the subject is included in the MRS data set, but not the SFM for some
% reason, then they will have a zero for corr_sfm_hz - to fix, replace
% all zeros with NaN
idx_no_thresh = corr_sfm_hz == 0;
output.subj_number_date_no_thresh = MRS_subj_date(idx_no_thresh,:);
warning([num2str(size(output.subj_number_date_no_thresh,1)) ' subjects with missing SFM data...']);
% if looking at OCC, should be the first 10 subjects before 20170906 who
% didn't get CSS psychophysics...

corr_sfm_hz(corr_sfm_hz == 0) = NaN;

%% toss outliers?
% if options.toss_css_outliers
%     find_out = abs(corr_css_params - repmat(nanmean(corr_css_params,1),...
%         [size(corr_css_params,1) 1] )) > repmat(3*nanstd(corr_css_params,0,1),...
%         [size(corr_css_params,1) 1] );
%     corr_css_params(find_out) = NaN;
% 
%     find_out = abs(corr_sfm_hz - nanmean(corr_sfm_hz,1) ) > 3*nanstd(corr_sfm_hz,0,1);
%     corr_sfm_hz(find_out) = NaN;
% end
%% calculate correlations and plot
for iM = 1:numel(options.which_metab)
    
    eval(['corr_metab = options.mrs_struct.' options.which_metab{iM}...
        '(1:numel(options.mrs_struct.subj_number));']);
    
    if ~isempty(wrong_grp_idx)
        warning(['Including only ' num2str(size(MRS_subj_date,1) - ...
            numel(wrong_grp_idx)) ' ' all_names{which_group_idx}{1} ...
            ', as requested...']);
        corr_metab(wrong_grp_idx) = NaN;
    end
    if ~isempty(toss_idx)
        warning(['Tossing ' num2str(numel(toss_idx)) ' data sets, as requested...']);
        corr_metab(toss_idx) = NaN;
    end

    % first correlate thresholds
    throw_out_nans_metab = isnan(corr_metab) | isnan(corr_sfm_hz);
    corr_sfm_hz_metab = corr_sfm_hz;
    corr_sfm_hz_metab(throw_out_nans_metab) = [];
    corr_metab_toss = corr_metab;
    corr_metab_toss(throw_out_nans_metab) = [];
    clear metab
    
    % Determine group index for colors
    comb_MRS_SFM_subj_idx = MRS_subj_date(~throw_out_nans_metab,1);
    group_color_idx(comb_MRS_SFM_subj_idx < 2000000) = 1;
    group_color_idx(comb_MRS_SFM_subj_idx >= 2000000 & comb_MRS_SFM_subj_idx < 6000000 ) = 2;
    group_color_idx(comb_MRS_SFM_subj_idx >= 6000000) = 3;
    
    % Make grouping array for unique individuals for seperate group correlations
    group_idx(unique(comb_MRS_SFM_subj_idx) < 2000000) = 1;
    group_idx(unique(comb_MRS_SFM_subj_idx) >= 2000000 & unique(comb_MRS_SFM_subj_idx) < 6000000 ) = 2;
    group_idx(unique(comb_MRS_SFM_subj_idx) >= 6000000) = 3;
    
    % Log10 normalize the SFM data
    plot_sfm_hz_metab = corr_sfm_hz_metab;
    plot_metab_toss = corr_metab_toss;
    if options.normalize == 1
        corr_sfm_hz_metab = log10(corr_sfm_hz_metab);
    end
    
    % Set array for correct values to use in LME vs corr (include repeats
    % in LME and not in spearman corr)
    lme_sfm_hz_metab = corr_sfm_hz_metab;
    lme_metab_toss = corr_metab_toss;
    
    % Find unique subjects to use with corr
    [comb_MRS_SFM_subj_idx_unique,unique_idx,~] = unique(comb_MRS_SFM_subj_idx);
    corr_sfm_hz_metab = corr_sfm_hz_metab(unique_idx);
    corr_metab_toss = corr_metab_toss(unique_idx);
    
    % Do correlation
    [metab.corr_hz.r, metab.corr_hz.p] = corr(corr_metab_toss, corr_sfm_hz_metab, ...
        'type',options.which_corr);
    metab.corr_hz.df = numel(corr_metab_toss)-2;
    
    % Do correlations for seperate groups
    % Controls
    [metab.corr_hz_group(1).r, metab.corr_hz_group(1).p] = corr(corr_metab_toss(group_idx==1), corr_sfm_hz_metab(group_idx==1), ...
        'type',options.which_corr);
    metab.corr_hz_group(1).df = numel(corr_metab_toss)-2;
    % Relatives
    [metab.corr_hz_group(2).r, metab.corr_hz_group(2).p] = corr(corr_metab_toss(group_idx==2), corr_sfm_hz_metab(group_idx==2), ...
        'type',options.which_corr);
    metab.corr_hz_group(2).df = numel(corr_metab_toss)-2;
    % Probands
    [metab.corr_hz_group(3).r, metab.corr_hz_group(3).p] = corr(corr_metab_toss(group_idx==3), corr_sfm_hz_metab(group_idx==3), ...
        'type',options.which_corr);
    metab.corr_hz_group(3).df = numel(corr_metab_toss)-2;
    
    % Do LME model
    % Make lme table
    metab.lme_hz.lmeTable = table(...
        lme_sfm_hz_metab,comb_MRS_SFM_subj_idx,group_color_idx',lme_metab_toss,...
        'VariableNames',{'Switch','Subj','Group',options.which_metab{iM}});
    
    % Turn the categorical variables into 'categorical' type
    metab.lme_hz.lmeTable.Subj = categorical(...
        metab.lme_hz.lmeTable.Subj);
    metab.lme_hz.lmeTable.Group = categorical(...
        metab.lme_hz.lmeTable.Group);
    
    metab.lme_hz.formula =...
        ['Switch ~ ' options.which_metab{iM} ' + (1 | Subj)'];
    metab.lme_hz.output =...
       fitlme(metab.lme_hz.lmeTable,...
       metab.lme_hz.formula);
   
    if options.displayFigs
        figure; hold on
        % Set figure size
        figSize.diag.baseSize = get(0,'Screensize');   % Base size in pixels
        figSize.diag.aspectRatio = [3 3];   % Aspect ratio
        figSize.diag.figSize = [0 0 ...
            (figSize.diag.baseSize(4)*...
            (figSize.diag.aspectRatio(1)/figSize.diag.aspectRatio(2)))+40 ...  % Make it just a bit bigger to fit the title
            figSize.diag.baseSize(4)];   % Size/postion of fig   % Size/postion of fig
        set(gcf,'Position',figSize.diag.figSize,'Units','inches')
        
        [poly_fit] = polyfit(plot_metab_toss, ...
            plot_sfm_hz_metab, 1);
        
        fit_x = [min(plot_metab_toss) max(plot_metab_toss)];
        fit_y = poly_fit(1).*fit_x + poly_fit(2);
        y_range = [min(plot_sfm_hz_metab) max(plot_sfm_hz_metab)];
        plot(fit_x,fit_y,'k-','linewidth',2)
        
        % Controls
        plot(plot_metab_toss(group_color_idx==1),...
            plot_sfm_hz_metab(group_color_idx==1),'go',...
            'linewidth',2,'MarkerSize',6,'MarkerFaceColor','w')
        % Relatives
        plot(plot_metab_toss(group_color_idx==2),...
            plot_sfm_hz_metab(group_color_idx==2),'bo',...
            'linewidth',2,'MarkerSize',6,'MarkerFaceColor','w')
        % Probands
        plot(plot_metab_toss(group_color_idx==3),...
            plot_sfm_hz_metab(group_color_idx==3),'ro',...
            'linewidth',2,'MarkerSize',6,'MarkerFaceColor','w')
        xlabel([options.which_metab{iM} ' (mM)'],'color','k')
        if strcmp(options.which_sfm_metric , 'Hz_flips')
            use_y_label = 'Switch rate (Hz)';
        elseif strcmp(options.which_sfm_metric , 'mean_percept_dur')
            use_y_label = 'Mean percept duration (s)';
        end
        ylabel(use_y_label,'color','k')
%         title(['SFM ' use_y_label ' vs. ' options.which_metab{iM}])
        set(gca,'YScale','log')
        set(gca,'ytick',[0.005 0.01 0.025 0.05 0.1 0.25 0.5])
        
        % Plot corr stats
        if options.plot_stats
            range_metab = max(plot_metab_toss) - min(plot_metab_toss);
            range_thresh = max(plot_sfm_hz_metab) - min(plot_sfm_hz_metab);
            text(max(plot_metab_toss)-range_metab*.2,max(plot_sfm_hz_metab)-range_thresh*.1,...
                ['n = ' num2str(numel(plot_metab_toss))],'fontsize',18)
            text(max(plot_metab_toss)-range_metab*.2,max(plot_sfm_hz_metab)-range_thresh*.3,...
                ['r = ' num2str(round(100*metab.corr_hz.r)/100)],'fontsize',18)
            text(max(plot_metab_toss)-range_metab*.2,max(plot_sfm_hz_metab)-range_thresh*.4,...
                ['p = ' num2str(round(100*metab.corr_hz.p)/100)],'fontsize',18)
            % Plot LME output
            lmeOutputIdx = strcmp(metab.lme_hz.output.Coefficients.Name,options.which_metab{iM});
            text(max(plot_metab_toss)-range_metab*.2,max(plot_sfm_hz_metab)-range_thresh*.55,...
                ['t(', sprintf('%d',metab.lme_hz.output.Coefficients.DF(1)),') = ' ...
                sprintf('%1.3f',metab.lme_hz.output.Coefficients.tStat(lmeOutputIdx)) ...
                ', p=' sprintf('%1.3f',metab.lme_hz.output.Coefficients.pValue(lmeOutputIdx))],'fontsize',15)
        end
        
        set(gcf,'color','w')
        set(gca,'FontSize',18,'XColor','k','YColor','k')
        x_span = fit_x(2) - fit_x(1);
        y_span = y_range(2) - y_range(1);
        axis([fit_x(1)-0.1*x_span fit_x(2)+0.1*x_span ...
            y_range(1)-0.1*y_span y_range(2)+0.1*y_span])  
    end
    
    eval(['output.' options.which_metab{iM} ' = metab;']);
end
%% out
output.options = options;
end