function SFMdata = analyze_SFM_data(options)
% usage: SFMdata = analyze_SFM_data(options)
%
% mps c. Feb 2019

%%
if ~exist('options','var')
    options = [];
end
if ~isfield(options,'top_dir')
    options.top_dir = '/home/shaw-raid1/data/pHCP.git/subjectResponseFiles/';
end
if ~isfield(options,'dateCutoff')
    options.dateCutoff = 1;   % Only take participants after a given date
    options.dateCutoffVal = 20210815;   % Don't take ppt after this date
    options.dateCutoffValConverted = datenum(num2str(options.dateCutoffVal),'yyyymmdd');
end
if ~isfield(options,'excludeRedcap')
    options.excludeRedcap = 1; % exclude subjects based on control task performance, 0 = no, 1 = yes
    addpath(genpath('/home/shaw-raid1/matlab_tools/COP_analysis.git')) % Add path for redcap excludion function (in COP_analysis.git)
end
if ~isfield(options,'displayFigs')
    options.displayFigs = 0; % 1 = on, 0 = off
end

%% pull in data

all_pickles = dir(fullfile(options.top_dir,'P*SFM*.pickle'));
all_txt = dir(fullfile(options.top_dir,'P*SFM*.txt'));
missing_txt = [];
excluded_data = {};

if numel(all_pickles) ~= numel(all_txt)
for iP = 1:numel(all_pickles)
    find_me = 0;
    for iM = 1:numel(all_txt)
        if strcmp(all_pickles(iP).name(1:end-7),all_txt(iM).name(1:end-4)) || ...
            sum(strcmp(all_pickles(iP).name(1:end-7),excluded_data))
            find_me = 1;
            break
        end
    end
    if ~find_me
        missing_txt = str2mat(missing_txt, all_pickles(iP).name(1:end-7));
    end
end
end
sfmData.missing_txt = missing_txt;
if ~isempty(missing_txt)
    warning('Quitting, because I found .pickle files that weren''t converted to .txt');
    return
else sfmData = rmfield(sfmData,'missing_txt');
end
%% sort out subjects
for iS = 1:numel(all_txt)
    data_files{iS} = fullfile(options.top_dir,all_txt(iS).name);
    string_idx = regexp(data_files{iS},'P\d\d\d\d\d\d\d');
    subj_name{iS} = data_files{iS}(string_idx:string_idx+7);
    subj_number(iS) = str2num(data_files{iS}(string_idx+1:string_idx+7));
    
    datetime_idx = regexp(data_files{iS},'\d\d\d\d\d\d\d\d');
    date_num(iS) = datenum(data_files{iS}(datetime_idx:datetime_idx+7),'yyyymmdd');

    type_idx = regexp(data_files{iS},'SFM_type') + length('SFM_type');
    SFM_type{iS} = data_files{iS}(type_idx);
end

% Get redcap exclusion list KWK-20210901
if options.excludeRedcap == 1
    exclusion_opts = [];
    exclusion_opts.subj_number = subj_number'; % numeric, not a string
    exclusion_opts.date_number = date_num'; % per the matlab function datenum, NOT a string!!
    exclusion_opts.overwrite_exclusion_csv = 0;
    
    redcap_exclusion_output = redcap_exclusion(exclusion_opts); % this function lives in the COP_analysis.git repo
elseif options.excludeRedcap == 0
    redcap_exclusion_output.excluded_binary = zeros([length(subj_number) 1]);
end

%% Filter by control performance
% KWK - 20200415
if options.excludeTypeA
    counter=0;
    for iI=1:length(subj_number)
        % Exclude subjects based on redcap
        if redcap_exclusion_output.excluded_binary(iI) ~= 1
            % Exclude subjects run after the cutoff date
            if options.dateCutoff & ...
                    options.dateCutoffValConverted >= date_num(iI)
                if ismember(subj_number(iI),options.illusoryCutoff(:,1))
                    if options.illusoryCutoff(find(subj_number(iI)==options.illusoryCutoff(:,1)),4)==1
                        counter=counter+1;
                        new_data_files{counter} = data_files{iI};
                        new_subj_name{counter} = subj_name{iI};
                        new_subj_number(counter) = subj_number(iI);
                        new_date_num(counter) = date_num(iI);
                        new_SFM_type{counter} = SFM_type{iI};
                    end
                else
                    new_data_files{counter} = data_files{iI};
                    new_subj_name{counter} = subj_name{iI};
                    new_subj_number(counter) = subj_number(iI);
                    new_date_num(counter) = date_num(iI);
                    new_SFM_type{counter} = SFM_type{iI};
                end
            elseif options.dateCutoff == 0
                if ismember(subj_number(iI),options.illusoryCutoff(:,1))
                    if options.illusoryCutoff(find(subj_number(iI)==options.illusoryCutoff(:,1)),4)==1
                        counter=counter+1;
                        new_data_files{counter} = data_files{iI};
                        new_subj_name{counter} = subj_name{iI};
                        new_subj_number(counter) = subj_number(iI);
                        new_date_num(counter) = date_num(iI);
                        new_SFM_type{counter} = SFM_type{iI};
                    end
                else
                    new_data_files{counter} = data_files{iI};
                    new_subj_name{counter} = subj_name{iI};
                    new_subj_number(counter) = subj_number(iI);
                    new_date_num(counter) = date_num(iI);
                    new_SFM_type{counter} = SFM_type{iI};
                end
            end
        end
    end
    % Reset arrays
    clear data_files subj_name subj_number date_num SFM_type
    data_files = new_data_files;
    subj_name = new_subj_name;
    subj_number = new_subj_number;
    date_num = new_date_num;
    SFM_type = new_SFM_type;
end


%% quantify metrics
n_blocks = [1 5]; % A & B respectively
max_blocks = max(n_blocks);
block_dur = 120;
run_types = {'A','B'};
n_run_types = numel(run_types); % task = A or B
run_indexing = 1:n_run_types;
max_num_sessions = 2; % A and B, or B and Z

[unique_subj_num, unique_idx, subj_idx] = unique(subj_number);

flip_time = cell(numel(unique_subj_num), max_num_sessions, n_run_types, max_blocks);
percept_dur = cell(numel(unique_subj_num), max_num_sessions, n_run_types, max_blocks);
n_flips = nan(numel(unique_subj_num), max_num_sessions, n_run_types, max_blocks);
Hz_flips = nan(numel(unique_subj_num), max_num_sessions, n_run_types, max_blocks);
mean_percept_dur = nan(numel(unique_subj_num), max_num_sessions, n_run_types, max_blocks);
coeff_var = nan(numel(unique_subj_num), max_num_sessions, n_run_types, max_blocks);

SFMdata.subj_number = nan(numel(unique_subj_num),1);
SFMdata.date_num = nan(numel(unique_subj_num),max_num_sessions);

for iS = 1:numel(unique_subj_num)
    subj_files = find(subj_idx == subj_idx(unique_idx(iS)));
    
    all_dates = date_num(subj_files);
    
    [unique_dates, unique_d_idx, date_idx] = unique(all_dates);
    
    if numel(unique_dates) > max_num_sessions
        error(['more than ' num2str(max_num_sessions) ' dates for this subject??'])
    end
    
    clear use_run_type_idx
    for iF = 1:numel(subj_files)
        use_date_idx = date_idx(iF);
        use_run_type_idx = run_indexing(strcmp(SFM_type{subj_files(iF)} , run_types));
    
        data = tdfread(data_files{subj_files(iF)},'tab');
        
        for iB = 1:n_blocks(use_run_type_idx)
            resp_idx = find(data.block == iB-1);
            all_resp = data.key(resp_idx);
            flip_resp_idx = resp_idx( find(all_resp(2:end) ~= all_resp(1:end-1)) + 1 ); 
            % add 1 because we want the 2nd of the pair, which differs from the 1st...
            % also need to index into resp_idx, so we're looking at the
            % corect block...
            
            flip_time{iS, use_date_idx, use_run_type_idx, iB} = ...
                data.time(flip_resp_idx);
            percept_dur{iS, use_date_idx, use_run_type_idx, iB} = ...
                [flip_time{iS, use_date_idx, use_run_type_idx, iB} - ...
                [0 ; flip_time{iS, use_date_idx, use_run_type_idx, iB}(1:end-1)]];
            n_flips(iS, use_date_idx, use_run_type_idx, iB) = ...
                numel(flip_time{iS, use_date_idx, use_run_type_idx, iB});
            Hz_flips(iS, use_date_idx, use_run_type_idx, iB) = ...
                n_flips(iS, use_date_idx, use_run_type_idx, iB)/block_dur;
            mean_percept_dur(iS, use_date_idx, use_run_type_idx, iB) = ...
                mean(percept_dur{iS, use_date_idx, use_run_type_idx, iB});
            coeff_var(iS, use_date_idx, use_run_type_idx, iB) = ...
                std(percept_dur{iS, use_date_idx, use_run_type_idx, iB}) / ...
                mean_percept_dur(iS, use_date_idx, use_run_type_idx, iB); % mps 20200519
        end
        
    end
    SFMdata.subj_number(iS,1) = subj_number(unique_idx(iS));
    SFMdata.date_num(iS,1:numel(unique_dates)) = unique_dates;
end

%% Import diagnostic/demographic info - KWK 20220123
% Must implement group level segmenting (probands, relatives, and
% controls).
options.target_file = '/home/shaw-raid1/data/7T/demographics/PHCP7TfMRIDemo.csv';
addpath(genpath('/home/shaw-raid1/matlab_tools/mpsCode/'))   % Add path for MPS function read_in_demog_data
demogDataHolder = read_in_demog_data(options);
for iI=1:length(SFMdata.subj_number)
    subjIdx = find((ismember(demogDataHolder.Record_ID,...
        ['P' num2str(SFMdata.subj_number(iI))])));
    if isempty(subjIdx)
        SFMdata.dxId(iI) = NaN;
    else
        SFMdata.dxId(iI) = demogDataHolder.Dx_code(subjIdx);
    end
end
options = rmfield(options,'target_file');
SFMdata.dxId = SFMdata.dxId';

%% output
SFMdata.n_flips = n_flips;
SFMdata.Hz_flips = Hz_flips;
SFMdata.mean_percept_dur = mean_percept_dur;
SFMdata.flip_time = flip_time;
SFMdata.percept_dur = percept_dur;
SFMdata.coeff_var = coeff_var; % mps 20200519
SFMdata.dim_key = 'subj # x session (1st or 2nd) x run type (type A or type B, i.e., control or main task) x block #';
end