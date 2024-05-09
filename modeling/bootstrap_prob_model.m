function output = bootstrap_prob_model(options)
% usage: output = boostrap_prob_model(options)
%
% options = structure with lots of fields, all of which have default values
%           if not provided
% output = structure with bootstrap results
%
% Adapted from:
% Said and Heeger (2013) A model of binocular rivalry and cross-orientation
% suppression. PLOS Computational Biology.
%
% mps 20191026

%% Default Options 
    %Vanessa Morgan, 2/8/23

%This section creates a struct of default options to compare the options passed to
defaultOptions = [];
defaultOptions.nLayers = 3; %1st = near depth, 2 = far depth, 3 = depth invariant 
defaultOptions.normGain = 1;% scale normalization factor
defaultOptions.tau = 3250;% set time constant -- 3250 matches controls, 2100 matches SZ
defaultOptions.nBoot = 100;% # times to boostrap
defaultOptions.displayFigs = 1;% 1 = yes, 0 = no
defaultOptions.sigma = 0.5;% semisaturation constant
defaultOptions.sigma_opp = 0.9;% semisaturation constant for opponency cells
defaultOptions.bar_color = [0.5 1 0.5];% color for histogram
defaultOptions.attnGain = 1;% attention gain
defaultOptions.noiseAmp = 0.2;% noise amplitude
defaultOptions.noisefilter_t = 800;% set noise filtering in time

%% options
if ~exist('options','var') % parameters
    options = [];
end
if ~isfield(options,'nLayers')
    options.nLayers = defaultOptions.nLayers; % let's start with modeling 3 layers by default
    % where 1st = near depth, 2 = far depth, 3 = depth invariant 
end
if ~isfield(options,'normGain')
    options.normGain = defaultOptions.normGain; % scale normalization factor
end
if ~isfield(options,'tau')
    options.tau = defaultOptions.tau; % set time constant -- 3250 matches controls, 2100 matches SZ
end
if ~isfield(options,'nBoot')
    options.nBoot = defaultOptions.nBoot; % # times to boostrap
end
if ~isfield(options,'displayFigs')
    options.displayFigs = defaultOptions.displayFigs; % 1 = yes, 0 = no
end
if ~isfield(options,'sigma')
    options.sigma = defaultOptions.sigma; % semisaturation constant
end
if ~isfield(options,'sigma_opp')
    options.sigma_opp = defaultOptions.sigma_opp; % semisaturation constant for opponency cells
end
if ~isfield(options,'bar_color')
    options.bar_color = defaultOptions.bar_color; % color for histogram
end
if ~isfield(options,'attnGain')
    options.attnGain = defaultOptions.attnGain; % attention gain
end
if ~isfield(options,'noiseAmp')
    options.noiseAmp = defaultOptions.noiseAmp; % noise amplitude
end
if ~isfield(options,'noisefilter_t')
    options.noisefilter_t = defaultOptions.noisefilter_t; %  noise filtering in time
end

%% run boot
use_opt = options;
use_opt.displayFigs = 0; % don't show figures...

h_bar = waitbar(0,'Bootstrapping prob. model, please wait...');
for iBoot = 1:options.nBoot
    boot(iBoot) = run_prob_model(use_opt);
    waitbar(iBoot/options.nBoot, h_bar);
end
close(h_bar);


%% calculate means
for iBoot = 1:options.nBoot
    nReversals(iBoot) = boot(iBoot).nReversals;
    Hz(iBoot) = boot(iBoot).Hz;
    avg_dur(iBoot)= boot(iBoot).avg_dur;
    CV_dur(iBoot)= boot(iBoot).CV_dur;
end
output.reversals.mean = mean(nReversals);
output.reversals.median = median(nReversals);
output.reversals.sd = std(nReversals);
output.reversals.range = [min(nReversals) max(nReversals)];

output.Hz.mean = mean(Hz);
output.Hz.median = median(Hz);
output.Hz.sd = std(Hz);
output.Hz.range = [min(Hz) max(Hz)];

output.duration.mean = mean(avg_dur);
output.duration.median = median(avg_dur);
output.duration.sd = std(avg_dur);
output.duration.range = [min(avg_dur) max(avg_dur)];

output.CV.mean = mean(CV_dur);
output.CV.median = median(CV_dur);
output.CV.sd = std(CV_dur);
output.CV.range = [min(CV_dur) max(CV_dur)];

%% Comparing Options
    %Vanessa Morgan, 2/8/23

%Allows for a quick view of changed options from the output variable in
%the workspace

%Takes names of all fields of struct defaultOptions
optionNames = fieldnames(defaultOptions);
output.changedOptions = defaultOptions;

%For each struct field
for opt = 1:numel(optionNames)
    %Subtracts default option value from specified options
        %Any variable that doesn't change will list as 0
    output.changedOptions.(optionNames{opt}) = options.(optionNames{opt}) - defaultOptions.(optionNames{opt});
end

%% plot
if options.displayFigs
    figure;
    subplot(2,2,1);
    hold on;
    hist(nReversals)
    
    hb2 = findobj(gca,'type','Patch');
    hb2.FaceColor = options.bar_color;
    hb2.EdgeColor = [0 0 0];
    
    title('# Reversals')
    xlabel('N')
    ax = axis;
    text(0.1*ax(2), 0.95*ax(4), ['median = ' num2str(output.reversals.median, 3)],...
        'fontsize', 14)

    
    subplot(2,2,2);
    hold on;
    hist(Hz)
    
    hb2 = findobj(gca,'type','Patch');
    hb2.FaceColor = options.bar_color;
    hb2.EdgeColor = [0 0 0];
    
    title('Switch Rate')
    xlabel('Hz')
    ax = axis;
    text(0.1*ax(2), 0.95*ax(4), ['mean = ' num2str(output.Hz.mean, 3) ' Hz'], ...
        'fontsize', 14)

    
    subplot(2,2,3);
    hold on;
    hist(avg_dur)
    
    hb2 = findobj(gca,'type','Patch');
    hb2.FaceColor = options.bar_color;
    hb2.EdgeColor = [0 0 0];
    
    title('Avg. Duration')
    xlabel('s')
    ax = axis;
    text(0.1*ax(2), 0.95*ax(4), ['median = ' num2str(output.duration.median, 3) ' s'], ...
        'fontsize', 14)
    

    subplot(2,2,4);
    hold on;
    hist(CV_dur)
    
    hb2 = findobj(gca,'type','Patch');
    hb2.FaceColor = options.bar_color;
    hb2.EdgeColor = [0 0 0];
    
    title('CV Duration')
    xlabel('arb. units')
    ax = axis;
    text(0.1*ax(2), 0.95*ax(4), ['median = ' num2str(output.CV.median, 3) ' s'], ...
        'fontsize', 14)

    set(gcf,'color','w')
end
%% out
output.boot = boot;
output.options = options;
output.date_run = datestr(now);
end