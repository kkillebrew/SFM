function output = run_prob_model(options)
% usage: output = run_prob_model(options)
%
% options = structure with lots of fields, all of which have default values
%           if not provided
% output = structure with model results
%
% Adapted from:
% Said and Heeger (2013) A model of binocular rivalry and cross-orientation
% suppression. PLOS Computational Biology.
%
% mps 20191024

%% Default Options
    %Vanessa Morgan 2/8/23

%Create a struct of default options to compare the options passed to
defaultOptions =[];
defaultOptions.displayFigs = 1;%1 = on, 0 = off
defaultOptions.nLayers = 3;%1st = near depth, 2 = far depth, 3 = depth invariant 
defaultOptions.normGain = 1;% scale normalization factor
defaultOptions.tau = 3250;% set time constant
defaultOptions.noiseAmp = 0.2;% arb. units
defaultOptions.noisefilter_t = 800;% set noise filtering in time
defaultOptions.sigma = 0.5;% semisaturation constant
defaultOptions.sigma_opp = 0.9;% semisaturation constant for opponency cells
defaultOptions.attnGain = 1;% scale attention factor

%% options
if ~exist('options','var') % parameters
    options = [];
end
if ~isfield(options,'displayFigs')
    options.displayFigs = defaultOptions.displayFigs; % 1 = on, 0 = off
end
if ~isfield(options,'nLayers')
    options.nLayers = defaultOptions.nLayers; % let's start with modeling 3 layers by default
    % where 1st = near depth, 2 = far depth, 3 = depth invariant 
end
if ~isfield(options,'normGain')
    options.normGain = defaultOptions.normGain; % scale normalization factor
end
if ~isfield(options,'tau')
    options.tau = defaultOptions.tau; % set time constant
end
if ~isfield(options,'noiseAmp')
    options.noiseAmp = defaultOptions.noiseAmp; % arb. units
end
if ~isfield(options,'noisefilter_t')
    options.noisefilter_t = defaultOptions.noisefilter_t; % set noise filtering in time
end
if ~isfield(options,'sigma')
    options.sigma = defaultOptions.sigma; % semisaturation constant
end
if ~isfield(options,'sigma_opp')
    options.sigma_opp = defaultOptions.sigma_opp; % semisaturation constant for opponency cells
end
if ~isfield(options,'attnGain')
    options.attnGain = defaultOptions.attnGain; % scale attention factor
end

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
%% set up params
c = 1; % stimulus contrast
input_amp_A = {[c c]}; % stimuli (input amplitude) for: SFM B
input_amp_B = {[c c]}; % stimuli for: SFM B

condnames =  {'SFM Type B (Illusory)'};
layernames = {'Near-Selective', 'Far-Selective', 'Depth Summation', ...
              'N-F Opponency', 'F-N Opponency'};

p.sigma         = options.sigma;         % semisaturation constant
p.sigma_opp     = options.sigma_opp;     % semisaturation constant for opponency cells
p.tau           = options.tau;           % time constant (ms)
p.dt            = 50;                    % time-step (ms)
p.sim_duration  = 120000;                % duration (ms)
p.noisefilter_t = options.noisefilter_t; % (ms)
p.noiseAmp      = options.noiseAmp;      % arb. units
p.nLayers       = options.nLayers;       % set to 3 for conventional model, 
                                         % 5 for opponency model
p.n_time_points = p.sim_duration/p.dt+1;
p.tlist         = 0:p.dt:p.sim_duration; % in ms
p.normGain      = options.normGain;      % scaling factor for normalization denominator
p.attnGain      = options.attnGain;      % scaling factor for attention

% Initializing time-courses for neuron (d)rives, (r)esponses, and (n)oise.
% Each neuron is tuned to either stimulus direction A or B.
for iLayer=1:5 % go through maximum possible layers. This way, if there are
               % < 5 layers, the feedback can be zero.
    p.dA{iLayer}   = zeros(1,p.n_time_points);
    p.dB{iLayer}   = zeros(1,p.n_time_points);
    p.rA{iLayer}   = zeros(1,p.n_time_points);
    p.rB{iLayer}   = zeros(1,p.n_time_points);
    p.nA{iLayer}   = makeNoise(p);
    p.nB{iLayer}   = makeNoise(p);
end

%% run model
% loop through conditions
for iCond = 1:size(input_amp_A,2)
    % stimulus inputs to depth-selective layers
    for iLayer = 1:2
        p.iA{iLayer} = input_amp_A{iCond}(iLayer)*ones(1,p.n_time_points);
        p.iB{iLayer} = input_amp_B{iCond}(iLayer)*ones(1,p.n_time_points);
    end
    
    % run the model
    p = prob_model(p);
    
    % compute winner-take-all index from summation layer
    wta(iCond) = nanmean(abs(p.rA{3}-p.rB{3})./(p.rA{3}+p.rB{3}));
    
    % find the reversals in the summation layer
    resp_diff = p.rA{3} - p.rB{3};
    CW_to_CCW_idx = find(resp_diff(2:end) > 0 & resp_diff(1:end-1) <= 0)+1;
    CCW_to_CW_idx = find(resp_diff(2:end) < 0 & resp_diff(1:end-1) >= 0)+1;
    switch_idx = find(resp_diff(2:end) > 0 & resp_diff(1:end-1) <= 0 | ...
        resp_diff(2:end) < 0 & resp_diff(1:end-1) >= 0)+1;
    if resp_diff(1) > 0 % resp_diff(1) == 0, so this doesn't matter, but
                        % let's be thorough...
        CW_to_CCW_idx = [1 CW_to_CCW_idx];
    elseif resp_diff(1) < 0
        CCW_to_CW_idx = [1 CCW_to_CW_idx];
    end
    CW_to_CCW_t = p.tlist(CW_to_CCW_idx);
    CCW_to_CW_t = p.tlist(CCW_to_CW_idx);
    
    % figure out dominance durations
    switch_t = p.tlist(switch_idx);
    switch_dur = [switch_t(2:end) - switch_t(1:end-1) p.sim_duration - ...
        switch_t(end)];
    
    % calculate # of reversals, frequency, avg. dominance duration
    output.nReversals = numel(switch_t);
    output.Hz = output.nReversals / (p.sim_duration / 1000); % convert from ms to s
    output.avg_dur = mean(switch_dur) / 1000; % in s
    output.CV_dur = std(switch_dur) / mean(switch_dur); % arb. units
    
    %If either CCW or CW == 0, then handle the exception

    % figure out indices for plotting patches
    if CW_to_CCW_idx(1) > CCW_to_CW_idx(1) % starts on Clockwise
        if CW_to_CCW_idx(end) > CCW_to_CW_idx(end) % end on Counter CW
            CW_to_CCW_idx_plotL = switch_idx(2:2:end);
            CCW_to_CW_idx_plotL = [switch_idx(3:2:end) numel(p.tlist)];
            CW_to_CCW_idx_plotR = switch_idx(2:2:end);
            CCW_to_CW_idx_plotR = switch_idx(1:2:end);
        else           
            CW_to_CCW_idx_plotL = switch_idx(2:2:end);
            CCW_to_CW_idx_plotL = switch_idx(3:2:end);
            CW_to_CCW_idx_plotR = [switch_idx(2:2:end) numel(p.tlist)];
            CCW_to_CW_idx_plotR = switch_idx(1:2:end);
        end
    else
        if CW_to_CCW_idx(end) > CCW_to_CW_idx(end) % end on Counter CW
            CW_to_CCW_idx_plotL = switch_idx(1:2:end);
            CCW_to_CW_idx_plotL = [switch_idx(2:2:end) numel(p.tlist)];
            CW_to_CCW_idx_plotR = switch_idx(2:2:end);
            CCW_to_CW_idx_plotR = switch_idx(3:2:end);
        else
            CW_to_CCW_idx_plotL = switch_idx(1:2:end);
            CCW_to_CW_idx_plotL = switch_idx(2:2:end);
            CW_to_CCW_idx_plotR = [switch_idx(3:2:end) numel(p.tlist)];
            CCW_to_CW_idx_plotR = switch_idx(2:2:end);
        end
    end
    
%% plot figures, if requested
    if options.displayFigs
        CW_color = {'b',[0.85 0.85 1]};
        CCW_color = {'r',[1 0.85 0.85]};
        
        figure

        if p.nLayers == 3
            subplotlocs = [1 3 2]; %on a 2x3 plot
            nRows = 1;
            text_y = [0.95 0.85 0.75 0.65];
        elseif p.nLayers == 5
            subplotlocs = [4 6 2 1 3]; %on a 2x3 plot
            nRows = 2;
            text_y = [-.67 -1 -1.33 -1.67];
        else
            error('I don''t know how to plot this # of layers!');
        end

        plot_x = p.tlist/1000;
        
        for iLayer = 1:p.nLayers
            subplot(nRows,3,subplotlocs(iLayer))
            cla; hold on;
            if iLayer == 3
                % Clockwise patch
                patch(plot_x([CW_to_CCW_idx_plotL ; CW_to_CCW_idx_plotL ; ...
                    CCW_to_CW_idx_plotL ; CCW_to_CW_idx_plotL]),...
                    [zeros(size(CW_to_CCW_idx_plotL)) ;...
                    ones(size(CW_to_CCW_idx_plotL)) ; ...
                    ones(size(CCW_to_CW_idx_plotL)) ; ...
                    zeros(size(CCW_to_CW_idx_plotL))],...
                    CW_color{2},'EdgeColor',CW_color{2})
                % Counter CW patch
                patch(plot_x([CW_to_CCW_idx_plotR ; CW_to_CCW_idx_plotR ; ...
                    CCW_to_CW_idx_plotR ; CCW_to_CW_idx_plotR]),...
                    [zeros(size(CW_to_CCW_idx_plotR)) ;...
                    ones(size(CW_to_CCW_idx_plotR)) ; ...
                    ones(size(CCW_to_CW_idx_plotR)) ; ...
                    zeros(size(CCW_to_CW_idx_plotR))],...
                    CCW_color{2},'EdgeColor',CCW_color{2})
                
                text(1, text_y(1), ['# Reversals: ' num2str(output.nReversals)],...
                    'fontsize',18)
                text(1, text_y(2), ['Freq.: ' num2str(output.Hz) ' Hz'],...
                    'fontsize',18)
                text(1, text_y(3), ['Avg. Dur.: ' num2str(output.avg_dur,3) ' s'],...
                    'fontsize',18)
                text(1, text_y(4), ['CV Dur.: ' num2str(output.CV_dur)],...
                    'fontsize',18)
            end
            p1 = plot(plot_x,p.rA{iLayer},'color',CW_color{1},...
                'Linewidth',2);
            p2 = plot(plot_x,p.rB{iLayer},'color',CCW_color{1},...
                'Linewidth',2);
            if iLayer == subplotlocs(1)
                ylabel('Firing rate')
                xlabel('Time (s)')
                legend([p1 p2], 'Clockwise','Counter CW')
            end
            axis([0 plot_x(end) 0 1])
            title(layernames(iLayer))
            set(gca,'YLim',[0 1],'fontsize',18);
        end
        set(gcf,'color','w','Name',condnames{iCond},'POS',[25 100 1025 420])
    end
end

%% output
output.reversals.CW_to_CCW_idx = CW_to_CCW_idx;
output.reversals.CCW_to_CW_idx = CCW_to_CW_idx;
output.reversals.CW_to_CCW_t = CW_to_CCW_t;
output.reversals.CCW_to_CW_t = CCW_to_CW_t;
output.reversals.switch_t = switch_t;
output.reversals.switch_dur = switch_dur;
output.p = p;
output.wta = wta;
output.options = options;
output.date_run = datestr(now);

end