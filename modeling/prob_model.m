function p = prob_model(options)
% usage: p = prob_model(options)
%
% options = structure with lots of fields, all of which have default values
%           if not provided.
% p = model parameters at the END of running a big for loop to calculate
%     differential equation results ("numerically approximated with Euler�s
%     Method")
%
% A = stimulus type #1, i.e., left-ward drifting
% B = stimulus type #2, i.e., right-ward drifting
%
% Layer 1 = selective for dimension #1, i.e., near-selective
% Layer 2 = selective for dimension #2, i.e., far-selective
% Layer 3 = summation across dimensions
% Layer 4 = opponency for dim. 1 - 2 (near - far)
% Layer 5 = opponency for dim. 2 - 1 (far - near)
%
% Adapted from:
% Said and Heeger (2013) A model of binocular rivalry and cross-orientation
% suppression. PLOS Computational Biology.
%
% mps 20191024
%% set parameters
if ~exist('options','var') % parameters
    options = [];
end
if ~isfield(options,'dt')
    options.dt = 50; % ms
end
if ~isfield(options,'sim_duration')
    options.sim_duration = 20000; % ms
end
options.n_time_points = options.sim_duration/options.dt+1;
if ~isfield(options,'tlist')
    options.tlist = 0:options.dt:options.sim_duration; % percent contrast
end
if ~isfield(options,'nLayers')
    options.nLayers = 3; % let's start with modeling 3 layers by default
    % where 1st = near depth, 2 = far depth, 3 = depth invariant
end
if ~isfield(options,'tau')
    options.tau = 3250; % ms
end
if ~isfield(options,'sigma')
    options.sigma = 0.5; % percent contrast
end
if ~isfield(options,'sigma_opp')
    options.sigma_opp = 0.9; % percent contrast
end
if ~isfield(options,'normGain')
    options.normGain = 1; % scale normalization factor
end
if ~isfield(options,'noiseAmp')
    options.noiseAmp = 0.2; % arb. units
end
if ~isfield(options,'attnGain')
    options.attnGain = 1; % scale attention factor
end

%% now pre-allocate model layer pieces
if ~isfield(options,'iA') % input A
    for iLayer = 1:nLayers
        options.iA{iLayer} = ones(1,options.n_time_points);
    end
end
if ~isfield(options,'iB') % input B
    for iLayer = 1:nLayers
        options.iB{iLayer} = ones(1,options.n_time_points);
    end
end
if ~isfield(options,'dA') % drive A
    for iLayer = 1:nLayers
        options.dA{iLayer} = zeros(1,options.n_time_points);
    end
end
if ~isfield(options,'dB') % drive B
    for iLayer = 1:nLayers
        options.dB{iLayer} = zeros(1,options.n_time_points);
    end
end
if ~isfield(options,'rA') % resp A
    for iLayer = 1:nLayers
        options.rA{iLayer} = zeros(1,options.n_time_points);
    end
end
if ~isfield(options,'rB') % resp B
    for iLayer = 1:nLayers
        options.rB{iLayer} = zeros(1,options.n_time_points);
    end
end
if ~isfield(options,'nA') % noise A
    for iLayer = 1:nLayers
        options.nA{iLayer} = makeNoise(options);
    end
end
if ~isfield(options,'nB') % noise B
    for iLayer = 1:nLayers
        options.nB{iLayer} = makeNoise(options);
    end
end

p = options; % change variable here to save input options later,
% without tranforming...
%% run model
start_idx = 1; % t = 0 ms, initial conditions only

for t_idx = (start_idx+1):p.n_time_points
    for iLayer = 1:p.nLayers
        % defining inputs (stimulus and recurrent connections)
        if iLayer == 1     % near selective      
            inpA = p.iA{iLayer}(t_idx) - p.rA{5}(t_idx-1) - p.rB{5}(t_idx-1);
            inpB = p.iB{iLayer}(t_idx) - p.rA{5}(t_idx-1) - p.rB{5}(t_idx-1);
        elseif iLayer == 2 % far selective
            inpA = p.iA{iLayer}(t_idx) - p.rA{4}(t_idx-1) - p.rB{4}(t_idx-1);
            inpB = p.iB{iLayer}(t_idx) - p.rA{4}(t_idx-1) - p.rB{4}(t_idx-1);
        elseif iLayer == 3 % summation layer
            inpA = p.rA{1}(t_idx) + p.rA{2}(t_idx);
            inpB = p.rB{1}(t_idx) + p.rB{2}(t_idx);
        elseif iLayer == 4 % opponency layer (near-far)
            inpA = p.rA{1}(t_idx) - p.rA{2}(t_idx);
            inpB = p.rB{1}(t_idx) - p.rB{2}(t_idx);
        elseif iLayer == 5 % opponency layer (far-near)
            inpA = p.rA{2}(t_idx) - p.rA{1}(t_idx);
            inpB = p.rB{2}(t_idx) - p.rB{1}(t_idx);
        end

        % define attention
        if iLayer == 3
            if p.rA{iLayer}(t_idx-1) > p.rB{iLayer}(t_idx-1) % A is higher, attend A
                scaleAttnA = p.attnGain;
                scaleAttnB = 1;
            elseif p.rA{iLayer}(t_idx-1) < p.rB{iLayer}(t_idx-1) % A is higher, attend A
                scaleAttnA = 1;
                scaleAttnB = p.attnGain;
            else % they are... equal?
                scaleAttnA = 1;
                scaleAttnB = 1;
            end
        else % only in layer 3
            scaleAttnA = 1;
            scaleAttnB = 1;
        end
        
        % updating drives
        p.dA{iLayer}(t_idx) = p.dA{iLayer}(t_idx-1) + (p.dt/p.tau)*...
            (-p.dA{iLayer}(t_idx-1) + inpA .* scaleAttnA + p.nA{iLayer}(t_idx));
        p.dB{iLayer}(t_idx) = p.dB{iLayer}(t_idx-1) + (p.dt/p.tau)*...
            (-p.dB{iLayer}(t_idx-1) + inpB .* scaleAttnB + p.nB{iLayer}(t_idx));

        % defining normalization pool for each layer
        if iLayer <= 2 % monocular
            pool = [p.sigma, p.dA{1}(t_idx-1), p.dB{1}(t_idx-1), ...
                p.dA{2}(t_idx-1), p.dB{2}(t_idx-1)];
        elseif iLayer == 3 % summation
            pool = [p.sigma, p.dA{3}(t_idx-1), p.dB{3}(t_idx-1)];
        elseif iLayer >= 4 % opponency
            pool = [p.sigma_opp, p.dA{iLayer}(t_idx-1), p.dB{iLayer}(t_idx-1)];
            % n.b., this differs from the code shipped by Said and Heeger,
            % because I think their code has a bug (was set to only layer
            % 4) -- mps 20200213
        end
        
        % normalization
        fA = halfExp(p.dA{iLayer}(t_idx),2) / ...
            (p.normGain*(sum(halfExp(pool,2))));
        fB = halfExp(p.dB{iLayer}(t_idx),2) / ...
            (p.normGain*(sum(halfExp(pool,2))));
        
        % update firing rates
        p.rA{iLayer}(t_idx) = p.rA{iLayer}(t_idx-1) + (p.dt/p.tau)*...
            (-p.rA{iLayer}(t_idx-1) + fA);
        p.rB{iLayer}(t_idx) = p.rB{iLayer}(t_idx-1) + (p.dt/p.tau)*...
            (-p.rB{iLayer}(t_idx-1) + fB);

    end
end

%% output
p.options = options;
p.date_run = datestr(now);

end