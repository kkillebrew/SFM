function [noise] = makeNoise(p)
% usage: [noise] = makeNoise(p)
%
% Adapted from:
% Said and Heeger (2013) A model of binocular rivalry and cross-orientation
% suppression. PLOS Computational Biology.
%
% mps 20191025
%% set parameters
if ~exist('p','var') % parameters
    p = [];
end
if ~isfield(p,'noiseAmp')
    p.noiseAmp = 0.2; % arb. units
end
if ~isfield(p,'noisefilter_t')
    p.noisefilter_t = 800; % (ms)
end
if ~isfield(p,'tlist')
    error('no p.tlist provided... I quit.');
end
%% run
noise = p.noiseAmp*randn(1,length(p.tlist));
tkern = normpdf(p.tlist,0,p.noisefilter_t);
tkern = tkern/norm(tkern);
noise = cconv2(noise,tkern);

end
