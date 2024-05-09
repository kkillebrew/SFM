function [] = n_runModel()
%Run this function.
%This function will loop through 3 conditions: Dichoptic gratings,
%monocular plaids, and binocular plaids
%The word 'layer' is used in the following sense:
%Layer 1 = Left monocular neurons
%Layer 2 = Right monocular neurons
%Layer 3 = Summation neurons
%Layer 4 = Left-Right opponency neurons
%Layer 5 = Right-Left opponency neurons
%A and B refer to the two stimulus orientations.
%
%If you use this code, please cite
%Said and Heeger (2013) A model of binocular rivalry and cross-orientation
%suppression. PLOS Computational Biology.

c = .5; %contrast
iA_amp_opts = {[0 c]; [c 0]; [c c]}; %dichoptic grating, monocular plaid, binocular plaids
iB_amp_opts = {[c 0]; [c 0]; [c c]}; %dichoptic grating, monocular plaid, binocular plaids
condnames =  {'Dich. Gratings', 'Mon. Plaid', 'Bin. Plaid'};
layernames =  {'L. Monocular', 'R. Monocular', 'Summation', 'L-R Opponency', 'R-L Opponency'};
p.sigma         = .5;       %semisaturation constant
p.sigma_opp     = .9;       %semisaturation constant for opponency cells
p.tau           = 50;       %time constant (ms)
p.dt            = 5;        %time-step (ms)
p.T             = 20000;    %duration (ms)
p.noisefilter_t = 800;      %(ms)
p.noiseamp      = .05;      
p.nLayers       = 5;        %set to 3 for conventional model, 5 for opponency model
p.nt            = p.T/p.dt+1;
p.tlist         = 0:p.dt:p.T;

%Initializing time-courses for neuron (d)rives, (r)esponses, and (n)oise.
%Each neuron is tuned to either orientation A or B.
for lay=1:5 %go through maximum possible layers. This way, if there are <5 layers, the feedback can be zero.
    p.dA{lay}   = zeros(1,p.nt);
    p.dB{lay}   = zeros(1,p.nt);
    p.rA{lay}   = zeros(1,p.nt);
    p.rB{lay}   = zeros(1,p.nt);
    p.nA{lay}   = n_makeNoise(p);
    p.nB{lay}   = n_makeNoise(p);
end


subplotlocs = [4 6 2 1 3]; %on a 2x3 plot

%loop through conditions
for cond = 1:3;
    %stimulus inputs to monocular layers
    for lay = 1:2
        p.iA{lay} = iA_amp_opts{cond}(lay)*ones(1,p.nt);
        p.iB{lay} = iB_amp_opts{cond}(lay)*ones(1,p.nt);
    end
     
    %run the model
    p = n_model(p);
    
    %compute WTA index from summation layer
    wta(cond) = nanmean(abs(p.rA{3}-p.rB{3})./(p.rA{3}+p.rB{3}));

    cpsFigure(3,1)
    set(gcf,'Name',condnames{cond})
    for lay = 1:p.nLayers
        subplot(2,3,subplotlocs(lay))
        cla; hold on;
        p1 = plot(p.tlist/1000,p.rA{lay},'color',[1 0 1]);
        p2 = plot(p.tlist/1000,p.rB{lay},'color',[0 0 1]);
        legend([p1 p2], 'A','B')
        ylabel('Firing rate')
        xlabel('Time (s)')
        title(layernames(lay))
        set(gca,'YLim',[0 1]);
    end
end


figure
cla; hold on;
barh(3,wta(1), 'FaceColor', [.6 .6 .6]);
barh(2,wta(2), 'FaceColor', [.6 .6 .6]);
barh(1,wta(3), 'FaceColor', [.6 .6 .6]);
xlabel('Winner-take-all index','FontSize',20)
set(gca,'YTick', [1 2 3], 'YLim', [.3 3.7], 'XTick', [0 .2 .4 .6 .8 1], 'XLim', [0 1], 'FontSize', 14)
set(gca,'YTickLabel', {'Bin. plaids', 'Mon. plaid', 'Dich. gratings'});
set(gca,'FontSize',20);

