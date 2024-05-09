function [p] = n_model(p)
%This function is called by n_runModel.m
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

idx = 1; %corresponds to t=0

for t = p.dt:p.dt:p.T
    idx = idx+1;
    for lay = 1:p.nLayers
        
        %defining inputs (stimulus and recurrent connections)
        if lay==1     %left eye      
            inpA = p.iA{lay}(idx) - p.rA{5}(idx-1) - p.rB{5}(idx-1);
            inpB = p.iB{lay}(idx) - p.rA{5}(idx-1) - p.rB{5}(idx-1);
        elseif lay==2 %right eye
            inpA = p.iA{lay}(idx) - p.rA{4}(idx-1) - p.rB{4}(idx-1);
            inpB = p.iB{lay}(idx) - p.rA{4}(idx-1) - p.rB{4}(idx-1);
        elseif lay==3 %summation layer
            inpA = p.rA{1}(idx) + p.rA{2}(idx);
            inpB = p.rB{1}(idx) + p.rB{2}(idx);
        elseif lay==4 %opponency layer (left-right)
            inpA = p.rA{1}(idx) - p.rA{2}(idx);
            inpB = p.rB{1}(idx) - p.rB{2}(idx);
        elseif lay==5 %opponency layer (right-left)
            inpA = p.rA{2}(idx) - p.rA{1}(idx);
            inpB = p.rB{2}(idx) - p.rB{1}(idx);
        end

        %updating drives
        p.dA{lay}(idx) = p.dA{lay}(idx-1) + (p.dt/p.tau)*(-p.dA{lay}(idx-1) + inpA + p.nA{lay}(idx));
        p.dB{lay}(idx) = p.dB{lay}(idx-1) + (p.dt/p.tau)*(-p.dB{lay}(idx-1) + inpB + p.nB{lay}(idx));
         
        %defining normalization pool for each layer
        if lay<=2 %monocular
            pool = [p.sigma, p.dA{1}(idx-1), p.dB{1}(idx-1), p.dA{2}(idx-1), p.dB{2}(idx-1)];
        elseif lay==3 %summation
            pool = [p.sigma, p.dA{3}(idx-1), p.dB{3}(idx-1)];
        elseif lay>=4 %opponency
            pool = [p.sigma_opp, p.dA{4}(idx-1), p.dB{4}(idx-1)]; % mps 20200213 this appears to be a bug
            % should be dA{lay} dB{lay}...
        end
        
        %normalization
        fA = halfExp(p.dA{lay}(idx),2) / (sum(halfExp(pool,2)));
        fB = halfExp(p.dB{lay}(idx),2) / (sum(halfExp(pool,2)));
    
        %update firing rates
        p.rA{lay}(idx) = p.rA{lay}(idx-1) + (p.dt/p.tau)*(-p.rA{lay}(idx-1) + fA);
        p.rB{lay}(idx) = p.rB{lay}(idx-1) + (p.dt/p.tau)*(-p.rB{lay}(idx-1) + fB);

    end
end

