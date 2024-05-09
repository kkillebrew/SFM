function determ_model
% usage: just run determ_model.m
%
% mps 20191024... doesn't work right
%% set up
% from Wilson 2003
nLayers = 3; % left eye, right eye, binocular

tau = 20;
tau_I = 11;
tau_H = 900;
g = 45;
h = 0.47;
H_Volt = 10;
scale_E = 100;

ex_gain_layer2 = .75;
scale_g_layer2 = 1.53;
scale_recurrent_fb = 0.002;

% initial conditions
E_v = [1 0 0]; % stimuli vertical left, horiz right, only in 1st layer
E_h = [0 1 0]; % if you bias one of the inputs slightly, e.g., 1.01, then
               % you get a dominant percept

I_v = [0 0 0];
I_h = [0 0 0];

H_v = [0 0 0];
H_h = [0 0 0];

input_Volt_v = [10 0];
input_Volt_h = [0 10];

start_idx = 1; % t = 0 ms

dt = 5; % time-step (ms)

sim_time = 10000; % ms

n_time_points = sim_time/dt+1;

output.input_v = [];

for t_idx = (start_idx+1):n_time_points; 
    for iLayer = 1:nLayers
%% Layer input
if iLayer == 1 % left eye
    % input is vertical in left eye
    input_v = input_Volt_v(iLayer) + scale_recurrent_fb .* E_v(t_idx-1, 3);
    input_h = input_Volt_h(iLayer) + scale_recurrent_fb .* E_h(t_idx-1, 3);
    
    % inhibition from the opposite layer & channel
    inhib_v = I_h(t_idx-1, 2);
    inhib_h = I_v(t_idx-1, 2);
    
    scale_g = 1.*g;
    
elseif iLayer == 2 % right eye
    % input is horizontal in right eye
    input_v = input_Volt_v(iLayer) + scale_recurrent_fb .* E_v(t_idx-1, 3);
    output.input_v = [output.input_v input_v];
    input_h = input_Volt_h(iLayer) + scale_recurrent_fb .* E_h(t_idx-1, 3);
    
    % inhibition from the opposite layer & channel
    inhib_v = I_h(t_idx-1, 1);
    inhib_h = I_v(t_idx-1, 1);
    
    scale_g = 1.*g;
    
elseif iLayer == 3 % binocular summation
    input_v = ex_gain_layer2 .* ( E_v(t_idx-1, 1) + E_v(t_idx-1, 2) );
    input_h = ex_gain_layer2 .* ( E_h(t_idx-1, 1) + E_h(t_idx-1, 2) );
    
    % inhibition from the same layer, opposite channel
    inhib_v = I_h(t_idx-1, 3);
    inhib_h = I_v(t_idx-1, 3);
    
    scale_g = scale_g_layer2.*g;
else
    error('What layer is this??');
end
%% E cells

E_v(t_idx, iLayer) = E_v(t_idx-1, iLayer) + ...
    ( dt / tau ) .* ( -E_v(t_idx-1, iLayer) + ...
    scale_E .* ( max([0, input_v - scale_g.*inhib_v]) ).^2 ...
    ./ ( (H_Volt + H_v(t_idx-1, iLayer)).^2 + ...
    ( max([0, input_v - scale_g.*inhib_v]) ).^2 ) );

E_h(t_idx, iLayer) = E_h(t_idx-1, iLayer) + ...
    ( dt / tau ) .* ( -E_h(t_idx-1, iLayer) + ...
    scale_E .* ( max([0, input_h - scale_g.*inhib_h]) ).^2 ...
    ./ ( (H_Volt + H_h(t_idx-1, iLayer)).^2 + ...
    ( max([0, input_h - scale_g.*inhib_h]) ).^2 ) );

%% I cells
I_v(t_idx, iLayer) = I_v(t_idx-1, iLayer) + ( dt / tau_I ) .* ...
    ( -I_v(t_idx-1, iLayer) + E_v(t_idx-1, iLayer) );

I_h(t_idx, iLayer) = I_h(t_idx-1, iLayer) + ( dt / tau_I ) .* ...
    ( -I_h(t_idx-1, iLayer) + E_h(t_idx-1, iLayer) );

%% H current (self hyper-polarizing)
H_v(t_idx, iLayer) = H_v(t_idx-1, iLayer) + ( dt / tau_H ) .* ...
    ( -H_v(t_idx-1, iLayer) + h.*E_v(t_idx-1, iLayer) );

H_h(t_idx, iLayer) = H_h(t_idx-1, iLayer) + ( dt / tau_H ) .* ...
    ( -H_h(t_idx-1, iLayer) + h.*E_h(t_idx-1, iLayer) );

    end
end

figure
subplot(3,1,1)
hold on
plot(E_v(:,3),'b')
title('binocular')
plot(E_h(:,3),'r')
%axis([0 10000 0 100])
legend('vert','horiz')

subplot(3,1,2)
hold on
plot(E_v(:,2),'b')
title('right')
plot(E_h(:,2),'r')
%axis([0 10000 0 50])

subplot(3,1,3)
hold on
plot(E_v(:,1),'b')
title('left')
plot(E_h(:,1),'r')
%axis([0 10000 0 50])

end