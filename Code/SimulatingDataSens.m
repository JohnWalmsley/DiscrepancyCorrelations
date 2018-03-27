function [I, Y, S, dIdp, flux, time, G ] = SimulatingDataSens(protocol,params,V,temperature)
%--------------------------------------------------------------------------------------------------------------------------------------
% In this function the initial conditions are identified for the model type and appropriate MexFunction for model is called (with model parameters
% passed to MexFunction) to simulate current in response to specified protocol.

%Temperature Converted to Kelvins

T = 273.15+temperature;

%Initial conditions for HH model

IC = [0.0,0.0,0.0];

%-----------------------------------------------------------------------------------------------------------------------------------------%
% Defines protocol number and protocol length (ms) for each protocol
%-----------------------------------------------------------------------------------------------------------------------------------------%
cd ..
% For non-sine wave protocol the protocol is also passed as additional parameters, avoiding detailed protocol description in Mex function
ProtocolLength = (length(V)/10);
if strcmp(protocol,'sine_wave')==1
    vv = -80; % First voltage for sensitivity initial conditions
    protocol_number=1;
end

if strcmp(protocol,'ap')==1
    protocol_number=7;
    cd Protocols
    vv=importdata('ap_protocol.mat');
    cd ..
    params(end)
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end

if strcmp(protocol,'ap_lag')==1
    protocol_number=7;
    cd Protocols
    vv=importdata('ap_lag_protocol.mat');
    cd ..
    params(end)
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end

if strcmp(protocol,'activation_kinetics_1')==1
    protocol_number=9;
    cd Protocols
    vv=importdata('activation_kinetics_1_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end

if strcmp(protocol,'activation_kinetics_2')==1
    protocol_number=11;
    cd Protocols
    vv=importdata('activation_kinetics_2_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end
if strcmp(protocol,'inactivation')==1
    protocol_number=12;
    cd Protocols
    vv=importdata('inactivation_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end
if strcmp(protocol,'deactivation')==1
    protocol_number=13;
    cd Protocols
    vv=importdata('deactivation_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end
if strcmp(protocol,'steady_activation')==1
    protocol_number=14;
    cd Protocols
    vv=importdata('steady_activation_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end

if strcmp(protocol,'grandi_2hz')==1
    
    protocol_number=18;
    cd Protocols
    vv=importdata('grandi_2hz_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end
if strcmp(protocol,'grandi_fast')==1
    
    protocol_number=19;
    cd Protocols
    vv=importdata('grandi_fast_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end
if strcmp(protocol,'repeated_activation_step')==1
    protocol_number=20;
    cd Protocols
    vv=importdata('repeated_activation_step_protocol.mat');
    cd ..
    temp_con = params(end);
    params = [params(1:end-1),vv',temp_con];
end
cd Code
%-----------------------------------------------------------------------------------------------------------------------------------------%
% Protocol number is added to start of parameter vector to be passed in mex function
%-----------------------------------------------------------------------------------------------------------------------------------------%
params = [protocol_number,params];

%-----------------------------------------------------------------------------------------------------------------------------------------%
% Using the approriate mex function as determined by the model type, the time vector, initial conditions and parameter vector (including the
% initial element being the protocol_number corresponding to the protocol being simulated) is passed to the mex function. The solution X
% is the probability of being in the open/activated state where there is just one vector, or the solutions of the probability of being in
% multiple different activated or the inactivated state where there are multiple solution outputs.
%-----------------------------------------------------------------------------------------------------------------------------------------%

% sensitivity initial condition is based on the steady state for the voltage at the first time
% point i nthe protocol

A_V = CalculateJacobianMatrix( vv(1), params );
k_V = CalculateSensitivityRhsVectors( vv(1), IC', params );
s_star = -inv(A_V)*k_V;
s_star = [ s_star(1,:) s_star(2,:) s_star(3,:) ]; % Unpack into row vector
time = 0:0.1:ProtocolLength;
[ Y, S, flux ] = MexHHSens( time, IC, params, s_star );

%Calculates reversal potential using Nernst equation to be used in current calculation equation. Temperature, extracellular and intracellular potassium concentration correspond to those used in experiment.

F = 96485;

R = 8314;
K_i = 130;
k_o = 4;

erev = ((R*T)/F)*log(k_o/K_i);

O = ones(length(Y(:,1)),1);

Vr = O.*erev;

%-----------------------------------------------------------------------------------------------------------------------------------------%

% Current calculations.
% For Markov models, current is calculated according to
%                         I = g*O*(V-Vr)
% where the conductance g is the final parameter in the parameter set for each model and O is the probability of being in the open state,
% which corresponds to the solution X for these models.

% Similarly, the current formulation for Hodgkin-Huxley style models is typically of the form
%                         I =  g*Xr*Rr*(V-Vr)
% so for model for which the output of the Mex function X represents Xr*Rr,
% the current can be expressed in the same way as for the Markov model formulations.
G = params( length( params ) );
I = G.*Y(:,3).*(V-Vr);
dIdp = params(length(params))*bsxfun( @times,S(:,17:24),(V-Vr));
time = time(2:end)';
