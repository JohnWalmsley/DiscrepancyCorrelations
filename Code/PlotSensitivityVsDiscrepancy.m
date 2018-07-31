function PlotSensitivityVsDiscrepancy()
% This script produces the data used to plot Figure 6
%Simulates response to action potential waveform protocol for literature models and new model (with best fitting parameters to sine wave for cell 5)
% and compares with experimental recording from cell 5

addpath Models
addpath MathworksFileExchange/export_fig/

model = 'hh';
fitting_protocol='sine_wave';
exp_ref = '16713110';


% Identifies temperature for appropriate experiment

if strcmp(exp_ref,'16713110')==1
    temperature = 21.4;
end

% Identifies best fitting parameters to sine wave
[chain,likelihood] = FindingBestFitsAfterMCMC(model,fitting_protocol,exp_ref);

[~,v]= max(likelihood);

MP= chain(v,:);

%Import action potential protocol
cd ../Protocols
V=importdata(['ap_protocol.mat']);
cd ..
cd Code
% Import experimental data
cd ../ExperimentalData
cd(exp_ref)
E=importdata(['ap_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
cd ..
cd ..
cd Code/
% Simulate model with parameters identified as providing bet fit to sine wave
[ I, Y, S, dIdp ] = SimulatingDataSens({'ap'},MP,V,temperature);
figure;plot(I)
hold on 
plot(E)
hold off
figure;plot(I-E, S(:,1));
figure;plot(V)

pacingSpikes = [2501,5700,7200,10561,17693,20721,24245,28097,29392,34398,36888,40508,43710,50420,56095,61254,65954,73246,78247];    
idx = 1 : pacingSpikes(1)-1;
for i = 2 : length(pacingSpikes)-1
    idx = [idx, pacingSpikes(i)+500 : pacingSpikes(i+1)-1];
end
idx = [ idx, pacingSpikes( end ) + 500 : length(V) ];

figure;plot( I(idx))
hold on
plot(E(idx))
hold off

discrepancy = I-E;


F = 96485;
T = 273.15+temperature;
R = 8314;
K_i = 130;
k_o = 4;

erev = ((R*T)/F)*log(k_o/K_i);

O = ones(length(Y(:,1)),1);

Vr = O.*erev;

sensFig = figure;
sensDiscFig = figure;
dIdPfig = figure;
dIdPDiscFig = figure;
for var = 1 : 3
    for par = 1 : 8
       figure( sensDiscFig )
       subplot( 3, 8, par +(var-1)*8 )
       sensitivity = S(:,par +(var-1)*8);
       plot( discrepancy(idx), sensitivity(idx), 'k.' );
       xlim([-0.2,0.2])
       xlabel( 'Discrepancy' )
       ylabel(['$S_{' num2str(var-1) ',' num2str(par-1) '} $'], 'Interpreter', 'Latex')
       figure( sensFig )
       subplot( 3, 8, par +(var-1)*8 )
       sensitivity = S(:,par +(var-1)*8);
       plot( sensitivity, 'k.' );
       xlim([0 76050]);
       xlabel(['$S_{' num2str(var-1) ',' num2str(par-1) '} $'], 'Interpreter', 'Latex')
       
   if var == 3
       figure( dIdPfig )
       subplot( 2, 8, par )
       plot(sensitivity )
       subplot( 2, 8, par+8 )
       plot( MP(end)*sensitivity.*(V-Vr))
       figure(dIdPDiscFig);
       subplot( 2, 4, par )
       plot(discrepancy, MP(end)*sensitivity.*(V-Vr) );
   end
    end
end
% Prepare and save data
new_model_data=[[0:0.0001:(length(V)/10000)-0.0001]',I];
experimental_data=[[0:0.0001:(length(V)/10000)-0.0001]',E];
protocol_data=[[0:0.0001:(length(V)/10000)-0.0001]',V];
save ap_new_model_prediction.txt new_model_data -ascii
save ap_experimental_data.txt experimental_data -ascii
save ap_protocol_data.txt protocol_data -ascii
