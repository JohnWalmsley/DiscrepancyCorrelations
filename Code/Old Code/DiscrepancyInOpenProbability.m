%% Defining the discrepancy in the open probability:
% G * (O+d)*(V-V_E) = data
% Therefore data - Isim = G*d*(V-V_E)
% and so d = ( data-Isim ) / (G*(V-V_E)) outside V_E + eps

addpath( genpath( '../../SharedFunctions' ) )
addpath ../../../MathworksFileExchange/getLinearDependent/

%% Configuration

training_protocol = 'sine_wave';
prediction_protocol = 'ap';
exp_ref = '16713110';

method = 'lasso';
%conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'sensitivity', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise', 'thin' };
%conf = { 'const', 'variables', 'I', 'rates', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise', 'thin' };
%conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise' };
conf = { 'const', 'variables' };

% eps determines voltage proximity to reversal when calculating discrepancy
% in open probability
eps = 5;

%% Training Data

[ variables, I, rates, fluxes, sensitivity, dIdp, voltage, dVdt, time, G ]  ...
                            = CalculateVariables( training_protocol, exp_ref );
figure;plot(voltage)
[ discrepancy, exp_data ] = CalculateDiscrepancy( exp_ref, training_protocol, I );
voltage=importdata(['../Protocols/' training_protocol '_protocol.mat']);

% Set up constants 
F = 96485;
R = 8314;
K_i = 130;
k_o = 4;
T = GetTemperature( exp_ref ) + 273.5;
erev = ((R*T)/F)*log(k_o/K_i);

% Remove the pacing spikes from the simulated and experimental data
idx1 = GetNoSpikeIdx( training_protocol, length( time ) );

voltage = voltage( idx1 );
time = time( idx1 );
variables = variables( idx1, : );
I = I(idx1);
rates = rates( idx1, :) ;
fluxes = fluxes( idx1, : );
exp_data = exp_data( idx1 );
sensitivity = sensitivity( idx1, :);
dIdp = dIdp(idx1,:);
dVdt = dVdt(idx1, :);
discrepancy = discrepancy( idx1 );

% Remove the indices near the reversal potential for Odisc
idx = find( abs( voltage - erev ) > eps );
idx_nearvoltage = find( abs( voltage - erev ) < eps );
discrepancy_o = discrepancy( idx ) ./ ( G * ( voltage( idx ) - erev ) );

% Construct the discrepancy model

% calculate rates fluxes sensitivity
[data_matrix_sw, factors_sw ] = BuildDataMatrix( conf, variables(idx,:), I(idx), rates(idx,:), fluxes(idx,:), ...
                                                        sensitivity( idx,:), dIdp(idx,:), voltage(idx), dVdt(idx,:), time(idx) );
[ ~, data_matrix, model ] = ...
                        LinearModelOfDiscrepancyWithInput( method, conf, discrepancy_o, data_matrix_sw );
if ~isa( model, 'LinearModel' )
    included_idx = find( abs( model ) > 0 );
    names = GetNames( conf );
    included_names = names( included_idx );
end
% predicted discrepancy in O
if isa( model, 'LinearModel' )
    if any( strcmp( conf, 'const' ) )
        oDisc_predicted = predict( model, data_matrix_sw( :, 2 : end ) );
    else
       oDisc_predicted = predict( model, data_matrix_sw ); 
    end
    model.Coefficients
else
    oDisc_predicted = data_matrix_sw*model;
end
                 
fig=figure;plot( time(idx), I(idx), 'k-', time(idx), I(idx) +  G*oDisc_predicted.*(voltage(idx)-erev), 'r-', time(idx), exp_data(idx), 'b-' )
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ training_protocol '- Training' ], 'Interpreter', 'none' )
legend( {'I_sim', 'I_sim +G.D.(V-V_E)', 'I_exp'} )
saveas( fig, [ 'Figures/ODiscModel_Train_' training_protocol '.tif' ] )
saveas( fig, [ 'Figures/ODiscModel_Train_' training_protocol '.fig' ] )

disp( [ 'Original error in ' training_protocol ' (training): ' num2str( sqrt( sum( discrepancy_o.^2) ) ) ] )
disp( [ 'New error in ' training_protocol ' (training): ' num2str( sqrt( sum( ( discrepancy_o - G*oDisc_predicted.*(voltage(idx)-erev) ).^2 ) ) ) ] )
    
%% Prediction Data

[ variables, I, rates, fluxes, sensitivity, dIdp, voltage, dVdt, time, G ]  ...
                                    = CalculateVariables( prediction_protocol, exp_ref );
[ discrepancy, exp_data ] = CalculateDiscrepancy( exp_ref, prediction_protocol, I );

voltage=importdata(['../Protocols/' prediction_protocol '_protocol.mat']);

% Remove the pacing spikes from the simulated and experimental data
idx1 = GetNoSpikeIdx( prediction_protocol, length( time ) );

voltage = voltage( idx1 );
time = time( idx1 );
variables = variables( idx1, : );
I = I(idx1);
rates = rates( idx1, :) ;
fluxes = fluxes( idx1, : );
exp_data = exp_data( idx1 );
sensitivity = sensitivity( idx1, :);
dIdp = dIdp(idx1,:);
dVdt = dVdt(idx1, :);
discrepancy = discrepancy( idx1 );

% Remove the indices near the reversal potential for Odisc
idx = find( abs( voltage - erev ) > eps );
idx_nearvoltage = find( abs( voltage - erev ) < eps );
discrepancy_o = discrepancy( idx ) ./ ( G * ( voltage( idx ) - erev ) );

% build the data matrix used for prediction
[data_matrix_ap, ~ ] = BuildDataMatrix( conf, variables(idx,:), I(idx), rates(idx,:), fluxes(idx,:), ...
                                                        sensitivity( idx,:), dIdp(idx,:), voltage(idx), dVdt(idx,:), time(idx), factors_sw );% predicted discrepancy in O
if isa( model, 'LinearModel' )
    if any( strcmp( conf, 'const' ) )
        oDisc_predicted = predict( model, data_matrix_ap( :, 2 : end ) );
    else
        oDisc_predicted = predict( model, data_matrix_ap ); 
    end
else
    oDisc_predicted = data_matrix_ap*model;
end
                                                    
fig=figure;plot( time(idx), I(idx), 'k-', time(idx), I(idx) + G*oDisc_predicted.*(voltage(idx)-erev), 'r-', time(idx), exp_data(idx), 'b-' )
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ prediction_protocol ' - prediction' ], 'Interpreter', 'none' )
legend( {'I_sim', 'I_sim +G.D.(V-V_E)', 'I_exp'} )

saveas( fig, [ 'Figures/ODiscModel_Predict_' prediction_protocol '.tif' ] )
saveas( fig, [ 'Figures/ODiscModel_Predict_' prediction_protocol '.fig' ] )

% Calculate MSE errors
disp( [ 'Original error in ' prediction_protocol ' (prediction): ' num2str( sqrt( sum( discrepancy_o.^2) ) ) ] )
disp( [ 'New error in ' prediction_protocol ' (prediction): ' num2str( sqrt( sum( ( discrepancy_o - G*oDisc_predicted.*(voltage(idx)-erev) ).^2 ) ) ) ] )
                                                    