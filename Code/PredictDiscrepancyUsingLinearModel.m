function PredictDiscrepancyUsingLinearModel( exp_ref, training_protocol, prediction_protocol, method, only_core, conf )


%% Defining the discrepancy in the open probability:
% G * (O+d)*(V-V_E) = data
% Therefore data - Isim = G*d*(V-V_E)
% and so d = ( data-Isim ) / (G*(V-V_E)) outside V_E + eps

addpath( genpath( '../../SharedFunctions' ) )
addpath ../../../MathworksFileExchange/getLinearDependent/

%% Configuration

if nargin < 5
    only_core = 0;
end

if nargin < 6
    conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'sensitivity', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise', 'thin' };
    %conf = { 'const', 'variables', 'I', 'rates', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise', 'thin' };
    %conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise' };
    %conf = { 'const', 'variables' };
end

% eps determines voltage proximity to reversal when calculating discrepancy
% in open probability
eps = 5;

%% Training Data
% Calculate variables and discrepancies

[ variables_tr, I_tr, rates_tr, fluxes_tr, sensitivity_tr, dIdp_tr, ~, dVdt_tr, time_tr, G ]  ...
                            = CalculateVariables( training_protocol, exp_ref );
[ discrepancy_tr, exp_data_tr ] = CalculateDiscrepancy( exp_ref, training_protocol, I_tr );
voltage_tr=importdata(['../Protocols/' training_protocol '_protocol.mat']);

% Set up constants 
F = 96485;
R = 8314;
K_i = 130;
k_o = 4;
T = GetTemperature( exp_ref ) + 273.5;
erev = ((R*T)/F)*log(k_o/K_i);

%--------------------------------------------------------------------------
% Discrepancy
%--------------------------------------------------------------------------

% Remove the pacing spikes from the simulated and experimental data

idx_nospike = GetNoSpikeIdx( training_protocol, length( time_tr ) );
if only_core == 1
    idx_core = GetCoreOfProtocolIdx( training_protocol );
else
    idx_core = 1 : length( voltage_tr );
end
idx_d = intersect( idx_nospike, idx_core );

% variables for calculating discrepancy model
voltage_tr_d = voltage_tr( idx_d );
time_tr_d = time_tr( idx_d );
variables_tr_d = variables_tr( idx_d, : );
I_tr_d = I_tr(idx_d);
rates_tr_d = rates_tr( idx_d, :) ;
fluxes_tr_d = fluxes_tr( idx_d, : );
exp_data_tr_d = exp_data_tr( idx_d );
sensitivity_tr_d = sensitivity_tr( idx_d, :);
dIdp_tr_d = dIdp_tr(idx_d,:);
dVdt_tr_d = dVdt_tr(idx_d, :);
discrepancy_tr_d = discrepancy_tr( idx_d );

% calculate rates fluxes sensitivity
[data_matrix_tr_d, factors_tr_d, names_tr_d ] = BuildDataMatrix( conf, variables_tr_d, I_tr_d, rates_tr_d, fluxes_tr_d, ...
                                                       sensitivity_tr_d, dIdp_tr_d, voltage_tr_d, dVdt_tr_d, time_tr_d );                                                

% Are any members of the data matrix linearly dependent on one another?
[ ~, ind_idx, grps ] = getLinearIndependent(data_matrix_tr_d( :, 2 : end ), 1); % leave out the constant term.
dependent_grps = ( cellfun( 'length', grps ) > 1 );
if any ( dependent_grps )
    disp( 'Warning: Some predictors are linearly dependent between classes for additive discrepancy!' )
    nonind_grp_idxs = find( dependent_grps );
    for i=nonind_grp_idxs
        dep_idxs = grps{ i };
        namestring = sprintf( '%s ', names_tr_d{ dep_idxs+1 } ); % +1 due to removing const from mx.
        disp( [ namestring 'are linearly dependent' ] );
    end
end                                                  
                                                   
[ ~, data_matrix_tr_d_used, model_d, lassofig_d ] = ...
                        LinearModelOfDiscrepancyWithInput( method, conf, discrepancy_tr_d, data_matrix_tr_d, names_tr_d );

if ~isa( model_d, 'LinearModel' )
    included_idx = find( abs( model_d ) > 0 );
    included_names = names_tr_d( included_idx );
end
% predicted discrepancy
if isa( model_d, 'LinearModel' )
    if any( strcmp( conf, 'const' ) )
        disc_modelled = predict( model_d, data_matrix_tr_d( :, 2 : end ) );
    else
        disc_modelled = predict( model_d, data_matrix_tr_d ); 
    end
    %model.Coefficients
else
    disc_modelled = data_matrix_tr_d*model_d;
end

%--------------------------------------------------------------------------
% Discrepancy in Open
%--------------------------------------------------------------------------
                               
% Remove the indices near the reversal potential for Odisc
idx_erev = find( abs( voltage_tr - erev ) > eps );
idx_od = intersect( idx_d, idx_erev );

voltage_tr_od = voltage_tr( idx_od );
time_tr_od = time_tr( idx_od );
variables_tr_od = variables_tr( idx_od, : );
I_tr_od = I_tr(idx_od);
rates_tr_od = rates_tr( idx_od, :) ;
fluxes_tr_od = fluxes_tr( idx_od, : );
exp_data_tr_od = exp_data_tr( idx_od );
sensitivity_tr_od = sensitivity_tr( idx_od, :);
dIdp_tr_od = dIdp_tr(idx_od,:);
dVdt_tr_od = dVdt_tr(idx_od, :);

discrepancy_tr_od = discrepancy_tr( idx_od ) ./ ( G * ( voltage_tr( idx_od ) - erev ) );

% Construct the discrepancy model

% calculate rates fluxes sensitivity
[data_matrix_tr_od, factors_tr_od, names_tr_od ] = BuildDataMatrix( conf, variables_tr_od, I_tr_od, rates_tr_od, fluxes_tr_od, ...
            sensitivity_tr_od, dIdp_tr_od, voltage_tr_od, dVdt_tr_od, time_tr_od );
% Are any members of the data matrix linearly dependent on one another?
[ ~, ind_idx, grps ] = getLinearIndependent(data_matrix_tr_od( :, 2 : end ), 1); % leave out the constant term.
dependent_grps = ( cellfun( 'length', grps ) > 1 );
if any ( dependent_grps )
    disp( 'Warning: Some predictors are linearly dependent between classes for additive discrepancy!' )
    nonind_grp_idxs = find( dependent_grps );
    for i=nonind_grp_idxs
        dep_idxs = grps{ i };
        namestring = sprintf( '%s ', names_tr_od{ dep_idxs+1 } ); % +1 due to removing const from mx.
        disp( [ namestring 'are linearly dependent' ] );
    end
end    
                                                    
[ ~, data_matrix_tr_od_used, model_od, lassofig_od ] = ...
        LinearModelOfDiscrepancyWithInput( method, conf, discrepancy_tr_od, data_matrix_tr_od, names_tr_od );

                    
if ~isa( model_od, 'LinearModel' )
    included_idx = find( abs( model_od ) > 0 );
    included_names = names_tr_od( included_idx );
end
% predicted discrepancy in O
if isa( model_od, 'LinearModel' )
    if any( strcmp( conf, 'const' ) )
       odisc_modelled = predict( model_od, data_matrix_tr_od( :, 2 : end ) );
    else
       odisc_modelled = predict( model_od, data_matrix_tr_od ); 
    end
    %model.Coefficients
else
    odisc_modelled = data_matrix_tr_od*model_od;
end
                 
%% Prediction Data

[ variables_pr, I_pr, rates_pr, fluxes_pr, sensitivity_pr, dIdp_pr, voltage_pr, dVdt_pr, time_pr, ~ ]  ...
                                    = CalculateVariables( prediction_protocol, exp_ref );
[ discrepancy_pr, exp_data_pr ] = CalculateDiscrepancy( exp_ref, prediction_protocol, I_pr );
voltage_pr=importdata(['../Protocols/' prediction_protocol '_protocol.mat']);

%--------------------------------------------------------------------------
% Discrepancy
%--------------------------------------------------------------------------

% Remove the pacing spikes from the simulated and experimental data
idx_nospike = GetNoSpikeIdx( prediction_protocol, length( time_pr ) );
if only_core == 1
    idx_core = GetCoreOfProtocolIdx( prediction_protocol );
else
    idx_core = 1 : length( voltage_pr );
end
idx_d = intersect( idx_nospike, idx_core );

% variables for evaluating discrepancy model
voltage_pr_d = voltage_pr( idx_d );
time_pr_d = time_pr( idx_d );
variables_pr_d = variables_pr( idx_d, : );
I_pr_d = I_pr(idx_d);
rates_pr_d = rates_pr( idx_d, :) ;
fluxes_pr_d = fluxes_pr( idx_d, : );
exp_data_pr_d = exp_data_pr( idx_d );
sensitivity_pr_d = sensitivity_pr( idx_d, :);
dIdp_pr_d = dIdp_pr(idx_d,:);
dVdt_pr_d = dVdt_pr(idx_d, :);
discrepancy_pr_d = discrepancy_pr( idx_d );

% Build the data matrix used for prediction
[data_matrix_pr_d, factors_pr_d, names_pr_d ] = BuildDataMatrix( conf, variables_pr_d, I_pr_d, rates_pr_d, fluxes_pr_d, ...
                                                       sensitivity_pr_d, dIdp_pr_d, voltage_pr_d, dVdt_pr_d, time_pr_d );

if isa( model_d, 'LinearModel' )
    if any( strcmp( conf, 'const' ) )
        disc_predicted = predict( model_d, data_matrix_pr_d( :, 2 : end ) );
    else
        disc_predicted = predict( model_d, data_matrix_pr_d ); 
    end
else
    disc_predicted = data_matrix_pr_d*model_d;
    included_idxs_d = abs(model_d) > 0;
    num_included_d = sum( included_idxs_d );
    names_included_d = names_pr_d( included_idxs_d );
end
                                                   
%--------------------------------------------------------------------------
% Discrepancy in Open
%--------------------------------------------------------------------------
                                                   
% Remove the indices near the reversal potential for Odisc
idx_erev = find( abs( voltage_pr - erev ) > eps );
idx_od = intersect( idx_d, idx_erev );

voltage_pr_od = voltage_pr( idx_od );
time_pr_od = time_pr( idx_od );
variables_pr_od = variables_pr( idx_od, : );
I_pr_od = I_pr(idx_od);
rates_pr_od = rates_pr( idx_od, :) ;
fluxes_pr_od = fluxes_pr( idx_od, : );
exp_data_pr_od = exp_data_pr( idx_od );
sensitivity_pr_od = sensitivity_pr( idx_od, :);
dIdp_pr_od = dIdp_pr(idx_od,:);
dVdt_pr_od = dVdt_pr(idx_od, :);

discrepancy_pr_od = discrepancy_pr( idx_od ) ./ ( G * ( voltage_pr( idx_od ) - erev ) );

% Construct the discrepancy model

% Build the data matrix used for prediction
[data_matrix_pr_od, factors_pr_od, names_pr_od ] = BuildDataMatrix( conf, variables_pr_od, I_pr_od, rates_pr_od, fluxes_pr_od, ...
                                                        sensitivity_pr_od, dIdp_pr_od, voltage_pr_od, dVdt_pr_od, time_pr_od );

if isa( model_od, 'LinearModel' )
    if any( strcmp( conf, 'const' ) )
        odisc_predicted = predict( model_od, data_matrix_pr_od( :, 2 : end ) );
    else
        odisc_predicted = predict( model_od, data_matrix_pr_od ); 
    end
else
    odisc_predicted = data_matrix_pr_od*model_od;
    included_idxs_od = abs(model_od) > 0;
    num_included_od = sum( included_idxs_od );
    names_included_od = names_pr_od( included_idxs_od );
end

%% Save data and plot results
% some config for naming files
if strcmp( method, 'lasso' )
    methodName = 'LASSO';
elseif strcmp( method, 'stepwiselm' )
    methodName = 'StepwiseLM';
else
    methodName = method;
end

if only_core == 1
    dataTypeName = 'core';
else
    dataTypeName = 'full';
end

if strcmp( training_protocol, 'sine_wave' )
    trainName= 'SW';
elseif strcmp( training_protocol, 'ap' )
    trainName = 'AP';
elseif strcmp( training_protocol, 'maz_wang_div_diff' )
    trainName = 'MWDD';
elseif strcmp( training_protocol, 'max_diff' )
    trainName = 'MD';
elseif strcmp( training_protocol, 'equal_proportions' )
    trainName = 'EP';
elseif strcmp( training_protocol, 'original_sine' )
    trainName = 'OSW';
else
    trainName = 'TR';
end
    
if strcmp( prediction_protocol, 'sine_wave' )
    predName= 'SW';
elseif strcmp( prediction_protocol, 'ap' )
    predName = 'AP';
elseif strcmp( prediction_protocol, 'maz_wang_div_diff' )
    predName = 'MWDD';
elseif strcmp( prediction_protocol, 'max_diff' )
    predName = 'MD';
elseif strcmp( prediction_protocol, 'equal_proportions' )
    predName = 'EP';
elseif strcmp( prediction_protocol, 'original_sine' )
    predName = 'OSW';
else
    predName = 'PR';
end

current_time = clock;
current_time = fix( current_time );

dirname = [ 'FittedModels/FitModel_' method '_tr_' training_protocol '_pr_' prediction_protocol '_' date ...
                '_' num2str( current_time( 4 ) ) '-' num2str( current_time( 5 ) ) '-' num2str( current_time( 6 ) ) ];
mkdir( dirname )
cd( dirname )

% save the config files and setup
save( 'config.mat', 'conf', 'only_core', 'method', 'training_protocol', 'prediction_protocol' );

% save the model of discrepancy
save( 'model_d.mat', 'model_d', 'data_matrix_tr_d', 'data_matrix_pr_d' )
% save the model of open-discrepancy
save( 'model_od.mat', 'model_od', 'data_matrix_tr_od', 'data_matrix_pr_od' )

fig_currents = figure( 'Units', 'Normalized', 'OuterPosition', [ 0 0 1 1 ] );

% detemine axes
time_min = min( [ min(time_pr_d), min(time_tr_d), min(time_pr_od), min(time_pr_d) ] );
time_max = max( [ max(time_pr_d), max(time_tr_d), max(time_pr_od), max(time_pr_d) ] );

d_min = min( [ min(I_tr_d), min( exp_data_tr_d ), min( I_tr_d + disc_modelled ) ...
               min(I_pr_d), min( exp_data_pr_d ), min( I_pr_d + disc_predicted )] );
d_max = max( [ max(I_tr_d), max( exp_data_tr_d ), max( I_tr_d + disc_modelled ) ...
               max(I_pr_d), max( exp_data_pr_d ), max( I_pr_d + disc_predicted )] );
ylim_d = 0.5*( [ d_min d_max ] );

od_min = min( [ min(I_tr_od), min( exp_data_tr_od ), min( I_tr_od + G*odisc_modelled.*( voltage_tr_od - erev ) ) ...
               min(I_pr_od), min( exp_data_pr_od ), min( I_pr_od + G*odisc_predicted.*( voltage_pr_od - erev ) )] );
od_max = max( [ max(I_tr_od), max( exp_data_tr_od ), max( I_tr_od + G*odisc_modelled.*( voltage_tr_od - erev ) ) ...
               max(I_pr_od), max( exp_data_pr_od ), max( I_pr_od + G*odisc_predicted.*( voltage_pr_od - erev ) )] );
ylim_od = 0.5*( [ od_min od_max ] );

% make plot of current with discrepancy
subplot( 2,2,1 )
plot( time_tr_d, exp_data_tr_d, 'r' )
hold on
plot( time_tr_d, I_tr_d, 'b', 'LineWidth', 2 )
plot( time_tr_d, I_tr_d + disc_modelled, 'k', 'LineWidth', 2 ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ trainName, ' Training d' ], 'Interpreter', 'None' )
ylim(ylim_d);

disc_base_tr_d = sqrt( sum( ( exp_data_tr_d - I_tr_d ).^2 ) );
disc_model_tr_d = sqrt( sum( ( exp_data_tr_d - ( I_tr_d + disc_modelled ) ).^2 ) );

legend( { 'Iexp', [ 'Isim ' num2str( disc_base_tr_d ) ], [ 'Isim+d '  num2str( disc_model_tr_d ) ] })

subplot( 2,2,2 )
plot( time_tr_od, exp_data_tr_od, 'r' )
hold on
plot( time_tr_od, I_tr_od, 'b', 'LineWidth', 2 )
plot( time_tr_od, I_tr_od + G*odisc_modelled.*( voltage_tr_od - erev ), 'k', 'LineWidth', 2 ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ trainName, ' Training d{_o}' ], 'Interpreter', 'Tex' )
ylim(ylim_od);

disc_base_tr_od = sqrt( sum( ( exp_data_tr_od - I_tr_od ).^2 ) );
disc_model_tr_od = sqrt( sum( ( exp_data_tr_od - ( I_tr_od + G*odisc_modelled.*( voltage_tr_od - erev ) ) ).^2 ) );

legend( { 'Iexp', [ 'Isim ' num2str( disc_base_tr_od ) ], [ 'Isim+G.Od.(V-E) '  num2str( disc_model_tr_od ) ] })

subplot( 2,2,3 )
plot( time_pr_d, exp_data_pr_d, 'r' )
hold on
plot( time_pr_d, I_pr_d, 'b', 'LineWidth', 2 )
plot( time_pr_d, I_pr_d + disc_predicted, 'k', 'LineWidth', 2 )
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ predName, ' Prediction d' ], 'Interpreter', 'None' )
ylim(ylim_d);

disc_base_pr_d = sqrt( sum( ( exp_data_pr_d - I_pr_d ).^2 ) );
disc_model_pr_d = sqrt( sum( ( exp_data_pr_d - ( I_pr_d + disc_predicted ) ).^2 ) );

legend( { 'Iexp', [ 'Isim ' num2str( disc_base_pr_d ) ], [ 'Isim+d '  num2str( disc_model_pr_d ) ] })

subplot( 2,2,4 )
plot( time_pr_od, exp_data_pr_od, 'r' )
hold on
plot( time_pr_od, I_pr_od, 'b', 'LineWidth', 2 )
plot( time_pr_od, I_pr_od + G*odisc_predicted.*( voltage_pr_od - erev ), 'k', 'LineWidth', 2  ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ predName, ' prediction d{_o}' ], 'Interpreter', 'Tex' )
ylim(ylim_od);

disc_base_pr_od = sqrt( sum( ( exp_data_pr_od - I_pr_od ).^2 ) );
disc_model_pr_od = sqrt( sum( ( exp_data_pr_od - ( I_pr_od + G*odisc_predicted.*( voltage_pr_od - erev ) ) ).^2 ) );

legend( { 'Iexp', [ 'Isim ' num2str( disc_base_pr_od ) ], [ 'Isim+G.Od.(V-E) '  num2str( disc_model_pr_od ) ] })

% Beautify plot
set( fig_currents, 'Color', 'w' );
set(findall( fig_currents, 'type', 'axes'), 'Box', 'off' );
set(findall( fig_currents, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
set(findall( fig_currents, 'type', 'axes'), 'xLim', [ time_min time_max]);

export_fig( fig_currents, [ methodName '_' trainName '_' predName '_' dataTypeName '_currents.tif' ], '-tif' )
export_fig( fig_currents, [ methodName '_' trainName '_' predName '_' dataTypeName '_currents.png' ], '-png' )

%--------------------------------------------------------------------------
% make plot of discrepancy

d_min = min( [ min(discrepancy_tr_d), min( disc_modelled ) ...
               min( discrepancy_pr_d), min( disc_predicted )] );
d_max = max( [ max(discrepancy_tr_d), max( disc_modelled ) ...
               max(discrepancy_pr_d), max( disc_modelled )] );
ylim_d = 0.5*( [ d_min d_max ] );

od_min = min( [ min(discrepancy_tr_od), min( odisc_modelled ) ...
               min(discrepancy_pr_od), min( odisc_predicted )] );
od_max = max( [ max(discrepancy_tr_od), max( odisc_modelled ) ...
               max(discrepancy_pr_od), max(odisc_predicted )] );
ylim_od = ( [ od_min od_max ] );


subplot( 2,2,1 )
plot( time_tr_d, discrepancy_tr_d, 'b', 'LineWidth', 2  )
hold on
plot( time_tr_d, disc_modelled, 'k', 'LineWidth', 2  ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ trainName, ' Training d' ], 'Interpreter', 'None' )
ylim(ylim_d);

disc_base_tr_d = sqrt( sum( discrepancy_tr_d.^2 ) );
disc_model_tr_d = sqrt( sum( (discrepancy_tr_d-disc_modelled).^2 ) );

legend( { [ 'Iexp-Isim ' num2str( disc_base_tr_d ) ], [ 'Iexp-(Isim+d) '  num2str( disc_model_tr_d ) ] })

subplot( 2,2,2 )
plot( time_tr_od, discrepancy_tr_od, 'b', 'LineWidth', 2 )
hold on
plot( time_tr_od, odisc_modelled, 'k', 'LineWidth', 2  ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ trainName, ' training d{_o}' ], 'Interpreter', 'Tex' )
ylim(ylim_od);

disc_base_tr_od = sqrt( sum( discrepancy_tr_od.^2 ) );
disc_model_tr_od = sqrt( sum( (discrepancy_tr_od-odisc_modelled).^2 ) );

legend( { [ 'Odexp ' num2str( disc_base_tr_od ) ], [ 'Odexp-Od '  num2str( disc_model_tr_od ) ] })

subplot( 2,2,3 )
plot( time_pr_d, discrepancy_pr_d, 'b', 'LineWidth', 2  )
hold on
plot( time_pr_d, disc_predicted, 'k', 'LineWidth', 2  )
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ predName, ' prediction d' ], 'Interpreter', 'None' )
ylim(ylim_d);

disc_base_pr_d = sqrt( sum( discrepancy_pr_d.^2 ) );
disc_model_pr_d = sqrt( sum( (discrepancy_pr_d-disc_predicted).^2 ) );

legend( {[ 'Iexp-Isim ' num2str( disc_base_pr_d ) ], [ 'Iexp-(Isim+d) '  num2str( disc_model_pr_d ) ] })

subplot( 2,2,4 )
plot( time_pr_od, discrepancy_pr_od, 'b', 'LineWidth', 2  )
hold on
plot( time_pr_od, odisc_predicted, 'k', 'LineWidth', 2  ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Current (nA)' )
title( [ predName, ' Prediction d{_o}' ], 'Interpreter', 'Tex' )
ylim(ylim_od);

disc_base_pr_od = sqrt( sum( ( discrepancy_pr_od).^2 ) );
disc_model_pr_od = sqrt( sum( ( discrepancy_pr_od - odisc_predicted ).^2 ) );

legend( {[ 'Odexp ' num2str( disc_base_pr_od ) ], [ 'Odexp-Od '  num2str( disc_model_pr_od ) ] })

% Beautify plot
set( fig_currents, 'Color', 'w' );
set(findall( fig_currents, 'type', 'axes'), 'Box', 'off' );
set(findall( fig_currents, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
set(findall( fig_currents, 'type', 'axes'), 'xLim', [ time_min time_max]);

export_fig( fig_currents, [ methodName '_' trainName '_' predName '_' dataTypeName '_discrepancy.tif' ], '-tif' )
export_fig( fig_currents, [ methodName '_' trainName '_' predName '_' dataTypeName '_discrepancy.png' ], '-png' )

%--------------------------------------------------------------------------
if ~isempty( lassofig_d )
    figure(lassofig_d )
    title( 'Cross-validated MSE of Lasso fit: d')
    set( gcf, 'Color', 'w' );
    set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
    set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
    export_fig( lassofig_d, [ methodName '_' trainName '_' predName '_' dataTypeName '_lambda_d.tif' ], '-tif' )
    export_fig( lassofig_d, [ methodName '_' trainName '_' predName '_' dataTypeName '_lambda_d.png' ], '-png' )
end

if ~isempty( lassofig_od )
    figure(lassofig_od )
    title( 'Cross-validated MSE of Lasso fit: d{_o}')
    set( gcf, 'Color', 'w' );
    set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
    set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
    
    export_fig( lassofig_od, [ methodName '_' trainName '_' predName '_' dataTypeName '_lambda_od.tif' ], '-tif' )
    export_fig( lassofig_od, [ methodName '_' trainName '_' predName '_' dataTypeName '_lambda_od.png' ], '-png' )
end

%--------------------------------------------------------------------------

% Bar chart of coefficients

fig_coeffs = figure( 'Units', 'Normalized', 'OuterPosition', [ 0 0 1 1 ] );
subplot( 1,2,1)
if isa( model_d, 'LinearModel' )
    model_d_bar = model_d.Coefficients.Estimate;
    names_included_d = model_d.CoefficientNames;
else
    i = find( model_d == 0 );
    model_d_bar = model_d;
    model_d_bar( i ) = [];
end

num_included_d = length( model_d_bar );

bar( model_d_bar )
xtick=get(gca,'xtick'); 
xMax=max(xtick); 
xMin=min(xtick); 
newXTick=linspace(xMin,xMax,num_included_d); 
set(gca,'xtick', newXTick );
xticklabels( names_included_d )
xtickangle( 90 )
title( 'Coefficients d' )

subplot( 1,2,2)
if isa( model_d, 'LinearModel' )
    model_od_bar = model_od.Coefficients.Estimate;
    names_included_od = model_od.CoefficientNames;
else
    i = find( model_od == 0 );
    model_od_bar = model_od;
    model_od_bar( i ) = [];
end

num_included_od = length( model_od_bar );

bar( model_od_bar )
xtick=get(gca,'xtick'); 
xMax=max(xtick); 
xMin=min(xtick); 
newXTick=linspace(xMin,xMax,num_included_od); 
set(gca,'xtick', newXTick );
xticklabels( names_included_od )
xtickangle( 90 )
title( 'Coefficients d{_o}' )

set( gcf, 'Color', 'w' );
set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);

export_fig( fig_coeffs, [ methodName '_' trainName '_' predName '_' dataTypeName '_coefficients.tif' ], '-tif' )
export_fig( fig_coeffs, [ methodName '_' trainName '_' predName '_' dataTypeName '_coefficients.png' ], '-png' )

% Back to code
cd ../..

end