function TestPredictDiscrepancyUsingLinearModel( exp_ref, training_protocol, method, test_idx, conf  )


%% Defining the discrepancy in the open probability:
% G * (O+d)*(V-V_E) = data
% Therefore data - Isim = G*d*(V-V_E)
% and so d = ( data-Isim ) / (G*(V-V_E)) outside V_E + eps

addpath( genpath( '../../SharedFunctions' ) )
addpath ../../../MathworksFileExchange/getLinearDependent/

%% Configuration

if nargin < 5
    conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'sensitivity', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise' };
    %conf = { 'const', 'variables', 'I', 'rates', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise', 'thin' };
    %conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise' };
    %conf = { 'const', 'variables' };
end

if nargin < 4
    generate = 1;
end
    

% eps determines voltage proximity to reversal when calculating discrepancy
% in open probability
eps = 5;
only_core = 0;

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
voltage_tr_d    = voltage_tr( idx_d );
time_tr_d       = time_tr( idx_d );
variables_tr_d  = variables_tr( idx_d, : );
I_tr_d          = I_tr(idx_d);
rates_tr_d      = rates_tr( idx_d, :) ;
fluxes_tr_d     = fluxes_tr( idx_d, : );
exp_data_tr_d   = exp_data_tr( idx_d );
sensitivity_tr_d = sensitivity_tr( idx_d, :);
dIdp_tr_d       = dIdp_tr(idx_d,:);
dVdt_tr_d       = dVdt_tr(idx_d, :);
discrepancy_tr_d = discrepancy_tr( idx_d );

% calculate rates fluxes sensitivity
[data_matrix_tr_d, factors_tr_d, names_tr_d ] = BuildDataMatrix( conf, variables_tr_d, I_tr_d, rates_tr_d, fluxes_tr_d, ...
                                                       sensitivity_tr_d, dIdp_tr_d, voltage_tr_d, dVdt_tr_d, time_tr_d );                                             
% Are any members of the data matrix linearly dependent on one another?
[ ~, ind_idx, grps ] = getLinearIndependent(data_matrix_tr_d( :, 2 : end ), 1); % leave out the constant term.
dependent_grps = ( cellfun( 'length', grps ) > 1 );
if any ( dependent_grps )
    disp( 'Warning: Some predictors are linearly dependent between classes!' )
    nonind_grp_idxs = find( dependent_grps );
    for i=nonind_grp_idxs
        dep_idxs = grps{ i };
        namestring = sprintf( '%s ', names_tr_d{ dep_idxs+1 } ); % +1 due to removing const from mx.
        disp( [ namestring 'are linearly dependent' ] );
    end
end

% add noise of the order of the experimental data
sigma = sqrt( var( exp_data_tr_d( 1 : 1000) ) );
% Need to get SNR about right:
num_var = length( names_tr_d );
if generate == 1
    num_var_test = randi( [ 1 10 ], 1 );
    test_idx = unique( randi( [ 1 num_var ], [ 1 num_var_test ]  ) );
    num_var_test = length( test_idx ); % correct in case two random numbers are equal.
    mag_test = range( discrepancy_tr_d ) / length( test_idx );
    num_var = length( names_tr_d );

    nonzero_idx = zeros( 1, num_var );
    nonzero_idx( test_idx ) = 1;
    
else
    num_var_test = length( test_idx );
    mag_test = range( discrepancy_tr_d ) / num_var_test;
    num_var = length( names_tr_d );
    
    nonzero_idx = zeros( 1, num_var );
    nonzero_idx( test_idx, num_var_test ) = 1;
end
names_original = names_tr_d( test_idx );

% Determine which traces to use for the discrepancy
discrepancy_tr_d = sum(bsxfun( @times, mag_test*nonzero_idx, data_matrix_tr_d ),2) + sigma*randn( size(discrepancy_tr_d) );                                                   
figure;plot( discrepancy_tr_d)
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
if ~isa( model_d, 'LinearModel' )
    included_idxs_d = abs(model_d) > 0;
    num_included_d = sum( included_idxs_d );
    names_included_d = names_tr_d( included_idxs_d );

else
    model_d.Coefficients.Estimate;
end

%% Figures

%--------------------------------------------------------------------------
if ~isempty( lassofig_d )
    figure(lassofig_d )
    title( 'Cross-validated MSE of Lasso fit: d')
    set( gcf, 'Color', 'w' );
    set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
    set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
    
    export_fig( lassofig_d, 'LassoMSE_d.tif', '-tif' )
end


%--------------------------------------------------------------------------

% Bar chart of coefficients

fig_coeffs = figure( 'Units', 'Normalized', 'OuterPosition', [ 0 0 1 1 ] );
subplot( 1,2,1)
if isa( model_d, 'LinearModel' )
    model_d_bar = model_d.Coefficients.Estimate;
    names_included_d = model_d.CoefficientNames;
    num_included_d = length( model_d_bar );
else
    i = find( model_d == 0 );
    model_d_bar = model_d;
    model_d_bar( i ) = [];
end

data_max = max( max( model_d_bar ), max( mag_test ) );
data_min = min( min( model_d_bar ), min( mag_test ) );

b=bar( 1 : length( model_d_bar ), model_d_bar );
width = get( b, 'BarWidth' );

%xtick=get(gca,'xtick'); 
%xMax=max(xtick); 
%xMin=min(xtick); 
%newXTick=linspace(xMin,xMax,num_included_d); 
newXTick = 1 : length( model_d_bar );
set(gca,'xtick', newXTick );
xticklabels( names_included_d )
xtickangle( 90 )
title( 'Predicted Coefficients d' )
ylim( [ data_min, data_max ] );

subplot( 1,2,2)
model_original_bar = mag_test*ones( 1, num_var_test );

bar( 1 : num_var_test, model_original_bar, 'BarWidth', width )
%xtick=get(gca,'xtick'); 
%xMax=max(xtick); 
%xMin=min(xtick); 
%newXTick=linspace(xMin,xMax,num_var_test); 
newXTick = 1 : num_var_test;
set(gca,'xtick', newXTick );

xticklabels( names_original )
xtickangle( 90 )
title( 'Original Coefficients d' )
ylim( [ data_min, data_max ] );

set( gcf, 'Color', 'w' );
set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);

export_fig( fig_coeffs, 'coeffs.tif', '-tif' )

% Back to code
%cd ../..

end