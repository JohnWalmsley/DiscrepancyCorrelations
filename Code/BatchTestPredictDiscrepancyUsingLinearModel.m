function BatchTestPredictDiscrepancyUsingLinearModel( exp_ref, training_protocol, num_repeats, noise, conf  )


%% Defining the discrepancy in the open probability:
% G * (O+d)*(V-V_E) = data
% Therefore data - Isim = G*d*(V-V_E)
% and so d = ( data-Isim ) / (G*(V-V_E)) outside V_E + eps

addpath( genpath( '../../SharedFunctions' ) )
addpath ../../../MathworksFileExchange/getLinearDependent/

if isa( training_protocol, 'char' )
    tp_char = training_protocol;
    training_protocol = {};
    training_protocol{ 1 } = tp_char;
end

if length( training_protocol ) > 1
    disp( 'Note: Only first protocol used to fit parameters will be simulated' )
end

%% Configuration

if nargin < 3
    num_repeats = 10;
end

if nargin < 4
    noise = 1;
end

if nargin < 5
    conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'sensitivity', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise' };
    %conf = { 'const', 'variables', 'I', 'rates', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise', 'thin' };
    %conf = { 'const', 'variables', 'I', 'rates', 'fluxes', 'dIdp', 'voltage', 'dVdt', 'time', 'normalise' };
    %conf = { 'const', 'variables' };
end
   
% eps determines voltage proximity to reversal when calculating discrepancy
% in open probability
eps = 5;
only_core = 0;

%% Training Data

% Calculate variables and discrepancies: 
[ variables_tr, I_tr, rates_tr, fluxes_tr, sensitivity_tr, dIdp_tr, ~, dVdt_tr, time_tr, G ]  ...
                            = CalculateVariables( training_protocol{1}, exp_ref, training_protocol{1} );
[ discrepancy_tr, exp_data_tr ] = CalculateDiscrepancy( exp_ref, training_protocol{1}, I_tr );
voltage_tr=importdata(['../Protocols/' training_protocol{1} '_protocol.mat']);
size(variables_tr)
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
names_tr_d
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
if noise
    sigma = sqrt( var( exp_data_tr_d( 1 : 1000) ) );
else
    sigma = 0;
end

% Generate the test idxs:
num_var = length( names_tr_d );
nonzero_mx = zeros( num_repeats, num_var );
mag_test_vec = zeros( num_repeats, 1 );
for it = 1 : num_repeats
    num_var_test = randi( [ 1 20 ], 1 );
    test_idx = unique( randi( [ 1 num_var ], [ 1 num_var_test ]  ) );
    nonzero_mx( it, test_idx ) = 1;
    mag_test_vec( it ) = range( discrepancy_tr_d ) / length( test_idx ); % means SNR is about right.
end

% now loop over the number of entries
numPredictors = zeros( num_repeats, 1 );
MSE_lasso = zeros( num_repeats, 1 );
MSE_slm = zeros( num_repeats, 1 );

for it = 1 :num_repeats
    nonzero_idx = nonzero_mx( it, : );
    idxs = find( nonzero_idx );
    names_original = names_tr_d( idxs );
    names_original
    mag_test = mag_test_vec( it );
    numPredictors( it ) = sum( nonzero_idx );
    % Determine which traces to use for the discrepancy
    discrepancy_tr_d = sum(bsxfun( @times, mag_test*nonzero_idx, data_matrix_tr_d ),2) + sigma*randn( size(discrepancy_tr_d) );                                                   

    [ ~, ~, model_d_lasso, lassofig_d ] = ...
                        LinearModelOfDiscrepancyWithInput( 'lasso', conf, discrepancy_tr_d, data_matrix_tr_d, names_tr_d );
    close( lassofig_d ); % stops generation of loads of figures
    MSE_lasso( it ) = sqrt ( sum( ( model_d_lasso' - ( mag_test * nonzero_idx ) ).^2 ) );
    
    [ ~, ~, model_d_slm, ~ ] = ...
                        LinearModelOfDiscrepancyWithInput( 'stepwiselm', conf, discrepancy_tr_d, data_matrix_tr_d, names_tr_d );

    model_d_slm_coeffs = model_d_slm.Coefficients.Estimate;
    names_included_d_slm = model_d_slm.CoefficientNames;
    names_included_d_slm
    % loop over the names
    MSE_slm_sum = 0;
    for nidx = 1 : length( names_included_d_slm )
        name = names_included_d_slm(nidx);
        if strcmp( name, '(Intercept)' ) && any( strcmp( 'const', names_original ) )
            % A constant is in the model (has differnet name)
            MSE_slm_sum = MSE_slm_sum + ( mag_test - model_d_slm_coeffs( nidx ) )^2;
        elseif any( strcmp( name, names_original ) ) % variable is in original model
            MSE_slm_sum = MSE_slm_sum + ( mag_test - model_d_slm_coeffs( nidx ) )^2;
            % Value of coefficient is always mag_test
        else % variable is not in original model therefore origianl coeff is zero
            MSE_slm_sum = MSE_slm_sum + ( model_d_slm_coeffs( nidx ) )^2;
        end
    end
    MSE_slm( it ) = sqrt ( MSE_slm_sum );
    
end

%% Figures

%--------------------------------------------------------------------------

max_error = max( [ max( MSE_lasso ) max( MSE_slm ) ] );

fig=figure; 
subplot( 1,2,1)
title( 'Error (lasso)' ) 
plot( numPredictors, MSE_lasso, 'k.', 'MarkerSize', 10 );
p = polyfit( numPredictors, MSE_lasso, 1 );
hold on
plot( numPredictors, p(1)*numPredictors + p(2), 'r--', 'LineWidth', 2 )
hold off
xlabel( 'Number of predictors' )
ylabel( 'MSE in stat model coeffs' )
ylim( [ 0 max_error ] );
title( 'LASSO' )

subplot( 1,2,2)
plot( numPredictors, MSE_slm, 'k.', 'MarkerSize', 10 );
p = polyfit( numPredictors, MSE_slm, 1 );
hold on
plot( numPredictors, p(1)*numPredictors + p(2), 'r--', 'LineWidth', 2 )
hold off
xlabel( 'Number of predictors' )
ylabel( 'MSE in stat model coeffs' )
ylim( [ 0 max_error ] );
title( 'StepwiseLM' )

set( gcf, 'Color', 'w' );
set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);

export_fig( 'TrainingVsPredictionBatch.tif', '-tif')

% Back to code
%cd ../..

end