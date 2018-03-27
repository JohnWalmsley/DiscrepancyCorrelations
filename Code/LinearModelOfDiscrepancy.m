function [ discrepancy, data_matrix, model, exp_data, I, time, factors ] = LinearModelOfDiscrepancy( fitting_protocol, exp_ref, method, conf )
%LINEARMODELOFDISCREPANCY Summary of this function goes here
%   Detailed explanation goes here

% Discrepancy = model*data

[ variables, I, rates, fluxes, sensitivity, dIdp, voltage, dVdt, time ] ...
                            = CalculateVariables( fitting_protocol, exp_ref );

[ discrepancy, exp_data ] = CalculateDiscrepancy( exp_ref, fitting_protocol, I );
% find pacing spikes
idx = GetNoSpikeIdx( fitting_protocol, length( I ) );
[ data_matrix, factors ] = BuildDataMatrix( conf, variables, I, rates, fluxes, sensitivity, dIdp, voltage, dVdt, time );
if (~any( strcmp( conf, 'const' ) ) && strcmp( method, 'regress' ) )
    disp( 'Warning: Method regess requires a constant vector, adding to data matrix' )
    data_matrix = [ ones( size( time ) ) data_matrix ];
end

data_matrix_corr = data_matrix( idx', : );

if strcmp( method, 'least-squares' )
    model = data_matrix_corr \ discrepancy( idx );
elseif strcmp( method, 'regress' )
    model = regress( discrepancy( idx ), data_matrix_corr );
elseif strcmp( method, 'stepwisefit' )
    if any( strcmp( conf, 'const' ) )
        data_matrix_corr_sw = data_matrix_corr( :, 2 : end );
    else
        disp( 'Warning: Method stepwise requires a constant vector, adding to data matrix')
        data_matrix_corr_sw = data_matrix_corr;
        data_matrix = [ ones( size( time) ) data_matrix ];
    end
    [ b,~,~,inmodel,stats,~,~ ] = stepwisefit( data_matrix_corr_sw, discrepancy( idx ) );
    intercept = stats.intercept;
    model = zeros( size( b ) );
    model( inmodel ) = b( inmodel );
    model = [ intercept; model ];
elseif strcmp( method, 'lasso' )
    % note - lasso does not take the constant vector as input:
    if any( strcmp( conf, 'const' ) )
        data_matrix_corr_lasso = data_matrix_corr( :, 2 : end );
    else
        disp( 'Warning: Method lasso requires a constant vector, adding to data matrix')
        data_matrix = [ ones( size( time) ) data_matrix ];
        data_matrix_corr_lasso = data_matrix_corr;
    end
    [ models, fitinfo ] = lasso( data_matrix_corr_lasso, discrepancy( idx ), 'CV', 10, 'DFMax', 20 );
    lassoPlot(models,fitinfo,'PlotType','Lambda','XScale','log');
    lassoPlot(models,fitinfo,'PlotType','CV');
    intercepts = fitinfo.Intercept;
    intercept = intercepts( fitinfo.Index1SE );    
    model = [ intercept; models( :, fitinfo.Index1SE ) ];
else
    disp('Method not valid, proceeding with Least Squares' )
    model = data_matrix_corr \ discrepancy( idx' );
end

end

