function [ discrepancy, data_matrix, model ] = ...
                        LinearModelOfDiscrepancyWithInput( method, conf, discrepancy, data_matrix )
%LINEARMODELOFDISCREPANCY Summary of this function goes here
%   Detailed explanation goes here

% Discrepancy = model*data
if (~any( strcmp( conf, 'const' ) ) && strcmp( method, 'regress' ) )
    disp( 'Warning: Method regess requires a constant vector, adding to data matrix' )
    data_matrix = [ ones( size( time ) ) data_matrix ];
end

data_matrix_corr = data_matrix;

if strcmp( method, 'least-squares' )
    model = data_matrix_corr \ discrepancy;
elseif strcmp( method, 'regress' )
    model = regress( discrepancy, data_matrix_corr );
elseif strcmp( method, 'stepwisefit' )
    if any( strcmp( conf, 'const' ) )
        data_matrix_corr_sw = data_matrix_corr( :, 2 : end );
    else
        disp( 'Warning: Method stepwise requires a constant vector, adding to data matrix')
        data_matrix_corr_sw = data_matrix_corr;
        data_matrix = [ ones( size( time) ) data_matrix ];
    end
    [ b,~,~,inmodel,stats,~,~ ] = stepwisefit( data_matrix_corr_sw, discrepancy, 'penter', 0.0005, 'premove', 0.01 );
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
    [ models, fitinfo ] = lasso( data_matrix_corr_lasso, discrepancy, 'CV', 10, 'DFmax', 10 );
    lassoPlot(models,fitinfo,'PlotType','CV');
    intercepts = fitinfo.Intercept;
    intercept = intercepts( fitinfo.Index1SE );    
    model = [ intercept; models( :, fitinfo.Index1SE ) ];
elseif strcmp( method, 'stepwiselm' )
    % note - strcmp does not take the constant vector as input:
    if any( strcmp( conf, 'const' ) )
        data_matrix_corr_swlm = data_matrix_corr( :, 2 : end );
    end
    model = stepwiselm( data_matrix_corr_swlm, discrepancy, 'quadratic', 'Criterion', 'bic', 'Verbose', 2 );
    model
else
    disp('Method not valid, proceeding with Least Squares' )
    model = data_matrix_corr \ discrepancy;
end

end

