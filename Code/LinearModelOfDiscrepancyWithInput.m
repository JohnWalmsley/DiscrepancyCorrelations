function [ discrepancy, data_matrix, model, fig ] = ...
                        LinearModelOfDiscrepancyWithInput( method, conf, discrepancy, data_matrix )
%LINEARMODELOFDISCREPANCY Summary of this function goes here
%   Detailed explanation goes here

fig = [];

% Discrepancy = model*data
if (~any( strcmp( conf, 'const' ) ) && strcmp( method, 'regress' ) )
    disp( 'Warning: Method regess requires a constant vector, adding to data matrix' )
    data_matrix = [ ones( size( time ) ) data_matrix ];
end

if any( strcmp( conf, 'thin' ) )
    num_rows = size( data_matrix, 1 );
    idx = 1 : 10 : num_rows;
    data_matrix = data_matrix(idx,:);
    discrepancy = discrepancy( idx );
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
        data_matrix = [ ones( size( data_matrix, 1), 1 ) data_matrix ];
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
        data_matrix = [ ones( size( data_matrix, 1), 1 ) data_matrix ];
        data_matrix_corr_lasso = data_matrix_corr;
    end
    % Build a fixed partition of the data

    % NOTE: THIS WILL ONLY WORK IF /usr/local/MATLAB/R2017a/toolbox/stats/stats/+internal/+stats/cvpartitionImpl.m 
    % has been edited so that indices are no longer 'protected', but instead are 'public'!
    numPartitions = 10;
    partition = cvpartition( discrepancy, 'K', numPartitions );
    numDataPoints = length( discrepancy );
    partitionIdxs = ones( size(discrepancy ) );
    numPerGroup = ceil( length( discrepancy ) / numPartitions );
    numOver = numPerGroup * numPartitions - numDataPoints;
    
    % Sort into groups of size as close as possible to
    % #discrepancy/numPartitions
    for k = 1 : numPartitions - numOver
        partitionIdxs( (k-1)* numPerGroup + 1 : k*numPerGroup ) = k * ones( numPerGroup, 1 );  
    end
    for k = numPartitions - numOver + 1 : numPartitions
        init_idx = (numPartitions - numOver)*numPerGroup;
        npg = numPerGroup-1;
        k1 = k-(numPartitions - numOver);
        partitionIdxs(  init_idx + (k1-1)*npg + 1 : init_idx+k1*npg ) = k * ones( numPerGroup-1, 1 );  
    end
    partition.Impl.indices = partitionIdxs;
       
    size(partitionIdxs)
    size( discrepancy)
    
	[ models, fitinfo ] = lasso( data_matrix_corr_lasso, discrepancy, 'CV', partition );
    %[ models, fitinfo ] = lasso( data_matrix_corr_lasso, discrepancy, 'CV', 'resubstitution' );

    [ ~, fig ] = lassoPlot(models,fitinfo,'PlotType','CV');

    intercepts = fitinfo.Intercept;
    intercept = intercepts( fitinfo.Index1SE );    
    model = [ intercept; models( :, fitinfo.Index1SE ) ];
    
elseif strcmp( method, 'stepwiselm' )
    % note - strcmp does not take the constant vector as input:
    if any( strcmp( conf, 'const' ) )
        data_matrix_corr_swlm = data_matrix_corr( :, 2 : end );
    end
    [~, lin_indep_idx, ~] = getLinearIndependent( data_matrix_corr_swlm, 1 );
    exclude_idx = 1 : size( data_matrix_corr_swlm, 2 );
    exclude_idx(lin_indep_idx) = []; % just leaves the dependent idxs.
    model = stepwiselm( data_matrix_corr_swlm, discrepancy, 'constant', 'Criterion', 'aic', 'Upper', 'linear', 'Verbose', 1);% 'exclude', exclude_idx );
else
    disp('Method not valid, proceeding with Least Squares' )
    model = data_matrix_corr \ discrepancy;
end

end

