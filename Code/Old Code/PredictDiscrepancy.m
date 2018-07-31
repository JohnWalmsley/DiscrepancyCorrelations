function PredictDiscrepancy( protocol, exp_ref, method, model_config )
%PREDICTDISCREPANCY Summary of this function goes here
%   Detailed explanation goes here

%model_config = { 'variables', 'I', 'rates', 'fluxes', 'sensitivity', 'dIdp', 'voltage', 'dVdt' };
%model_config = { 'variables' };
%model_config = { 'variables', 'voltage' };

Names = GetNames( model_config ); 

[ ~, data_matrix, model, exp_data, I, time, factors ] = LinearModelOfDiscrepancy( 'sine_wave', exp_ref, method, model_config );
idx = GetNoSpikeIdx( 'sine_wave', length( I ) );
size(model)
% Plot the currents - sine wave
fig=figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] ); plot( time(idx), exp_data(idx), 'r' )
hold on
plot( time(idx), I(idx), 'b' )
plot( time(idx), I(idx) + data_matrix(idx,:)*model, 'k' );
hold off
xlabel( 'Time(ms' )
ylabel( 'Current (mv)')
legend('Exp. Data', 'Model', 'Corrected Model')
hold off
title( 'Sine Wave Protocol - Training' )
saveas( fig,'CurrentSineWave.fig' )
saveas( fig,'CurrentSineWave.tif' )

% Plot the discrepancies - sine wave
fig=figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] ); plot( time(idx), exp_data(idx)-I(idx), 'r')
hold on 
plot( time(idx), data_matrix(idx,:)*model, 'k', 'LineWidth', 2)
xlabel( 'Time(ms' )
ylabel( 'Discrepancy (mv)')
legend('Actual Discrepancy', 'Predicted Discrepancy')
hold off
title( 'Sine Wave Protocol - Training' )
saveas( fig,'DiscSineWave.fig' )
saveas( fig,'DiscSineWave.tif' )

contribution_matrix = data_matrix.*model';
fig=figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] ); plot( time(idx), contribution_matrix(idx,:) ) 
hold on
plot( time(idx), data_matrix(idx,:)*model, 'k', 'LineWidth', 4 ) 
hold off
title( 'Contributions to est. discrepancy, AP Protocol' )
saveas( fig,'DiscContributionsSineWave.fig' )
saveas( fig,'DiscContributionsSineWave.tif' )

limits = [ min(min(contribution_matrix( idx, : ))) max(max(contribution_matrix( idx, : ))) ];
fig=figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] );
for i = 1 : length( model )
    subplot( 4, ceil( length(model)/4), i )
    plot( time(idx), contribution_matrix( idx, i ) );
    ylim( limits )
    title(Names{i});
end
saveas( fig,'DiscContributionsByParamSineWave.fig' )
saveas( fig,'DiscContributionsByParamSineWave.tif' )

%% Now use the model to try and predict the discrepancy in the AP protocol

[ variables_new, I_new, rates_new, fluxes_new, sensitivity_new, dIdp_new, voltage_new, dVdt_new, time_new ] = CalculateVariables( protocol, exp_ref );
%data_matrix_new = BuildDataMatrix( model_config, variables_new, I_new, rates_new, fluxes_new, sensitivity_new, dIdp_new, voltage_new, dVdt_new );
data_matrix_new = BuildDataMatrix( model_config, variables_new, I_new, rates_new, fluxes_new, sensitivity_new, dIdp_new, voltage_new, dVdt_new, time_new, factors );
[ discrepancy_new, exp_data_new ] = CalculateDiscrepancy( exp_ref, protocol, I_new );

idx_ap = GetNoSpikeIdx( 'ap', length(I_new) );
% Plot the currents - AP
figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] ); plot( time_new(idx_ap), exp_data_new(idx_ap), 'r' )
hold on
plot( time_new(idx_ap), I_new(idx_ap), 'b' )
plot( time_new(idx_ap), I_new(idx_ap) + data_matrix_new(idx_ap,:)*model, 'k' );
xlabel( 'Time(ms' )
ylabel( 'Current (mv)')
legend('Exp. Data', 'Model', 'Corrected Model')
hold off
title( 'AP Protocol - Prediction' )
% Plot the discrepancies - AP
fig = figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] ); plot( time_new(idx_ap), discrepancy_new(idx_ap), 'r')
hold on 
plot( time_new(idx_ap), data_matrix_new(idx_ap,:)*model, 'k', 'LineWidth', 2)
xlabel( 'Time(ms' )
ylabel( 'Discrepancy (mv)')
legend('Actual Discrepancy', 'Predicted Discrepancy')
hold off
title( 'AP Protocol - Prediction' )
saveas( fig,'DiscrepancyAP.fig' )
saveas( fig,'DiscrepancyAP.tif' )

discrepancy_new_corrected = exp_data_new - ( I_new + data_matrix_new*model );

% Sum of least squares: data-model discrepancy
sls_new = sqrt( sum( discrepancy_new(idx_ap).^2 ) );
% Sum of least squares: data-corrected model discrepancy
sls_new_corrected = sqrt( sum( discrepancy_new_corrected(idx_ap).^2 ) );

disp( [ 'Old error: ' num2str( sls_new ) ] )
disp( [ 'New error: ' num2str( sls_new_corrected ) ] ) 

contribution_matrix_new = data_matrix_new.*model';
fig=figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] ); plot( time_new(idx_ap), contribution_matrix_new(idx_ap,:) ) 
hold on
plot( time_new(idx_ap), data_matrix_new(idx_ap,:)*model, 'k', 'LineWidth', 4 ) 
hold off
xlabel( 'Time (ms)' )
ylabel( 'Contribution (mv)' )
title( 'Contributions to est. discrepancy, AP Protocol' )
saveas( fig,'DiscContributionsAP.fig' )
saveas( fig,'DiscContributionsAP.tif' )

fig=figure( 'Units', 'Normalized', 'OuterPosition', [ 0.1 0.1 0.8 0.8 ] );
limits = [ min(min(contribution_matrix_new( idx_ap, : ))) max(max(contribution_matrix_new( idx_ap, : ))) ];
for i = 1 : length( model )
    subplot( 4, ceil( length(model)/4), i )
    plot( time_new(idx_ap), contribution_matrix_new( idx_ap, i ) );
    ylim( limits )
    title(Names{i});
end
saveas( fig,'DiscContributionsByParamAP.fig' )
saveas( fig,'DiscContributionsByParamAP.tif' )

end
