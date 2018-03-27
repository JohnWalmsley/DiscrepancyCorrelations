protocol = 'ap';
exp_ref = '16713110';
[ ~, I, ~, ~, ~, ~, ~, ~, time, ~ ]  ...
                            = CalculateVariables( protocol, exp_ref );
[ discrepancy, exp_data ] = CalculateDiscrepancy( exp_ref, protocol, I );

idx = GetNoSpikeIdx( protocol, length( time ) );

error = sum( discrepancy(idx).^2 );

figure; plot( time(idx), exp_data(idx), 'k' )
hold on 
plot( time, I, 'r' )

protocol = 'ap_lag';

[ ~, I, ~, ~, ~, ~, ~, ~, time, ~ ]  ...
                            = CalculateVariables( protocol, exp_ref );
[ discrepancy_lag, exp_data ] = CalculateDiscrepancy( exp_ref, protocol, I );

plot( time, I, 'b' )

xlabel( 'Time (ms)' )
ylabel( 'Voltage (mV)' )

error_lag = sum( discrepancy_lag(idx).^2 );

disp( [ 'Standard voltage error: ' num2str(error) ])
disp( [ 'Lag voltage error: ' num2str(error_lag) ])