fig = figure( 'Units', 'Normalized', 'OuterPosition', [ 0 0 1 1 ] );

% sine wave
% Protocol
subplot( 5, 2, 1 )
V=importdata( '../Protocols/sine_wave_protocol.mat' );
time = 0 : 0.1:(length( V )-1)/10;
plot( time, V, 'LineWidth', 2 );
sw_core_idx = GetCoreOfProtocolIdx( 'sine_wave' );
hold on
plot( time( sw_core_idx ), V(sw_core_idx ), 'r-', 'LineWidth', 2 );
hold off
xlim( [ 0 9500 ] );
ylabel( 'Voltage (mV)' )
title( 'SW Protocol' )
%Data
subplot( 5, 2, 2 )
E=importdata( '../ExperimentalData/16713110/sine_wave_16713110_dofetilide_subtracted_leak_subtracted.mat');
sw_nospike_idx = GetNoSpikeIdx( 'sine_wave', length( E ) );
plot( time( sw_nospike_idx ), E( sw_nospike_idx ), 'LineWidth', 2 );
hold on
sw_core_idx_nospike = intersect( sw_nospike_idx, sw_core_idx );
plot( time( sw_core_idx_nospike ), E(sw_core_idx_nospike ), 'r', 'LineWidth', 2 );
hold off
xlim( [ 0 9500 ] );
ylim([-2 2 ] )
ylabel( 'Current (nA)' )
title( 'SW Current' )
% Action potential
% Protocol
subplot( 5, 2, 3 )
V=importdata( '../Protocols/ap_protocol.mat' );
time = 0 : 0.1:(length( V )-1)/10;
plot( time, V, 'LineWidth', 2 );
ap_core_idx = GetCoreOfProtocolIdx( 'ap' );
hold on
plot( time( ap_core_idx ), V(ap_core_idx ), 'r-', 'LineWidth', 2 );
hold off
xlim( [ 0 9500 ] );
ylabel( 'Voltage (mV)' )
title( 'AP Protocol' )
% Data
subplot( 5, 2, 4 )
E=importdata( '../ExperimentalData/16713110/ap_16713110_dofetilide_subtracted_leak_subtracted.mat');
ap_nospike_idx = GetNoSpikeIdx( 'ap', length( E ) );
plot( time( ap_nospike_idx ), E( ap_nospike_idx ), 'LineWidth', 2 );
hold on
ap_core_idx_nospike = intersect( ap_nospike_idx, ap_core_idx );
plot( time( ap_core_idx_nospike ), E(ap_core_idx_nospike ), 'r', 'LineWidth', 2 );
hold off
xlim( [ 0 9500 ] );
ylim([-2 2 ] )
ylabel( 'Current (nA)' )
title( 'AP Current' )
% original sine wave
% Protocol
subplot( 5, 2, 5 )
V=importdata( '../Protocols/original_sine_protocol.mat' );
time = 0 : 0.1:(length( V )-1)/10;
plot( time, V, 'LineWidth', 2 );
xlim( [ 0 9500 ] );
ylim([-200 150])
ylabel( 'Voltage (mV)' )
title( 'OSW Protocol' )
%Data
subplot( 5, 2, 6 )
E=importdata( '../ExperimentalData/16713110/original_sine_16713110_dofetilide_subtracted_leak_subtracted.mat');
osw_nospike_idx = GetNoSpikeIdx( 'original_sine', length( E ) );
plot( time( osw_nospike_idx ), E( osw_nospike_idx ), 'LineWidth', 2 );
xlim( [ 0 9500 ] );
ylim([-2 5 ] )
ylabel( 'Current (nA)' )
title( 'AP Current' )
% Equal Proportions
% Protocol
subplot( 5, 2, 7 )
V=importdata( '../Protocols/equal_proportions_protocol.mat' );
time = 0 : 0.1:(length( V )-1)/10;
plot( time, V, 'LineWidth', 2 );
xlim( [ 0 9500 ] );
ylabel( 'Voltage (mV)' )
title( 'EP Protocol' )
% Data
subplot( 5, 2, 8 )
E=importdata( '../ExperimentalData/16713110/equal_proportions_16713110_dofetilide_subtracted_leak_subtracted.mat');
ep_nospike_idx = GetNoSpikeIdx( 'equal_proportions', length( E ) );
plot( time( ep_nospike_idx ), E( ep_nospike_idx ), 'LineWidth', 2 );
xlim( [ 0 9500 ] );
ylim([-2 2 ] )
ylabel( 'Current (nA)' )
title( 'EP Current' )
% Maz Wang Div Diff
% Protocol
subplot( 5, 2, 9 )
V=importdata( '../Protocols/maz_wang_div_diff_protocol.mat' );
time = 0 : 0.1:(length( V )-1)/10;
plot( time, V, 'LineWidth', 2 );
xlim( [ 0 9500 ] );
ylabel( 'Voltage (mV)' )
xlabel( 'Time (ms)' )
title( 'MWDD Protocol' )
%Data
subplot( 5, 2, 10 )
E=importdata( '../ExperimentalData/16713110/maz_wang_div_diff_16713110_dofetilide_subtracted_leak_subtracted.mat');
mwdd_nospike_idx = GetNoSpikeIdx( 'maz_wang_div_diff', length( E ) );
plot( time( mwdd_nospike_idx ), E( mwdd_nospike_idx ), 'LineWidth', 2 );
xlim( [ 0 9500 ] );
ylim([-2 2 ] )
ylabel( 'Current (nA)' )
xlabel( 'Time (ms)' )
title( 'MWDD Current' )

set( gcf, 'Color', 'w' );
set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);

export_fig( fig, 'Protocols.png', '-png' )