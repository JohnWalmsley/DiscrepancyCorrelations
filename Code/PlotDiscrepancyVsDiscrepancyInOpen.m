function PlotDiscrepancyVsDiscrepancyInOpen( exp_ref, fitting_protocols, prediction_protocols )

addpath ../../SharedFunctions/
addpath ../../SharedFunctions/Models/
addpath ~/MathworksFileExchange/export_fig/

model = 'hh';
seed = '30082148';

F = 96485;
R = 8314;
K_i = 130;
k_o = 4;
T = GetTemperature( exp_ref )+273.5;
erev = ((R*T)/F)*log(k_o/K_i);

eps = 5;

data_fig = figure( 'Units', 'Normalized', 'OuterPosition', [ 0 0 1 1 ] );

% Get number of protocols used in the fitting process
num_fit = length( fitting_protocols );
num_pred = length( prediction_protocols );

max_num_plots = 2*max( num_fit, num_pred );

% load parameters for fit
fitting_protocol = fitting_protocols{ 1 };
l = 1;
while l < num_fit
    l = l+1;
    fitting_protocol = [ fitting_protocol '_' fitting_protocols{ l } ];
end

prediction_protocol = prediction_protocols{ 1 };
m = 1;
while m < num_pred
    m = m+1;
    prediction_protocol = [ prediction_protocol '_' prediction_protocols{ m } ];
end

time_max = 0;
for pr = 1 : num_fit
    
    protocol = fitting_protocols{ pr };
    % Simulate data
    [ ~, I, ~, ~, ~, ~, V, ~, time, G ] = CalculateVariables( protocol, exp_ref, fitting_protocol );

    % load experimental data
    I_exp = importdata([ '../ExperimentalData/' exp_ref '/' protocol '_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
        
    % get indices for spike removal
    idx = GetNoSpikeIdx( protocol, length( time ) );
    
    % plot discrepancy
    subplot( 2, max_num_plots, 2*pr-1 )
    
    disc = I_exp-I;
    plot( time( idx ), disc( idx ), 'k.' );
    
    xlabel( 'Time (ms)' );
    ylabel( 'Current (nA)' );
    title( [protocol ' discrepancy'], 'Interpreter', 'None' )
    ylim([-0.5 0.5])
   
    sse_disc = sqrt( sum( disc.^2));

    leg = legend( [ 'Error: ' num2str( sse_disc) ] );
    set( leg, 'Box', 'Off' );
    
    % plot discrepancy in open

    near_erev_idx = find( abs( V-erev) < eps );
    idx1 = setdiff( idx, near_erev_idx );
    odisc = disc(idx1) ./ ( G * ( V( idx1 )- erev ) );
    
    subplot( 2, max_num_plots, 2*pr )
    plot( time( idx1 ), odisc, 'k.' );
    
    xlabel( 'Time (ms)' );
    ylabel( 'P (-)' );
    title( [protocol ' O-discrepancy'], 'Interpreter', 'None' )
    ylim([-0.1 0.1])
       
    sse_odisc = sqrt( sum( odisc.^2));
    
    leg = legend( [ 'Error: ' num2str( sse_odisc) ] );
    set( leg, 'Box', 'Off' );
    
    time_max = max(time_max,max(time));
        
end

for pr = 1 : num_pred
    protocol = prediction_protocols{ pr };
    % Simulate data
    [ ~, I, ~, ~, ~, ~, V, ~, time, G ] = CalculateVariables( protocol, exp_ref, fitting_protocol );
    
    % load experimental data
    I_exp = importdata([ '../ExperimentalData/' exp_ref '/' protocol '_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
        
    % get indices for spike removal
    idx = GetNoSpikeIdx( protocol, length( time ) );
    
    % plot data
    disc = I_exp-I;
    
    subplot( 2, max_num_plots, max_num_plots + 2*pr-1 )
    plot( time( idx ), disc( idx ), 'k.' );
    
    xlabel( 'Time (ms)' );
    ylabel( 'Current (nA)' );
    title( [ protocol ' discrepancy' ], 'Interpreter', 'None' )
    ylim([-0.5 0.5 ])
    
    sse_disc = sqrt( sum( disc.^2 ) );

    leg = legend( [ 'Error: ' num2str( sse_disc) ] );
    set( leg, 'Box', 'Off' );
    
    % plot data
    
    near_erev_idx = find( abs( V-erev) < eps );
    idx1 = setdiff( idx, near_erev_idx );
    odisc = disc(idx1) ./ ( G * ( V( idx1 )- erev ) );
    
    subplot( 2, max_num_plots, max_num_plots+2*pr )
    plot( time( idx1 ), odisc, 'k.' );
    
    xlabel( 'Time (ms)' );
    ylabel( 'Current (nA)' );
    title( [protocol ' O-discrepancy'], 'Interpreter', 'None' )
    ylim([-0.1 0.1])
    
    sse_odisc = sqrt( sum( odisc.^2));
    
    leg = legend( [ 'Error: ' num2str( sse_odisc) ] );
    set( leg, 'Box', 'Off' );
    
    time_max = max(time_max,max(time));

end

% Beautify plot
set( gcf, 'Color', 'w' );
set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
set(findall( gcf, 'type', 'axes'), 'xLim', [ 0 time_max]);

export_fig( gcf, [ 'Figures/PlotDiscrepancyVsDiscrepancyInOpen_' fitting_protocol '_' prediction_protocol '.tif' ], '-tif' )

end