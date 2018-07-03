function PlotVoltageTraces( fitting_protocols, prediction_protocols )

addpath ../../SharedFunctions/
addpath ../../SharedFunctions/Models/
addpath ~/MathworksFileExchange/export_fig/

data_fig = figure( 'Units', 'Normalized', 'OuterPosition', [ 0 0 1 1 ] );

% Get number of protocols used in the fitting process
num_fit = length( fitting_protocols );
num_pred = length( prediction_protocols );

max_num_plots = max( num_fit, num_pred );

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

    % load voltage data
    V = importdata([ '../Protocols/' protocol '_protocol.mat']);
    time = 0: 0.1 : 0.1*(length(V)-1);
    
    % plot data
    subplot( 2, max_num_plots, pr )
    plot( time, V, 'k-', 'LineWidth', 2.0 );
    
    xlabel( 'Time (ms)' );
    ylabel( 'Membrane Potential (mV)' );
    title( protocol, 'Interpreter', 'None' )
    ylim( [ -140 60 ] )
    
    time_max = max( time_max, max(time));
    
end

for pr = 1 : num_pred
    protocol = prediction_protocols{ pr };
    
    % load voltage data
    V = importdata([ '../Protocols/' protocol '_protocol.mat']);
    time = 0: 0.1 : 0.1*(length(V)-1);

    % plot data
    subplot( 2, max_num_plots, max_num_plots + pr )
    plot( time, V, 'k-', 'LineWidth', 2.0 );
    
    xlabel( 'Time (ms)' );
    ylabel( 'Membrane Potential (mV)' );
    title( protocol, 'Interpreter', 'None' )
    ylim( [ -140 60 ] )
    
    time_max = max( time_max, max(time));
            
end

% Beautify plot
set( gcf, 'Color', 'w' );
set(findall( gcf, 'type', 'axes'), 'Box', 'off' );
set(findall( gcf, 'type', 'axes'), 'FontName','Arial','FontSize',12,'FontWeight','Bold',  'LineWidth', 2);
set(findall( gcf, 'type', 'axes'), 'xLim', [ 0 time_max]);

export_fig( gcf, [ 'Figures/PlotVoltage_' fitting_protocol '_' prediction_protocol '.tif' ], '-tif' )

end