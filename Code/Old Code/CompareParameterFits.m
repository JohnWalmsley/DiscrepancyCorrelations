% Compare parameter fits between AP and Sine Wave Protocols
exp_ref = '16713110';
model = 'hh';
hist_fig = figure;
fit_fig = figure;

% Load the parameters fitted to the sine wave protocol
protocol = 'sine_wave';
[ params_sine_wave, v_sine_wave ] = FindingBestFitsAfterMCMC( model, protocol, exp_ref );

chain_sine_wave = importdata( [ '../MCMCResults/MCMCChain_' exp_ref '_' model '_' protocol '_30082148.mat' ] );
likelihood_1=importdata( [ '../MCMCResults/MCMCLikelihood_' exp_ref '_' model '_' protocol '_30082148.mat' ] );
[~,v]=max(likelihood_1);
max_likelihood_sine_wave = chain_sine_wave(v,:);
% Discard first 50000 points (burn in)
chain_sine_wave=chain_sine_wave(50001:250000,:);
figure( hist_fig )
for i = 1 : size( chain_sine_wave, 2 )
    subplot(3,3,i)
    hist(chain_sine_wave(:,i),100)
%     h = histogram(chain_sine_wave(:,i),100);
%     set(h,'FaceColor','r','LineStyle','none');
%     hold on
%     plot(max_likelihood_sine_wave(i),0,'kx','MarkerSize',16,'LineWidth',4)
%     hold off
end
  
% load the parameters for the AP protocol
protocol = 'sine_wave_ap_protocol';
[ params_ap, v_ap ] = FindingBestFitsAfterMCMC( model, protocol, exp_ref );

chain_ap = importdata( [ '../MCMCResults/MCMCChain_' exp_ref '_' model '_' protocol '_29082349.mat' ] );
likelihood_1=importdata( [ '../MCMCResults/MCMCLikelihood_' exp_ref '_' model '_' protocol '_29082349.mat' ] );
[~,v]=max(likelihood_1);
max_likelihood_ap = chain_ap(v,:);
% Discard first 50000 points (burn in)
chain_ap=chain_ap(50001:250000,:);
figure( hist_fig )
for i = 1 : size( chain_ap, 2 )
    subplot(3,3,i)
    hold on
    %h = histogram(chain_ap(:,i),100);
    hist(chain_ap(:,i),100);
    %set(h,'FaceColor','b','LineStyle','none');
%     plot(max_likelihood_ap(i),0,'kx','MarkerSize',16,'LineWidth',4)
    hold off
%     axis('tight')
%     % Plots parameters with maximum likelihood from experimental MCMC chain
%     % (which corresponds to parameter set used to generate simulated data).
%     % first chain is ref and second is blue
%     set(gca, 'FontSize',20)
%     xlabel(['Parameter ',num2str(i)],'FontSize',20)
%     
%     if i==9
%         xlabel('Conductance','FontSize',20)
%     end
%     if i==1||i==4||i==7
%         ylabel(['Probability Density'],'FontSize',20)
%     end
end

% Plot the sine wave 

protocol = 'sine_wave';
exp_ref = '16713110';
fitting_protocol = 'sine_wave';
[ ~, I_sw_sw, ~, ~, ~, ~, ~, ~, ~, ~ ]=CalculateVariables( protocol, exp_ref, fitting_protocol );
[ discrepancy_sw_sw, ~ ] = CalculateDiscrepancy( exp_ref, protocol, I_sw_sw );
fitting_protocol = 'ap';
[ ~, I_sw_ap, ~, ~, ~, ~, ~, ~, ~, ~ ]=CalculateVariables( protocol, exp_ref, fitting_protocol );
[ discrepancy_sw_ap, ~ ] = CalculateDiscrepancy( exp_ref, protocol, I_sw_ap );
fitting_protocol = 'sine_wave_ap';
[ ~, I_sw_swap, ~, ~, ~, ~, ~, ~, time_sw, ~ ]=CalculateVariables( protocol, exp_ref, fitting_protocol );
[ discrepancy_sw_swap, exp_data_sw ] = CalculateDiscrepancy( exp_ref, protocol, I_sw_swap );

idx_sw = GetNoSpikeIdx( protocol, length(I_sw_ap) );
figure;
hold on
plot( time_sw( idx_sw ), exp_data_sw( idx_sw ), 'k' )
plot( time_sw( idx_sw ), I_sw_sw( idx_sw ), 'r', 'LineWidth', 2 )
plot( time_sw( idx_sw ), I_sw_ap( idx_sw ), 'b', 'LineWidth', 2 )
plot( time_sw( idx_sw ), I_sw_swap( idx_sw ), 'g', 'LineWidth', 2 )

title( 'Sine Wave Protocol: Current' )
xlabel( 'Time (ms)' );
ylabel( 'Current (nA)' );
legend( { 'Data', 'Fit to SW', 'Fit to AP', 'Fit to SW+AP' } );
savefig( gcf, 'Figures/SWprotocol_SWvsSWAPfit_Current.fig' );

figure;
hold on
plot( time_sw( idx_sw ), discrepancy_sw_sw( idx_sw ), 'r', 'LineWidth', 2 )
plot( time_sw( idx_sw ), discrepancy_sw_ap( idx_sw ), 'b', 'LineWidth', 2 )
plot( time_sw( idx_sw ), discrepancy_sw_swap( idx_sw ), 'g', 'LineWidth', 2 )
plot( time_sw, zeros(size(time_sw)), 'k--' )

title( 'Sine Wave Protocol: Discrepancy' )
xlabel( 'Time (ms)' );
ylabel( 'Discrepancy (nA)' );
legend( { 'Fit to SW', 'Fit to AP', 'Fit to SW+AP' } );
savefig( gcf, 'Figures/SWprotocol_SWvsSWAPfit_Discrepancy.fig' );

protocol = 'ap';
exp_ref = '16713110';
fitting_protocol = 'sine_wave';
[ ~, I_ap_sw, ~, ~, ~, ~, ~, ~, ~, ~ ]=CalculateVariables( protocol, exp_ref, fitting_protocol );
[ discrepancy_ap_sw, ~ ] = CalculateDiscrepancy( exp_ref, protocol, I_ap_sw );
fitting_protocol = 'ap';
[ ~, I_ap_ap, ~, ~, ~, ~, ~, ~, ~, ~ ]=CalculateVariables( protocol, exp_ref, fitting_protocol );
[ discrepancy_ap_ap, ~ ] = CalculateDiscrepancy( exp_ref, protocol, I_ap_ap );
fitting_protocol = 'sine_wave_ap';
[ ~, I_ap_swap, ~, ~, ~, ~, ~, ~, time_ap, ~ ]=CalculateVariables( protocol, exp_ref, fitting_protocol );
[ discrepancy_ap_swap, exp_data_ap ] = CalculateDiscrepancy( exp_ref, protocol, I_ap_swap );

idx_ap = GetNoSpikeIdx( protocol, length(I_ap_ap) );
figure;
hold on
plot( time_ap( idx_ap ), exp_data_ap( idx_ap ), 'k' )
plot( time_ap( idx_ap ), I_ap_sw( idx_ap ), 'r', 'LineWidth', 2 )
plot( time_ap( idx_ap ), I_ap_ap( idx_ap ), 'b', 'LineWidth', 2 )
plot( time_ap( idx_ap ), I_ap_swap( idx_ap ), 'g', 'LineWidth', 2 )
title( 'AP Protocol: Current' )
xlabel( 'Time (ms)' );
ylabel( 'Current (nA)' );
legend( { 'Data', 'Fit to SW', 'Fit to AP', 'Fit to AP+SW' } );
savefig( gcf, 'Figures/APprotocol_SWvsSWAPfit_Current.fig' );

figure;
hold on
plot( time_ap( idx_ap ), discrepancy_ap_sw( idx_ap ), 'r', 'LineWidth', 2 )
plot( time_ap( idx_ap ), discrepancy_ap_ap( idx_ap ), 'b', 'LineWidth', 2 )
plot( time_ap( idx_ap ), discrepancy_ap_swap( idx_ap ), 'g', 'LineWidth', 2 )
plot( time_ap, zeros(size(time_ap)), 'k--' )
title( 'AP Protocol: Discrepancy' )
xlabel( 'Time (ms)' );
ylabel( 'Discrepancy (nA)' );
legend( { 'Fit to SW', 'Fit to AP', 'Fit to AP+SW' } );
savefig( gcf, 'Figures/APprotocol_SWvsSWAPfit_Discrepancy.fig' );


disp( [ 'Sine Wave trace discrepancy sine wave fit: ' num2str( sqrt( sum( discrepancy_sw_sw( idx_sw ).^2 ) ) ) ] )
disp( [ 'Sine Wave trace discrepancy AP fit: ' num2str( sqrt( sum( discrepancy_sw_ap( idx_sw ).^2 ) ) ) ] )
disp( [ 'Sine Wave trace discrepancy sine wave + AP fit: ' num2str( sqrt( sum( discrepancy_sw_swap( idx_sw ).^2 ) ) ) ] )
disp( [ 'AP trace discrepancy sine wave fit: ' num2str( sqrt( sum( discrepancy_ap_sw( idx_sw ).^2 ) ) ) ] )
disp( [ 'AP trace discrepancy AP fit: ' num2str(sqrt( sum( discrepancy_ap_ap( idx_sw ).^2 ) )) ] )
disp( [ 'AP trace discrepancy sine wave + AP fit: ' num2str( sqrt( sum( discrepancy_ap_swap( idx_sw ).^2 ) ) ) ] )