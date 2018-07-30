function CompareParameterFits_SummaryPlots()

%% Setup
close all

% Generate data

model='hh';
protocol_list={'sine_wave', 'ap'};
plot_strings ={ 'b-', 'g-' };
exp_ref='16713110';
temperature=21.4; % temperature for cell 5 recordings

% loop over protocols used for fitting:
fig_protocols = figure;
fig_summary = figure;
for prot = 1 : length( protocol_list )
    protocol = protocol_list{ prot };
    protocol
    %% Inactivation protocol
    
    %Inactivation protocol
    %[~,model_type] =modeldata(model);
    
    %Import experimental data
    cd ../ExperimentalData
    cd(exp_ref)
    D=importdata(['inactivation_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
    cd ..
    cd ..
    l=length(D)/16;
    % Reformat experimental data to plot in traditional voltage steps
    inactivation_experimental_data(:,1)=[0:0.0001:-0.0001+l/10000];
    for i=1:16
        inactivation_experimental_data(:,i+1) = D(1+(i-1)*l:i*l);
        
    end
    cd Code
    [chain,likelihood] = FindingBestFitsAfterMCMC(model,protocol,exp_ref);
    
    [i,v]= max(likelihood);
    P= chain(v,:);
        
    % Import inactivation protocol
    cd ../Protocols
    V=importdata('inactivation_protocol.mat');
    cd ..
    cd Code
    % Simulate inactivation protocol with the parameters identified as providing best fit to sine wave protocol
    I=SimulatingData(35,'inactivation',P,V,temperature);
    inactivation_protocol_data(:,1)=[0:0.0001:-0.0001+l/10000];
    inactivation_prediction_data(:,1)=[0:0.0001:-0.0001+l/10000];
    for i=1:16
        inactivation_protocol_data(:,i+1) = V(1+(i-1)*l:i*l);
        inactivation_prediction_data(:,i+1) = I(1+(i-1)*l:i*l);     
    end
    % Plot results
    figure(fig_protocols);
    if prot == 1
        subplot( 3,3,1 )
        plot( inactivation_protocol_data(:,1),inactivation_protocol_data(:,2:end), 'k-' )
        title( 'Inactivation' )
        xlabel( 'Time (s)' )
        ylabel( 'Voltage (mV)' )
        ylim( [-125 65] )
        xlim( [ 0 3 ] )
        
        subplot( 3,3,4 )
        plot( inactivation_experimental_data(:,1),inactivation_experimental_data(:,2:end), 'r-' )
        xlabel( 'Time (s)' )
        ylabel( 'Exp. Current (nA)' )
        ylim( [-4 10 ] )
        xlim( [ 1.2 1.5 ] )
    end
    subplot( 3,3,7)
    hold on
    plot( inactivation_prediction_data(:,1),inactivation_prediction_data(:,2:end), plot_strings{prot} )
    hold off
    if prot == 2
        xlabel( 'Time (s)' )
        ylabel( 'Model Current (nA)' )
        ylim( [-4 10 ] )
        xlim( [ 1.2 1.5 ] )
    end
    
    D=[];
    I=[];
    
    %% Deactivation Protocol
    %[~,model_type] =modeldata(model);
    % Import experimental data
    cd ../ExperimentalData
    cd(exp_ref)
    D=importdata(['deactivation_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
    cd ..
    cd ..
    
    % Reformat experimental data to plot in traditional voltage steps
    l=length(D)/9;
    deactivation_experimental_data(:,1)=[0:0.0001:-0.0001+l/10000];
    deactivation_protocol_data(:,1)=[0:0.0001:-0.0001+l/10000];
    deactivation_prediction_data(:,1)=[0:0.0001:-0.0001+l/10000];
    for i=1:9
        
        deactivation_experimental_data(:,i+1) = D(1+(i-1)*l:i*l);
        
    end
    
    %Import deactivation protocol
    cd Protocols
    V=importdata('deactivation_protocol.mat');
    cd ..
    cd Code
    % Simulate deactivation protocol with parameters which provided best fit to sine wave
    I=SimulatingData(35,'deactivation',P,V,temperature);
    
    for i=1:9
        
        deactivation_protocol_data(:,i+1) = V(1+(i-1)*l:i*l);
        deactivation_prediction_data(:,i+1) = I(1+(i-1)*l:i*l);
        
    end
    
    % Plot data
        figure(fig_protocols);
    if prot == 1
        subplot( 3,3,2 )
        plot( deactivation_protocol_data(:,1),deactivation_protocol_data(:,2:end), 'k-' )
        title( 'Deactivation' )
        ylim( [-125 65] )
        xlim( [ 0 10 ] )
        subplot( 3,3,5 )
        plot( deactivation_experimental_data(:,1),deactivation_experimental_data(:,2:end), 'r-' )
        ylim( [-3.5 2 ] )
        xlim( [ 2.5 8 ] )
    end
    subplot( 3,3,8)
    hold on
    plot( deactivation_prediction_data(:,1),deactivation_prediction_data(:,2:end),plot_strings{prot} )
    hold off
    if prot == 2
        ylim( [-3.5 2 ] )
        xlim( [ 2.5 8 ] )   
        xlabel( 'Time (s)' )
    end
    D=[];
    I=[];
        
    %% Steady state activation protocol
    % Import experimental data
    cd ..
    cd ExperimentalData
    cd(exp_ref)
    D=importdata(['steady_activation_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
    cd ..
    cd ..
    
    % Reformat experimental data to plot in traditional voltage steps
    l=length(D)/7;
    steady_activation_experimental_data(:,1)=[0:0.0001:-0.0001+l/10000];
    steady_activation_protocol_data(:,1)=[0:0.0001:-0.0001+l/10000];
    steady_activation_prediction_data(:,1)=[0:0.0001:-0.0001+l/10000];
    for i=1:7
        
        steady_activation_experimental_data(:,i+1) = D(1+(i-1)*l:i*l);
    end
    
    % Import protocol
    cd Protocols
    V=importdata('steady_activation_protocol.mat');
    cd ..
    cd Code
    % Simulate steady state activation protocol with parameters which were identified as providing best fit to sine wave
    I=SimulatingData(35,'steady_activation',P,V,temperature);
    
    for i=1:7
        
        steady_activation_protocol_data(:,i+1) = V(1+(i-1)*l:i*l);
        steady_activation_prediction_data(:,i+1) = I(1+(i-1)*l:i*l);
    end
    
    % Plot data
    figure(fig_protocols);
    if prot == 1
        
        subplot( 3,3,3 )
        plot( steady_activation_protocol_data(:,1),steady_activation_protocol_data(:,2:end), 'k-' )
        title( 'Steady Activation' )
        ylim( [-125 65] )
        xlim( [ 0 8 ] )
        subplot( 3,3,6 )
        plot( steady_activation_experimental_data(:,1),steady_activation_experimental_data(:,2:end), 'r-' )
        ylim( [-1 2 ] )
        xlim( [ 0.5 6.5 ] )
        
    end
    subplot( 3,3,9)
    hold on
    plot( steady_activation_prediction_data(:,1),steady_activation_prediction_data(:,2:end),plot_strings{prot} )
    hold off
    if prot == 2
        ylim( [-1 2 ] )   
        xlim( [ 0.5 6.5 ] )
        xlabel( 'Time (s)' )
    end
    D=[];
    I=[];
    
    %% steady state activation protocol
    
    % Import protocol
    cd ../Protocols
    V=importdata(['steady_activation_protocol.mat']);
    cd ..
    cd Code
    % Identify model type for simulation
    %[~,model_type] =modeldata(model);
    
    % Import experimental data
    cd ../ExperimentalData
    cd(exp_ref)
    E=importdata(['steady_activation_',exp_ref,'_dofetilide_subtracted_leak_subtracted.mat']);
    cd ..
    cd ..
    cd Code
    % Simulate cell-specific and literature models
    I = SimulatingData(35,{'steady_activation'},P,V,temperature);
    
    V=[-60,-40,-20,0,20,40,60];
    
    for i=1:7
        
        D=I(56292+(82580*(i-1)):57292+82580*(i-1));
        
        S(i) = max(abs(D));
        if min(D) == -S(i)
            
            S(i) = -S(i);
        end
        
        D=[];
    end
    
    % Experimental peak currents identified manually
    if strcmp(exp_ref,'16713110')==1
        
        SExp = [0.0201,0.0357,0.2412,1.21,1.6253,1.6373,1.6473];
        
    end
    if strcmp(exp_ref,'16713003')==1
        
        SExp = [0.0286,0.0344,0.1205,0.6948,1.1962,1.3079,1.3802];
        
    end
    if strcmp(exp_ref,'16715049')==1
        
        SExp = [0.0363,0.0445,0.1645,0.5908,0.7596,0.7946,0.813];
        
    end
    
    if strcmp(exp_ref,'16704007')==1
        
        SExp = [0.0235,0.0764,0.6929,1.793,2.2334,2.2947,2.3034];
        
    end
    if strcmp(exp_ref,'16704047')==1
        
        SExp = [0.0792,0.1058,0.2722,0.7578,1.1566,1.382,1.4495];
        
    end
    
    if strcmp(exp_ref,'16708016')==1
        
        SExp = [0.05,0.0587,0.3118,0.6876,0.7373,0.7395,0.7423];
        
    end
    
    if strcmp(exp_ref,'16707014')==1
        
        SExp = [0.0539,0.0586,0.1294,0.3797,0.5509,0.5684,0.571];
        
    end
    
    if strcmp(exp_ref,'16708118')==1
        
        SExp = [0.0356,0.06,0.1799,0.3834,0.4283,0.4343,0.4459];
        
    end
    
    if strcmp(exp_ref,'16708060')==1
        
        SExp = [0.0387,0.0555,0.3067,0.6727,0.7875,0.7973,0.8142];
        
    end
    
    % Prepare and save data
    steady_activation_iv_experimental_data(:,1)=V;
    steady_activation_iv_experimental_data(:,2)=SExp./max(SExp);
    steady_activation_iv_prediction_data(:,1)=V;
    steady_activation_iv_prediction_data(:,2) =S./max(S);
    
    figure( fig_summary )
    if prot == 1
        plot( steady_activation_iv_experimental_data(:,1), steady_activation_iv_experimental_data(:,2), 'r--' )
    end
    hold on
    plot( steady_activation_iv_prediction_data(:,1), steady_activation_iv_prediction_data(:,2), plot_strings{ prot } )
    hold off
    xlabel( 'Voltage (mV)' )
    xlim( [ -60 60 ] )
    ylabel( 'Normalized Current (-)' )
    ylim( [ 0 1 ] )
        
end
