load P % Best result from Kylie's MCMC
IC = [0.0,0.0,0.0];
vv=importdata('sine_wave_protocol.mat');
I_exp = importdata('sine_wave_16713110_dofetilide_subtracted_leak_subtracted.mat');
ProtocolLength = length( vv )/10;
params = [ 1, P ];
X=MexHHSens([0:0.1:ProtocolLength],IC,params);
temperature = 21.4;
T = 273.15+temperature;
F = 96485;
R = 8314;
K_i = 130;
k_o = 4;
erev = ((R*T)/F)*log(k_o/K_i);
O = ones(length(X),1);
Vr = O.*erev;

I = params(length(params)).*X.*(vv-Vr);
O=X;
figure; plot( 0.1:0.1:ProtocolLength, I, 0.1 :0.1:ProtocolLength, I_exp );
figure; plot( 0.1:0.1:ProtocolLength, I-I_exp );