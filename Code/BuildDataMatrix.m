function [ mx, factors, names ] = BuildDataMatrix( conf, variables, I, rates, fluxes, ...
                                            sensitivity, dIdp, voltage, dVdt, time, input_factors )
%BUILDDATAMATRIX Summary of this function goes here
%   Detailed explanation goes here

mx = [];
normalise = any(strcmp( conf, 'normalise' ) );

% Note - others only have one entry
[~, var_ind_idx, ~]= getLinearIndependent(variables, 1);
[~, rate_ind_idx, ~]= getLinearIndependent(rates, 1);
[~, flux_ind_idx, ~]= getLinearIndependent(fluxes, 1);
[~, sens_ind_idx, ~]= getLinearIndependent(sensitivity, 1);
[~, dIdp_ind_idx, ~]= getLinearIndependent(dIdp, 1);

if normalise
    if nargin == 11 % useful when making predictions from a model built on different data
        % unpack the factors vector
        min_fac_var  = input_factors( 1, 1 : 4 );
        max_fac_var  = input_factors( 2, 1 : 4 );
        min_fac_I    = input_factors( 1, 5 );
        max_fac_I    = input_factors( 2, 5 );
        min_fac_rate = input_factors( 1, 6 : 13 );
        max_fac_rate = input_factors( 2, 6 : 13 );
        min_fac_flux = input_factors( 1, 14 : 21 );
        max_fac_flux = input_factors( 2, 14 : 21 );
        min_fac_sens = input_factors( 1, 22 : 45 );
        max_fac_sens = input_factors( 2, 22 : 45 );
        min_fac_dIdp = input_factors( 1, 46 : 53 );
        max_fac_dIdp = input_factors( 2, 46 : 53 );
        min_fac_V    = input_factors( 1, 54 );
        max_fac_V    = input_factors( 2, 54 );
        min_fac_dVdt = input_factors( 1, 55 );
        max_fac_dVdt = input_factors( 2, 55 );
        min_fac_t    = input_factors( 1, 56 );
        max_fac_t    = input_factors( 2, 56 );
    else % Generate the factors vector
        min_fac_var  = min( variables );
        max_fac_var  = max( variables );
        min_fac_I    = min( I );
        max_fac_I    = max( I );
        min_fac_rate = min( rates );
        max_fac_rate = max( rates );
        min_fac_flux = min( fluxes );
        max_fac_flux = max( fluxes );
        min_fac_sens = min( sensitivity );
        max_fac_sens = max( sensitivity );
        min_fac_dIdp = min( dIdp );
        max_fac_dIdp = max( dIdp );
        min_fac_V    = min( voltage );
        max_fac_V    = max( voltage );
        min_fac_dVdt = min( dVdt );
        max_fac_dVdt = max( dVdt );
        min_fac_t    = min( time );
        max_fac_t    = max( time );
    end
else
        min_fac_var  = zeros( 1, 4 );
        max_fac_var  = ones( 1, 4 );
        min_fac_I    = 0;
        max_fac_I    = 1;
        min_fac_rate = zeros( 1, 8 );
        max_fac_rate = ones( 1, 8 );
        min_fac_flux = zeros( 1, 8 );
        max_fac_flux = ones( 1, 8 );
        min_fac_sens = zeros( 1, 24 );
        max_fac_sens = ones( 1, 24 );
        min_fac_dIdp = 0;
        max_fac_dIdp = 1;
        min_fac_V    = 0;
        max_fac_V    = 1;
        min_fac_dVdt = 0;
        max_fac_dVdt = 1;
        min_fac_t    = 0;
        max_fac_t    = 1;
end

factors = [ min_fac_var min_fac_I, min_fac_rate, min_fac_flux, min_fac_sens, min_fac_dIdp, min_fac_V, min_fac_dVdt, min_fac_t;... 
            max_fac_var max_fac_I, max_fac_rate, max_fac_flux, max_fac_sens, max_fac_dIdp, max_fac_V, max_fac_dVdt, max_fac_t];

names = {};        
        
if any(strcmp( conf, 'const' ) )
    mx = [ mx ones( length( I ), 1 ) ];
    name = GetNames( 'const' );
    names = [ names, name ];
end
if any(strcmp( conf, 'variables' ) )    
    mx = [ mx ( variables( :, var_ind_idx ) - min_fac_var( var_ind_idx ) ) ./ (max_fac_var( var_ind_idx ) - min_fac_var( var_ind_idx ) ) ];
    name = GetNames( 'variables' );
    names = [ names, name( var_ind_idx ) ];
end
if any(strcmp( conf, 'I' ) )
    mx = [ mx ( I - min_fac_I ) ./ ( max_fac_I - min_fac_I ) ];
    name = GetNames( 'I' );
    names = [ names, name ];
end
if any(strcmp( conf, 'rates' ) )
    mx = [ mx ( rates( :, rate_ind_idx ) - min_fac_rate( rate_ind_idx ) ) ./ ( max_fac_rate( rate_ind_idx ) - min_fac_rate( rate_ind_idx ) ) ];
    name = GetNames( 'rates' );
    names = [ names, name( rate_ind_idx ) ];
end
if any(strcmp( conf, 'fluxes' ) )
    mx = [ mx ( fluxes( :, flux_ind_idx ) - min_fac_flux( flux_ind_idx ) ) ./ (max_fac_flux( flux_ind_idx ) - min_fac_flux( flux_ind_idx ) ) ];
    name = GetNames( 'fluxes' );
    names = [ names, name( flux_ind_idx ) ];
end
if any(strcmp( conf, 'sensitivity' ) )
    mx = [ mx ( sensitivity( :, sens_ind_idx ) - min_fac_sens( sens_ind_idx ) ) ./ ( max_fac_sens( sens_ind_idx ) - min_fac_sens( sens_ind_idx ) ) ];
    name = GetNames( 'sensitivity' );
    names = [ names, name( sens_ind_idx ) ];
end
if any(strcmp( conf, 'dIdp' ) )
    mx = [ mx ( dIdp( :, dIdp_ind_idx ) - min_fac_dIdp( dIdp_ind_idx ) ) ./ ( max_fac_dIdp( dIdp_ind_idx ) - min_fac_dIdp( dIdp_ind_idx ) ) ];
    name = GetNames( 'dIdp' );
    names = [ names, name( dIdp_ind_idx ) ];
end
if any(strcmp( conf, 'voltage' ) )
    mx = [ mx ( voltage - min_fac_V ) ./ ( max_fac_V - min_fac_V ) ];
    name = GetNames( 'voltage' );
    names = [ names, name ];
end
if any(strcmp( conf, 'dVdt' ) )
    mx = [ mx ( dVdt - min_fac_dVdt ) ./ ( max_fac_dVdt - min_fac_dVdt ) ];
    name = GetNames( 'dVdt' );
    names = [ names, name ];
end
if any(strcmp( conf, 'time' ) )
    mx = [ mx ( time - min_fac_t ) ./ ( max_fac_t - min_fac_t ) ];
    name = GetNames( 'time' );
    names = [ names, name ];
end
    
end

