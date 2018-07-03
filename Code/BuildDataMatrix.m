function [ mx, factors ] = BuildDataMatrix( conf, variables, I, rates, fluxes, ...
                                            sensitivity, dIdp, voltage, dVdt, time, input_factors )
%BUILDDATAMATRIX Summary of this function goes here
%   Detailed explanation goes here

mx = [];
normalise = any(strcmp( conf, 'normalise' ) );

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

if any(strcmp( conf, 'const' ) )
    mx = [ mx ones( length( I ), 1 ) ];
end
if any(strcmp( conf, 'variables' ) )
    mx = [ mx ( variables - min_fac_var ) ./ (max_fac_var - min_fac_var ) ];
end
if any(strcmp( conf, 'I' ) )
    mx = [ mx ( I - min_fac_I ) ./ ( max_fac_I - min_fac_I ) ];
end
if any(strcmp( conf, 'rates' ) )
    mx = [ mx ( rates - min_fac_rate ) ./ ( max_fac_rate - min_fac_rate ) ];
end
if any(strcmp( conf, 'fluxes' ) )
    mx = [ mx ( fluxes - min_fac_flux ) ./ (max_fac_flux - min_fac_flux ) ];
end
if any(strcmp( conf, 'sensitivity' ) )
    mx = [ mx ( sensitivity - min_fac_sens ) ./ ( max_fac_sens - min_fac_sens ) ];
end
if any(strcmp( conf, 'dIdp' ) )
    mx = [ mx ( dIdp - min_fac_dIdp ) ./ ( max_fac_dIdp - min_fac_dIdp ) ];
end
if any(strcmp( conf, 'voltage' ) )
    mx = [ mx ( voltage - min_fac_V ) ./ ( max_fac_V - min_fac_V ) ];
end
if any(strcmp( conf, 'dVdt' ) )
    mx = [ mx ( dVdt - min_fac_dVdt ) ./ ( max_fac_dVdt - min_fac_dVdt ) ];
end
if any(strcmp( conf, 'time' ) )
    mx = [ mx ( time - min_fac_t ) ./ ( max_fac_t - min_fac_t ) ];
end
    
end

