function [ mx, factors ] = BuildDataMatrix( conf, variables, I, rates, fluxes, ...
                                            sensitivity, dIdp, voltage, dVdt, time, input_factors )
%BUILDDATAMATRIX Summary of this function goes here
%   Detailed explanation goes here

mx = [];
normalise = any(strcmp( conf, 'normalise' ) );

if normalise
    if nargin == 11 % useful when making predictions from a model built on different data
        % unpack the factors vector
        fac_var  = input_factors( 1 );
        fac_I    = input_factors( 2 );
        fac_rate = input_factors( 3 );
        fac_flux = input_factors( 4 );
        fac_sens = input_factors( 5 );
        fac_dIdp = input_factors( 6 );
        fac_V    = input_factors( 7 );
        fac_dVdt = input_factors( 8 );
        fac_t    = input_factors( 9 );
    else % Generate the factors vector
        fac_var  = max(max( abs( variables ) ) );
        fac_I    = max(max( abs( I ) ) );
        fac_rate = max(max( abs( rates ) ) );
        fac_flux = max(max( abs( fluxes ) ) );
        fac_sens = max(max( abs( sensitivity ) ) );
        fac_dIdp = max(max( abs( dIdp ) ) );
        fac_V    = max(max( abs( voltage ) ) );
        fac_dVdt = max(max( abs( dVdt ) ) );
        fac_t    = max(max( abs( time ) ) );
    end
else
    fac_var  = 1;
    fac_I    = 1;
    fac_rate = 1;
    fac_flux = 1;
    fac_sens = 1;
    fac_dIdp = 1;
    fac_V    = 1;
    fac_dVdt = 1;
    fac_t    = 1;
end

factors = [ fac_var fac_I, fac_rate, fac_flux, fac_sens, fac_dIdp, fac_V, fac_dVdt, fac_t ];

if any(strcmp( conf, 'const' ) )
    mx = [ mx ones( length( I ), 1 ) ];
end
if any(strcmp( conf, 'variables' ) )
    mx = [ mx variables/fac_var ];
end
if any(strcmp( conf, 'I' ) )
    mx = [ mx I/fac_I ];
end
if any(strcmp( conf, 'rates' ) )
    mx = [ mx rates/fac_rate ];
end
if any(strcmp( conf, 'fluxes' ) )
    mx = [ mx fluxes/fac_flux ];
end
if any(strcmp( conf, 'sensitivity' ) )
    mx = [ mx sensitivity/fac_sens ];
end
if any(strcmp( conf, 'dIdp' ) )
    mx = [ mx dIdp/fac_dIdp ];
end
if any(strcmp( conf, 'voltage' ) )
    mx = [ mx voltage/fac_V ];
end
if any(strcmp( conf, 'dVdt' ) )
    mx = [ mx dVdt/fac_dVdt ];
end
if any(strcmp( conf, 'time' ) )
    mx = [ mx time/fac_t ];
end
    
end

