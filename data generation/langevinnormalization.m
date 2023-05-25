function c = langevinnormalization(n, kappa, scaling)
% Normalization coefficients for the isotropic Langevin distribution on
% SO(n), with concentration parameters kappa >= 0.
%
% Returns a vector c of the same length as kappa such that for all k,
% 
%  c(k) = c_n(kappa(k)),
%
% where the function c_n is defined in the paper
%
% In particular, c_n(0) = 1 for all n.
%
% If scaling is set to 1, the outputs are scaled down by exp(-n*kappa).
% By default, scaling is set to 0, which drastically limits the range of
% acceptable values of kappa, for numerical reasons. In most places, it is
% (way) preferable to activate scaling, as done in infoweight for example.
%
% This implementation is valid for n = 2, 3 or 4.

    assert(numel(n) == 1);
    
    if ~exist('scaling', 'var') || isempty(scaling)
        scaling = 0;
    end
    scaling = double(scaling);
    assert(scaling == 0 || scaling == 1);

    switch n
        
        case 2
            c = besseli(0, 2*kappa, scaling);
            
        case 3
            I0 = besseli(0, 2*kappa, scaling);
            I1 = besseli(1, 2*kappa, scaling);
            if scaling
                c = I0-I1;
            else
                c = exp(kappa).*(I0-I1);
            end 
            
        case 4
            I0 = besseli(0, 2*kappa, scaling);
            I1 = besseli(1, 2*kappa, scaling);
            I2 = besseli(2, 2*kappa, scaling);
            c = I0.^2 - 2*I1.^2 + I0.*I2;
            
        otherwise
            error('Invalid value for n. Should be 2, 3 or 4.');
        
    end


end
