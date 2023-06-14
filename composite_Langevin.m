function p = composite_Langevin(k1,k2,c1,c2,q,Z)
    % kappa: float
    % c: float, normalizing constant
    % Z: size (d,d,M)
    p = exp(multitrace(Z));
    p = q/c1 * p.^k1 + (1-q)/c2 * p.^k2;
end


function p = Langevin(k,c,Z)
    % k: float
    % c: float, normalizing constant
    % Z: size (d,d,M)
    p = exp(k * multitrace(Z)) / c;
end

