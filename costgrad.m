function [f,g] = costgrad(data,X)
    % X: size (d,d,m)
    % H: size (d,d,M) with M=|E|
    % I,J: size (M)
    % with pdf notation: H(:,:,k) = H_(I(k),J(k))
    Z = pagemtimes(X(data.I(data.A)),'transpose', data.H,'none');
    Z = pagemtimes(Z, X(data.J(data.A)));
    [p,P] = composite_Langevin_frac(data.k1,data.k2,data.c1,data.c2,data.q,Z);
    f = - sum(log(p),'all');
    g = - P(:)' * Z(:);
end

function [p,P] = composite_Langevin_frac(k1,k2,c1,c2,q,Z)
    etz = exp(multitrace(Z));
    p = (1-q)/c1*etz^k1 + q/c2*etz^k2;
    P = ((1-q)*k1/c1*etz^k1 + q*k2/c2*etz^k2) ./ p;
end