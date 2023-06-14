function [f,g] = costgrad(data,X)
    % X: size (d,d,m)
    % H: size (d,d,M) with M=|E|
    % I,J: size (M)
    % with pdf notation: H(:,:,k) = H_(I(k),J(k))
    [d,~,m] = size(X);
    ma = length(data.A);
    Z = pagemtimes(X(:,:,data.I),'transpose', data.H,'none');
    Z = pagemtimes(Z, X(:,:,data.J));
    [p,P] = composite_Langevin_frac(data.k1,data.k2,data.c1,data.c2,data.q,Z);
    f = -sum(log(p),'all');
    somme = permute(P,[2,3,1]) .* Z;
    [~,~,M] = size(data.H);
    g = zeros(d,d,m-ma);
    for k=1:M
        if data.J(k) > ma
            g(:,:,data.J(k)-ma) = g(:,:,data.J(k)-ma) + somme(:,:,k)';
            if data.I(k) > ma
                g(:,:,data.I(k)-ma) = g(:,:,data.I(k)-ma) + somme(:,:,k);
            end
        end
    end
    g=multiskew(g);
    g=-g;
end

function [p,P] = composite_Langevin_frac(k1,k2,c1,c2,q,Z)
    etz = exp(multitrace(Z));
    p = q/c1 * etz.^k1 + (1-q)/c2 * etz.^k2;
    P = (q*k1/c1 * etz.^k1 + (1-q)*k2/c2 * etz.^k2) ./ p;
end