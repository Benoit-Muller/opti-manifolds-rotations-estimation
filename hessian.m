function h = hessian(data,X,S)
    % X: size (d,d,m)
    % S: size(d,d,m), is the skew part Omega of a tangent vector
    % H: size (d,d,M) with M=|E|
    % I,J: size (M)
    % with pdf notation: H(:,:,k) = H_(I(k),J(k))
    
    [d,~,m] = size(X);
    [~,~,M] = size(data.H);
    ma = length(data.A);
    Z = pagemtimes(X(:,:,data.I),'transpose', data.H,'none');
    Z = pagemtimes(Z, X(:,:,data.J)); % (d,d,M)
    [P1,P2] = frac(data.k1,data.k2,data.c1,data.c2,data.q,Z);  % (M)
    h = zeros(d,d,m-ma);
    for k=1:M
        if data.J(k) > ma
            Omega = Z(:,:,k)'*S(:,:,data.I(k)) - S(:,:,data.J(k))*Z(:,:,k)';
            h(:,:,data.J(k)-ma) = h(:,:,data.J(k)-ma) ...
                    + P1(k) * S(:,:,data.J(k)) * multiskew(Z(:,:,k)')...
                    + (P2(k)-P1(k)^2) * trace(Omega) * multiskew(Z(:,:,k)') ...
                    + P1(k) * multiskew(Omega);
            if data.I(k) > ma
                Omega = Z(:,:,k)*S(:,:,data.J(k)) - S(:,:,data.I(k))*Z(:,:,k);
                h(:,:,data.I(k)-ma) = h(:,:,data.I(k)-ma) ...
                    + P1(k) * S(:,:,data.I(k)) * multiskew(Z(:,:,k))...
                    + (P2(k)-P1(k)^2) * trace(Omega) * multiskew(Z(:,:,k)) ...
                    + P1(k) * multiskew(Omega);
            end
        end
    end
    h = - multiskew(h);
end

function [P1,P2] = frac(k1,k2,c1,c2,q,Z)
    etz = exp(multitrace(Z));
    p = q/c1 * etz.^k1 + (1-q)/c2 * etz.^k2;
    P1 = (q*k1/c1 * etz.^k1 + (1-q)*k2/c2 * etz.^k2) ./ p;
    P2 = (q*k1^2/c1 * etz.^k1 + (1-q)*k2^2/c2 * etz.^k2) ./ p;
end