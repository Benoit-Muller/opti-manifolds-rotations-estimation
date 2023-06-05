function f = cost(data,X)
    % X: size (d,d,m)
    % H: size (d,d,M) with M=|E|
    % I,J: size (M)
    % with pdf notation: H(:,:,k) = H_(I(k),J(k))
    
    Z = pagemtimes(X(:,:,data.I),'transpose', data.H,'none');
    Z = pagemtimes(Z, X(:,:,data.J));
    f = -log(composite_Langevin(data.k1,data.k2,data.c1,data.c2,data.q,Z));
end