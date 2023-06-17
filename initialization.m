function X = initialization(data)
% compute a first initializaition guess X of size (d,d,m-ma)
    m = data.m;
    d = data.d;
    I=data.I;
    J=data.J;
    D = spdiags(data.deg, 0, m, m);
    D1 = kron(D,speye(d));
    W = zeros(d*m, d*m);
    period = 1:d;
    for k=1:data.cardE
        W((I(k)-1)*d + period,(J(k)-1)*d + period) = data.H(:,:,k);
    end
    W = W + W';
    [Y,S] = eigs(W, D1, d);
    Y = Y*sqrt(trace(D));
    tol = 1e-6; % tolerance for assertions
    assert(norm(W*Y-D1*Y*S) < tol) % GEP
    assert(norm(Y'*D1*Y - trace(D)*eye(d)) < tol) % feasibility
    Xa = zeros(d,d,m);
    Xb = zeros(d,d,m);
    for i=1:m
        Xa(:,:,i) = proj(Y((i-1)*d + period,:));
        Xb(:,:,i) = Y((i-1)*d + period,:);
        Xb(:,end,i) = - Xb(:,end,i);
        Xb(:,:,i) = proj(Xb(:,:,i)) ;
    end
    if cost(data,Xa) < cost(data,Xb)
        Xtilde=Xa;
    else
        Xtilde=Xb;
    end
    Q = pagemtimes(Xtilde(:,:,data.A),'transpose', data.Ra,'none');
    Q = proj(sum(Q,3));
    X = pagemtimes(Xtilde(:,:,(data.ma+1):end), Q);
end

function R = proj(A)
    % compute the projection of A to SO(d)
    tol = 1e-3;
    [U,~,V] = svd(A);
    if abs(det(U)-det(V)) > 0.1
        U(:,end) = - U(:,end);
    end
    R = U*V';
    assert(norm(R'*R - eye(size(A))) < tol);
    assert(abs(det(R) - 1) < tol);
end