function M = manyrotationsfactory(d,m,J,R)
    m_J = length(J);
    N = rotationsfactory(d, m - m_J);
    N.retr = N.retr_polar;

    M.dim = N.dim();
    M.inner = N.inner;
    M.norm = N.norm(x,d)^2;
    M.proj = N.proj;
    M.retr = N.retr;
    M.rand = N.rand;
    M.randvec =  N.randvec;
    M.zerovec = N.zerovec;
    M.lincomb = M.lincomb();
    M.tangent = N.tangent();

    M.tangent2ambient(X, S)
end