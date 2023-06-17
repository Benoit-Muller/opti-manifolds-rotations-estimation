function M = manyrotationsfactory(d,m,J,Ra)
% anchors must be in the first indices, i.e. J = 1:m_a
    assert(all(J==(1:max(J))));
    m_a = length(J);
    M = rotationsfactory(d, m - m_a);
    M.retr = M.retr_polar;
    M.J = J;
    M.Ra = Ra;
    M.add_anchors = @add_anchors;
    function X = add_anchors(X)
        X= cat(3,M.Ra,X);
    end
    M.add_zeros = @add_zeros;
    function S = add_zeros(S)
        S= cat(3,zeros(size(M.Ra)),S);
    end
end