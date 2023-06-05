function M = manyrotationsfactory(d,m,J,Ra)
% anchors must be in the first indices, i.e. J = 1:m_J
    m_J = length(J);
    M = rotationsfactory(d, m - m_J);
    M.retr = M.retr_polar;
    M.J = J;
    M.Ra = Ra;
    M.add_anchors=@add_anchors;
    function X = add_anchors(X)
        X= cat(3,M.Ra,X);
    end
end