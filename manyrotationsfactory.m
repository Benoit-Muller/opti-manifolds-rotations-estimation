function M = manyrotationsfactory(d,m,J,Ra)
    m_J = length(J);
    M = rotationsfactory(d, m - m_J);
    M.retr = N.retr_polar;
    M.J = J;
    M.Ra = Ra;
    M.add_anchors=@add_anchors;
    function X = add_anchors(X)
        cat(3,Ra,X)
    end
end