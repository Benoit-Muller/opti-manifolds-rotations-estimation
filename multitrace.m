function t = multitrace(Z)
    % Z: size (d,d,M)
    % output t of size (1,1,M) with page-wise trace: t(1,1,k) = trace(Z(:,:,k))
    [d,~,~] = size(Z);
    t = zeros(1,1,d);
    for i=1:d
        t = t + Z(i,i,:);
    end
    t = squeeze(t);
end