function t = multitrace(Z)
    % Z: size (d,d,M)
    % output t of size (1,1,M) with page-wise trace: t(1,1,k) = trace(Z(:,:,k))
    [d,~,M] = size(Z);
    t = zeros(1,1,M);
    for i=1:d
        t = t + Z(i,i,:);
    end
    t = squeeze(t);
end