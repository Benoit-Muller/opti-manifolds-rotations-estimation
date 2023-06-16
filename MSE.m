function mse = MSE(data,X)
    % X of size (d,d,m-ma)
    m = data.m;
    ma = data.ma;
    mse = 0;
    for i=1:m-ma
        arg = data.R(:,:,ma-1+i)' * X(:,:,i);
        mse = mse + norm(logm(arg),"fro")^2;
    end
    mse = mse/(m-ma);
end