function A = lowrankpluspolydecay(n,r,p)
    %
    % n - size of matrix
    % r - Rank of "significant part"
    % q - parameter that controls rate of decay
    
    
    um = randn(n,n);    [u,~] = qr(um);
    vm = randn(n,n);    [v,~] = qr(vm);

    % Generate the matrix A
    Sigma(1:r) = 1;
    Sigma((r+1):n) = (2:(n-r+1)).^(-p);
    A = u*diag(Sigma)*v';
end