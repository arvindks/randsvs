function [u,s,v] = randsvd(A,Q,k)
    % 
    % Randomized SVD 
    % 
    % Input
    % A       : m x  n matrix of interest
    % Q       : n x ell matrix with orthonormal columns 
    % k       : target rank (<= min(m,n)
    %  
    % Output
    % u       : m x k approximate left singular vectors
    % s       : k x 1 approximate singular values
    % v       : n x k approximate right singular vectors
    
    
    B = Q'*A;
    [ub,s,v] = svd(B,0);
    u = Q*ub;       u = u(:,1:k);
    s = diag(s);    s = s(1:k);
    v = v(:,1:k);
end