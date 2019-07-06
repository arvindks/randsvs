function Q = subspaceiteration(A,Omega,q)
    % Applies a stable version of subspace iteration  
    % Mathematically equivalent to range(A^q*Omega) 
    % 
    % Input
    % A       :  m x n matrix   (matrix of interest)
    % Omega   :  n x ell matrix (starting guess, ell <= min(m,n))
    % q       :  =>0 number of subspace  iterations
    % 
    % Output  
    % Q      :   m x ell matrix with orthonormal columns
    % 
    
    Y = A*Omega;
    [Q,~] = qr(Y,0);
    for i = 1:q
        Y = (A'*Q);  [Q,~] = qr(Y,0);
        Y = A*Q;     [Q,~] = qr(Y,0);
    end
    
end