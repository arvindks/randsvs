function A = decayingeigenvalues(n,gamma)
    %
    % Input
    % n - size of matrix
    % gamma   - eigenvalue gap
    
    % Output
    % A   - resulting matrix
    
    %Geometric distribution of eigenvalues

    kappa = (1/gamma)^(n-1);
    A = gallery('randsvd',n,kappa,5);


