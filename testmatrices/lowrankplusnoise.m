function A = lowrankplusnoise(n,r,gamma)
    %
    % Input
    % m,n   - size of matrix
    % r     - significant part of spectrum
    % gamma - SNR
    % 
    % Output
    % A   - resulting matrix


    I = [eye(r),zeros(r,n-r); zeros(n-r,n)];
    fact = sqrt(gamma*r./(2*n^2));
    
    G = randn(n,n);
    A = I + fact*(G+G');

end