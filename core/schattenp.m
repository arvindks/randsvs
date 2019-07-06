function Anrm = schattenp(A,p)
    %
    % Computes the schatten-p norm of the matrix
    % 
    % Input:
    % A - matrix
    % p - integer
    % 
    % Output
    % Anrm - Schatten-p norm
    
    s = svd(A,0);
    Anrm = norm(s,p);
end