function sintheta = subspace_angles(u,uh)
    % 
    % Computes the canonical angles between two subspaces U and Uh
    % Input
    % u     :   n x k matrix with orthonormal columns (basis for U) 
    % uh    :   n x k matrix with orthonormal columns (basis for Uh)
    % 
    % Output
    % sintheta :  k x 1 sin of the canonical angles 
    %             (sin theta_1, ..., sin  theta_k)
    
    A = u'*uh;
    s = svd(A);
    sintheta = sqrt(1-s.^2);
end