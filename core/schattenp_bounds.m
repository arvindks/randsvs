function bound = schattenp_bounds(V,Omega,s,k,q,p)
    %
    % Computes lower and upper bounds for the singular values
    % 
    % Inputs
    % V     - right singular vectors
    % Omega - random starting guess
    % s     - singular values
    % k     - target rank
    % q     - number of subspace iterations
    % p     - parameter of schatten norm 
    % 
    % Outputs
    % lb    - lower bound for singular values
    % ub    - upper bound for singular values
    
    
    Omega  = V'*Omega;
    Omega1 = Omega(1:k,:);  Omega2 = Omega(k+1:end,:);
    
    s1 = s(1:k);    s2 = s(k+1:end);    
    gamma = s2(1)/s1(end); mat = diag(s2)*(Omega2/Omega1);
    SOnrm = schattenp(mat,p);
    bound =  sqrt(norm(s2,p)+ gamma.^(4*q)*SOnrm);
     
end