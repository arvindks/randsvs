function [lb, ub] = sv_bounds(V,Omega,s,k,q)
    %
    % Computes lower and upper bounds for the singular values
    % 
    % Inputs
    % V     - right singular vectors
    % Omega - random starting guess
    % s     - singular values
    % k     - target rank
    % q     - number of subspace iterations
    
    % Outputs
    % lb    - lower bound for singular values
    % ub    - upper bound for singular values
    
    
    Omega  = V'*Omega;
    Omega1 = Omega(1:k,:);  Omega2 = Omega(k+1:end,:);
    Onrm = norm(Omega2/Omega1);
    
    s1 = s(1:k);    s2 = s(k+1:end);    
    ub = s1;                        %Upper bound
    
    tau = (s2(1)./s1).^(4*q+2);
    lb = s1./sqrt(1+tau*Onrm.^2);   %Lower bound  
    
end