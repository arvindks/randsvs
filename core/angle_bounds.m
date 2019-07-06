function [st_u, st_v] = angle_bounds(V,Omega,s,k,q)
    
    % Inputs
    % V     - right singular vectors
    % Omega - random starting guess
    % s     - singular values
    % k     - target rank
    % q     - number of subspace iterations
    
    % Outputs
    % st_u  - bounds for sin Theta(U_1,Uh)
    % st_v  - bounds for sin Theta(V_1,Vh)
    

    Omega  = V'*Omega;
    Omega1 = Omega(1:k,:);  Omega2 = Omega(k+1:end,:);
    Onrm = norm(Omega2/Omega1);
    
    s1 = s(1:k);    s2 = s(k+1:end);    %partition singular values
    tau = (s2(1)./s1).^(2*q+1); %eigenvalue gap
    
    
    tO  = tau*Onrm;
    st_u = tO./sqrt(1+tO.^2);
    
    t1 = (s2(1)./s1).^(2*q+2)*Onrm;
    st_v = t1./sqrt(1+t1.^2);
    
end