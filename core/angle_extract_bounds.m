function st_u = angle_extract_bounds(V,Omega,s,k,q)
    
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
    
    s1 = s(1:k);    s2 = s(k+1:end);    %partition singular value
    gamma = s2(1)/s1(end);
    
    tO  = (s2(1)./s1).*(gamma).^(2*q)*Onrm;
    st_u = tO./(1-gamma);

    
end