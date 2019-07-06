close all
clc


% Plots the bounds for the extraction of subspaces
% Generates plots corresponding to Figure 3 in the paper
% 
% Saibaba, Arvind K. "Randomized subspace iteration: 
% Analysis of canonical angles and unitarily invariant norms." 
% SIAM Journal on Matrix Analysis and Applications 40.1 (2019): 23-48.



%% Controlled gap
rng(100) 
m = 3000; n = 300;  % Matrix dimensions
kmax = 15;          % Target rank 
r    = 15;          % Location of the eigenvalue gap


% Random matrices
p = 20;
Omega_g = randn(n,kmax+p);      %Gaussian
gap = [1,2,10];

for g = 1:3
    
    A = controlledgap(m,n,r,gap(g));  
    [u,s,v] = svds(A,kmax+1);
    s = diag(s);
    
    figure(g),
    color = 'rbk';
    for q = 0:2
        

        % Subspace iteration
        Q = subspaceiteration(A,Omega_g,q);


        % Singular vector angle bounds
        [uq,sq,vq] = randsvd(A,Q,kmax);
        stu = angle_extract_bounds(v,Omega_g,s,kmax,q);
        sinthetau = subspace_angles(u(:,1:kmax),uq);   
        sinthetav = subspace_angles(v(:,1:kmax),vq);
        maxtheta  = max(real(sinthetau),real(sinthetav));
        semilogy(1:kmax,maxtheta,strcat('-',color(q+1)),...
                1:kmax,stu,strcat('--',color(q+1)),'LineWidth',3)
        hold on

    end
    if g == 1
        legend({'Computed q = 0', 'Estimate q = 0', ...
            'Computed q = 1', 'Estimate q = 1', ...
            'Computed q = 2', 'Estimate q = 2'},'Location','southeast');
        ylabel('$\max\{\sin\theta_j^\prime,\sin\nu_j^\prime\}$',...
            'Interpreter','LaTeX', 'FontSize',24)
    
    end
    xlabel('Index')
    set(gca,'FontSize',20)
    
    fname = strcat('figs/extract_gap_',num2str(g));
    print(fname, '-depsc') 
end


%% Low rank plus noise
rng(100)
n = 300;  % Matrix dimensions
noise = [1.e-2,1.e-1,1]; 

color = 'rbk';
for g = 1:3
    
    A = lowrankplusnoise(n,r,noise(g));   
    [u,s,v] = svd(A);
    uk = u(:,1:kmax);
    s = diag(s); s = s(1:(kmax+1));
    figure(g+3),
    
    for q = 0:2
        
        % Subspace iteration
        Q = subspaceiteration(A,Omega_g,q);


        % Singular vector angle bounds
        [uq,sq,vq] = randsvd(A,Q,kmax);
        stu = angle_extract_bounds(v,Omega_g,s,kmax,q);
        sinthetau = subspace_angles(u(:,1:kmax),uq);   
        sinthetav = subspace_angles(v(:,1:kmax),vq);
        maxtheta  = max(real(sinthetau),real(sinthetav));
        semilogy(1:kmax,maxtheta,strcat('-',color(q+1)),...
                1:kmax,stu,strcat('--',color(q+1)),'LineWidth',3)
        hold on

    end
    if g == 1
        ylabel('$\max\{\sin\theta_j^\prime,\sin\nu_j^\prime\}$',...
            'Interpreter','LaTeX','FontSize',24)
    end
    xlabel('Index')
    
    set(gca,'FontSize',20)
    
    fname = strcat('figs/extract_noise_',num2str(g));
    print(fname, '-depsc')
end



%% Low rank plus decay
rng(100)
n = 300;  % Matrix dimensions
decay = [0.5,1,2];

color = 'rbk';
for g = 1:3
    
    A = lowrankpluspolydecay(n,r,decay(g));   
    [u,s,v] = svd(A);
    uk = u(:,1:kmax);
    s = diag(s);    s = s(1:(kmax+1));
    figure(g+6),
    
    for q = 0:2
        
        % Subspace iteration
        Q = subspaceiteration(A,Omega_g,q);


        % Singular vector angle bounds
        [uq,sq,vq] = randsvd(A,Q,kmax);
        stu = angle_extract_bounds(v,Omega_g,s,kmax,q);
        sinthetau = subspace_angles(u(:,1:kmax),uq);   
        sinthetav = subspace_angles(v(:,1:kmax),vq);
        maxtheta  = max(real(sinthetau),real(sinthetav));
        semilogy(1:kmax,maxtheta,strcat('-',color(q+1)),...
                1:kmax,stu,strcat('--',color(q+1)),'LineWidth',3)
        hold on

    end
    if g == 1
        ylabel('$\max\{\sin\theta_j^\prime,\sin\nu_j^\prime\}$',...
            'Interpreter','LaTeX', 'FontSize',24)
    end
    xlabel('Index')
    
    set(gca,'FontSize',20)
    
    fname = strcat('figs/extract_decay_',num2str(g));
    print(fname, '-depsc')
end

