close all
clc


% Plots the bounds for the singular values
% Generates plots corresponding to Figure 4 in the paper
% 
% Saibaba, Arvind K. "Randomized subspace iteration: 
% Analysis of canonical angles and unitarily invariant norms." 
% SIAM Journal on Matrix Analysis and Applications 40.1 (2019): 23-48.



%% Controlled gap
m = 3000; n = 300;  % Matrix dimensions
kmax = 25;          % Target rank 
r    = 15;          % Location of the eigenvalue gap


% Random matrices
p = 20;
Omega_g = randn(n,kmax+p);      %Gaussian

gap = 1;  % small gap
A = controlledgap(m,n,r,gap);
[~,s,v] = svds(A,kmax+1);
s = diag(s); 
for q = 0:2
    
    % Subspace iteration
    Q = subspaceiteration(A,Omega_g,q);
    
    
    % Singular vector angle bounds
    [uq,sq,vq] = randsvd(A,Q,kmax+p);
    [lb,ub] = sv_bounds(v,Omega_g,s,kmax,q);
    
    figure, 
    semilogy(1:kmax,ub,'k--','LineWidth',7.0), hold on
    semilogy(1:kmax,sq(1:kmax),'r-','LineWidth',3.0)
    semilogy(1:kmax,lb,'b--','LineWidth',3.0);
    hold on
    xlabel('Index')
    if q == 0
        ylabel('Singular values', 'FontSize', 24)
        legend('Upper bound', 'Computed', 'Lower bound')
    end
    set(gca,'FontSize',20)
    fname = strcat('figs/svs_gap_q_',num2str(q));
    print(fname, '-depsc')

    
end


%% Low-rank plus noise example
rng(100)
n = 300;  % Matrix dimensions
noise = 1.e-2; %Medium noise

A = lowrankplusnoise(n,r,noise);

[~,s,v] = svd(A);
s = diag(s);s = s(1:(kmax+1));
for q = 0:2
    % Subspace iteration
    Q = subspaceiteration(A,Omega_g,q);
    
    
    % Singular vector angle bounds
    [uq,sq,vq] = randsvd(A,Q,kmax+p);
    [lb,ub] = sv_bounds(v,Omega_g,s,kmax,q);
    
    figure, 
    semilogy(1:kmax,ub,'k--','LineWidth',7.0), hold on
    semilogy(1:kmax,sq(1:kmax),'r-','LineWidth',3.0)
    semilogy(1:kmax,lb,'b--','LineWidth',3.0);
    hold on
    xlabel('Index')
    if q == 0
        ylabel('Singular values', 'FontSize', 24)
        legend('Upper bound', 'Computed', 'Lower bound')
    end
    set(gca,'FontSize',20)
    fname = strcat('figs/svs_noise_q_',num2str(q));
    print(fname, '-depsc')
    
end


%% Low-rank plus decay example
rng(100)
n = 300;  % Matrix dimensions
decay = 1.; %Medium decay

A = lowrankpluspolydecay(n,r,decay);

[u,s,v] = svd(A);
s = diag(s);s = s(1:(kmax+1));
for q = 0:2
    % Subspace iteration
    Q = subspaceiteration(A,Omega_g,q);
    
    
    % Singular vector angle bounds
    [uq,sq,vq] = randsvd(A,Q,kmax+p);
    [lb,ub] = sv_bounds(v,Omega_g,s,kmax,q);
    
    figure, 
    semilogy(1:kmax,ub,'k--','LineWidth',7.0), hold on
    semilogy(1:kmax,sq(1:kmax),'r-','LineWidth',3.0)
    semilogy(1:kmax,lb,'b--','LineWidth',3.0);
    hold on
    xlabel('Index')
    if q == 0
        ylabel('Singular values', 'FontSize', 24)
        legend({'Upper bound', 'Computed', 'Lower bound'},...
                    'Location','southwest')
    end
    set(gca,'FontSize',20)
    fname = strcat('figs/svs_decay_q_',num2str(q));
    print(fname, '-depsc')
    
end

     
