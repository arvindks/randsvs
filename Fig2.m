close all
clc

% Plots the bounds for the canonical angles b/w subspaces
% Generates plots corresponding to Figure 2 in the paper
% 
% Saibaba, Arvind K. "Randomized subspace iteration: 
% Analysis of canonical angles and unitarily invariant norms." 
% SIAM Journal on Matrix Analysis and Applications 40.1 (2019): 23-48.
% 


%% Controlled gap
rng(100) 
m = 3000; n = 300;  % Matrix dimensions
kmax = 25;          % Target rank 
r    = 15;          % Location of the eigenvalue gap

% Random matrices
p = 20;             % Oversampling parameter
Omega_g = randn(n,kmax+p);      %Gaussian

gap = [1,2,10];

for g = 1:3
    
    A = controlledgap(m,n,r,gap(g));   
    figure,
    color = 'rbk';
    [u,s,v] = svds(A,kmax+1);
    s = diag(s);
    for q = 0:2
       
        % Subspace iteration
        Q = subspaceiteration(A,Omega_g,q);


        % Singular vector angle bounds
        [uq,sq,vq] = randsvd(A,Q,kmax+p);
        [stu,stv] = angle_bounds(v,Omega_g,s,kmax,q);
        sinthetau = subspace_angles(u(:,1:kmax),uq);   

        semilogy(1:kmax,real(sinthetau),strcat('-',color(q+1)),...
                1:kmax,stu + 1.e-16,strcat('--',color(q+1)),'LineWidth',3)
        hold on

    end
    
    
    if g == 1
        legend({'Computed q = 0', 'Estimate q = 0', ...
            'Computed q = 1', 'Estimate q = 1', ...
            'Computed q = 2', 'Estimate q = 2'},'Location','southeast');
        ylabel('$\sin\theta_j$','Interpreter','LaTeX','FontSize',24)
    end
    xlabel('Index')
    set(gca,'FontSize',20)
    
    

    
    
    fname = strcat('figs/canon_gap_',num2str(g));
    print(fname, '-depsc') 
end

%% Low rank plus noise
rng(100);

n = 300;            % Matrix dimensions
r    = 15;          % Location of the eigenvalue gap

noise = [1.e-2,1.e-1,1];
color = 'rbk';
for g = 1:3
    
    A = lowrankplusnoise(n,r,noise(g));   
    [u,s,v] = svd(A);
    uk = u(:,1:kmax);
    s = diag(s);    s = s(1:(kmax+1));
    figure(g+3),
    
    for q = 0:2
        
        % Subspace iteration
        Q = subspaceiteration(A,Omega_g,q);

        % Singular vector angle bounds
        [uq,sq,vq] = randsvd(A,Q,kmax+p);
        [stu,stv] = angle_bounds(v,Omega_g,s,kmax,q);
        sinthetau = subspace_angles(uk,uq);   

        semilogy(1:kmax,real(sinthetau),strcat('-',color(q+1)),...
                1:kmax,stu + 1.e-16,strcat('--',color(q+1)),'LineWidth',3)
        hold on

    end

    if g == 1
        ylabel('$\sin\theta_j$','Interpreter','LaTeX','FontSize',24)
    end
    
    xlabel('Index')
    set(gca,'FontSize',20)
    
    
    fname = strcat('figs/canon_noise_',num2str(g));
    print(fname, '-depsc') 
end


%% Low rank plus poly decay
rng(100);

decay = [0.5,1,2];
color = 'rbk';
for g = 1:length(decay)
    
    A = lowrankpluspolydecay(n,r,decay(g));   
    [u,s,v] = svd(A);
    s = diag(s); s = s(1:(kmax+1));
    
    
    figure(g+6),
    sinthetau = zeros(kmax,3);
    for q = 0:2
        
        % Subspace iteration
        Q = subspaceiteration(A,Omega_g,q);

        % Singular vector angle bounds
        [uq,~,~] = randsvd(A,Q,kmax+p);
        [stu,stv] = angle_bounds(v,Omega_g,s,kmax,q);
        sinthetau(:,q+1) = real(subspace_angles(u(:,1:kmax),uq)); 
        
        semilogy(1:kmax,sinthetau(:,q+1) + 1.e-16,strcat('-',color(q+1)),...
                1:kmax,stu + 1.e-16,strcat('--',color(q+1)),'LineWidth',3)
        hold on
        
    end
    
    mins = min(min(sinthetau)) + 1.e-16;
    maxs = max(max(sinthetau));
    
    axis([1,kmax,0.01*mins,10*maxs])
 
    
    if g == 1
        ylabel('$\sin\theta_j$','Interpreter','LaTeX','FontSize',24)
    end

    xlabel('Index')
    set(gca,'FontSize',20)
    
    fname = strcat('figs/canon_decay_',num2str(g));
    print(fname, '-depsc') 
end

    
