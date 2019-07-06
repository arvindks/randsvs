close all
clc


% Plots the singular values of the various matrices
% Generates plots corresponding to Figure 1 in the paper
% 
% Saibaba, Arvind K. "Randomized subspace iteration: 
% Analysis of canonical angles and unitarily invariant norms." 
% SIAM Journal on Matrix Analysis and Applications 40.1 (2019): 23-48.


rng(100)
kmax = 35;          % Target rank

%% Plot controlled gap
m = 3000; n = 300;  % Matrix dimensions

color = 'rgk';
gap = [1,2,10];
figure,
for g = 1:3
    A = controlledgap(m,n,r,gap(g));   
    s = svds(A,kmax+1);
    semilogy(1:kmax,s(1:kmax),color(g),'Linewidth',3.0); hold on
end
legend('gap = 1', 'gap = 2', 'gap = 10')
xlabel('Index')
ylabel('Singular Values')
set(gca,'FontSize',20)
print('figs/sv_gap','-depsc')

%% Plot low-rank-plus-noise singular values.
rng(100)
n = 300;            % Matrix dimensions
noise = [1.e-2,1.e-1,1];

figure,
for g = 1:length(decay)
    
    A = lowrankplusnoise(n,r,noise(g));   
    [u,s,v] = svd(A);
    s = diag(s);
    
    semilogy(1:kmax,s(1:kmax),color(g),'Linewidth',3.0), hold on
end

legend({'noise = $10^{-2}$', 'noise = $10^{-1}$', 'noise = 1.0'},...
    'Location','southwest','Interpreter','LaTeX')
xlabel('Index')
set(gca,'FontSize',20)
print('figs/sv_noise','-depsc')

%% Plot low-rank-plus-poly-decay singular values.
rng(100)
n = 300;            % Matrix dimensions

decay = [0.5,1,2];

figure,
for g = 1:length(decay)
    
    A = lowrankpluspolydecay(n,r,decay(g));   
    [u,s,v] = svd(A);
    s = diag(s);
    
    semilogy(1:kmax,s(1:kmax),color(g),'Linewidth',3.0), hold on
end

legend({'decay = 0.5', 'decay = 1.0', 'decay = 2.0'},...
    'Location','southwest')
xlabel('Index')
set(gca,'FontSize',20)
print('figs/sv_decay','-depsc')
