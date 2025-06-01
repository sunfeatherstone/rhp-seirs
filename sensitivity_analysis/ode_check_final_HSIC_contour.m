rng(20250530);  


clc; close all;
format long g

S  = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp = mkpp(S.breaks, S.coefs);
F_base = @(tau) ppval(pp, mod(tau, 84));
F_hour = @(t_hr) F_base(t_hr/2);

paramFiles = { ...
    'best_parameter/EventA.mat', ...
    'best_parameter/EventB.mat', ...
    'best_parameter/EventC.mat' };

nSample = 10000;     
D_feat  = 200;       

for f = 1:numel(paramFiles)
    data = load(paramFiles{f}, 'theta_refined');
    [ tbl, raw ]  = run_target_hsic(data.theta_refined, F_hour, nSample, D_feat);
    
    fprintf('\n=== %s : Target-oriented HSIC (ρ>1) ===\n', ...
            erase(paramFiles{f}, 'best_parameter/'));
    for r = 1:size(tbl,1)
        fprintf('%8s :  %.4f\n', tbl{r,1}, tbl{r,2});
    end
    score_stack(f,:) = raw; 

    theta_hat=data.theta_refined;

%========== (r_beta, r_gamma) ρ≈1 ==========
beta_hat  = theta_hat(1);
gamma_hat = theta_hat(6);

rB = logspace(log10(0.25), log10(4), 101);
rG = logspace(log10(0.25), log10(4), 101);
[RB, RG] = meshgrid(rB, rG);
rho_grid = zeros(size(RB));

for ii = 1:numel(RB)
    ii
    th_tmp      = theta_hat;
    th_tmp(1)   = log(RB(ii)) + (beta_hat);
    th_tmp(6)   = log(RG(ii)) + (gamma_hat);
    rho_grid(ii)= floquet_multiplier(th_tmp, F_hour);
end

fig2 = figure('Visible','off',...
              'Units','centimeters','Position',[2 2 14 12]);
ax  = axes(fig2);

Lgrid = log10(rho_grid);  

contourf(ax, RB, RG, Lgrid, 24, 'LineColor','none'); hold(ax,'on');
contour (ax, RB, RG, rho_grid, [1 1], 'k', 'LineWidth',2.0); 

plot    (ax, 1, 1, 'r^', 'MarkerFaceColor','r','MarkerSize',8);

ticks = [0.25 ,0.5,1,2, 4];          % 0.5, 0.6, 0.7, …, 2.0
set(ax, 'XScale','log', 'YScale','log', ...
        'XLim',[0.25 4], 'YLim',[0.25 4], ...
        'XTick',ticks, 'YTick',ticks, ...
        'FontSize',9,'Box','on', ...
        'GridAlpha',0.25,'LineWidth',0.7);
set(gca,'XMinorTick','off','YMinorTick','off','ZMinorTick','off')

% set(ax,'XTickLabel',arrayfun(@(v)sprintf('%.2f',v),ticks,'uni',0));
% set(ax,'YTickLabel',arrayfun(@(v)sprintf('%.2f',v),ticks,'uni',0));

grid(ax,'on');

% xtickangle(ax,45);   % ytickangle 

xlabel(ax,'$r_{\beta}=\beta_0/\hat{\beta}_0$',...
           'Interpreter','latex','FontSize',11);
ylabel(ax,'$r_{\gamma}=\gamma/\hat{\gamma}$',...
           'Interpreter','latex','FontSize',11);

[~, city] = fileparts(paramFiles{f});
title(ax, sprintf('%s: $\\rho(r_{\\beta},r_{\\gamma})$', city),...
          'Interpreter','latex');

cb = colorbar(ax);                 
cb.Title.String      = 'lg\rho';         
cb.Title.Interpreter = 'tex';            
cb.Ruler.TickLabelFormat = '%.1f'; 
% cb.Label.FontSize    = 9;
% cb.Label.Rotation    = 0;
% cb.Label.VerticalAlignment   = 'bottom';

pdfName = sprintf('rho_contour_%s.pdf', city);
exportgraphics(fig2, pdfName,...
               'ContentType','vector','BackgroundColor','none');
close(fig2);
end


paramNames = {'\beta_0','\sigma_0','\kappa','\gamma','\delta'};
cityNames  = {'EventA','EventB','EventC'};

fig = figure('Units','centimeters','Position',[2 2 24 10]);

bh = bar(paramNames, score_stack.', 'grouped');   

set(bh(1),'FaceColor',[0.04 0.28 0.62],'EdgeColor','none');   
set(bh(2),'FaceColor',[0.10 0.45 0.80],'EdgeColor','none');   
set(bh(3),'FaceColor',[0.55 0.75 0.92],'EdgeColor','none');   

ylabel('Target-oriented HSIC score','FontSize',11)
legend(cityNames,'Location','northeast','FontSize',9)
ylim([0 0.42]); grid on; box on;

exportgraphics(fig,'HSIC_grouped_bar.pdf', ...
               'ContentType','vector','BackgroundColor','none');








function [HSIC_table, raw_scores] = run_target_hsic(theta_refined, F_hour, nSample, D_feat)
    rng(0);                                
    names = {'beta0','sigma0','kappa','gamma','delta'};
    idx_use = [1 2 5 6 7];               
    theta_hat = theta_refined(:);
    p_hat  = theta_hat(idx_use);
    k      = numel(idx_use);

    % --- r = exp( log0.5 + U*(log2-log0.5) )
    S_unit = lhsdesign(nSample,k,'criterion','maximin','iterations',50);
    log_lo = log(0.25);   log_hi = log(4);
    Theta_s = p_hat.' + (log_lo + S_unit*(log_hi-log_lo));  % n×k

    % ---------- rho(Θ) ----------
    rho_vec = zeros(nSample,1);
    for ii = 1:nSample
        ii
        th_tmp          = theta_hat;
        th_tmp(idx_use) = Theta_s(ii,:);
        rho_vec(ii)     = floquet_multiplier(th_tmp, F_hour);
    end
    Z      = double(rho_vec > 1);          
    Zc     = Z - mean(Z);                
    Zc2    = sum(Zc.^2);
    
    S_HSIC = zeros(1,k);
    for j = 1:k
        Xj = Theta_s(:,j);
        sig = median(pdist(Xj));
        W   = randn(1,D_feat) / sig;
        b   = 2*pi*rand(1,D_feat);
        Phi = sqrt(2/D_feat) * cos( Xj*W + b );
        Phi = Phi - mean(Phi,1);
        c   = Phi' * Zc;
        hsicXZ = sum(c.^2) / nSample^2;
        M     = Phi' * Phi;
        hsicXX= sum(M(:).^2) / nSample^2;
        hsicZZ= (Zc2^2) / nSample^2;
        S_HSIC(j) = hsicXZ / sqrt(hsicXX*hsicZZ + eps);
    end
    
    raw_scores = S_HSIC; 
    [score_sorted, ord] = sort(S_HSIC,'descend');
    HSIC_table = [ names(ord).' , num2cell(score_sorted(:)) ];
end












%====================== floquet_multiplier.m ================
function rho = floquet_multiplier(theta20,F_hour)
    % -----  unpack ( beta0 sigma0 kappa gamma delta)
    beta0 = exp(theta20(1));
    sigma0= exp(theta20(2));
    kappa = exp(theta20(5));
    gamma = exp(theta20(6));
    delta = exp(theta20(7));

    T_tot = 168;   dt = 0.5;  nStep = T_tot/dt;
    t_hr  = (0:nStep-1).' * dt;
    F     = F_hour(t_hr);          
    Fp    = F.^kappa;  F_star = Fp / mean(Fp);

    y = eye(2);
    for kk = 1:nStep
        b = beta0 * F_star(kk);
        s = sigma0* F_star(kk);
        J = [- (s+delta) ,  b ;
               s         , -gamma];
        y = y + dt * J * y;         
    end
    rho = max(abs(eig(y)));
end


function score = hsic_rff_1d(X, Z, D)
n  = numel(X);
Zc = Z - mean(Z);
sig = median(pdist(X));
W   = randn(1,D) / sig;
b   = 2*pi*rand(1,D);
Phi = sqrt(2/D) * cos( X*W + b );
Phi = Phi - mean(Phi,1);

c      = Phi' * Zc;                 
hsicXZ = sum(c.^2)       / n^2;
M      = Phi' * Phi;              
hsicXX = sum(M(:).^2)    / n^2;
hsicZZ = (Zc'*Zc)^2      / n^2;

score  = hsicXZ / sqrt(hsicXX*hsicZZ + eps);
end

