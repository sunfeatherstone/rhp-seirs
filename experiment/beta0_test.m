
clear; clc;

dt_hr       = 0.5;       
T_period_hr = 168; 
nStepPerPer = T_period_hr/dt_hr;
S           = load('rhythm/preprocessed_data/weibo_spline_pp_98.mat');
pp          = mkpp(S.breaks, S.coefs);
F_hour      = @(t) ppval(pp, mod(t/2,84));
tGrid       = (0:nStepPerPer-1)' * dt_hr;
F_vec       = F_hour(tGrid);

beta0_seed = 8.1607e-01;
sigma0     = 7.6450e-01;
kappa      = 8.8634e+00;
delta      = 3.5503e-05;
gamma      = 7.2033e-02;

floquet_rho = @(b0) local_rho(b0, sigma0, gamma, delta, kappa, F_vec, dt_hr);

beta_min = beta0_seed * 0.05;   
beta_max = beta0_seed * 0.1005;  
nBeta    = 300;           
beta_vec = logspace(log10(beta_min), log10(beta_max), nBeta);

rho_vec = arrayfun(floquet_rho, beta_vec);
idx = find(rho_vec < 1 , 1, 'last');  
idx2 = idx + 1;                     
if isempty(idx) || rho_vec(idx2) < 1
    error('The beta_vec interval is not enough to cross rho=1, so expand beta_max');
end

beta_est = interp1(rho_vec([idx idx2]), beta_vec([idx idx2]), 1, 'linear');
fun      = @(b) floquet_rho(b) - 1; 
opts     = optimset('TolX',1e-12,'Display','off');
beta_star= fzero(fun, beta_vec([idx idx2]), opts);

fprintf('\nThreshold beta0*: %.10g (linear interpolation=%.4g)\n', beta_star, beta_est);


figure; semilogx(beta_vec, rho_vec,'k.-'); hold on;
plot(beta_star, 1,'ro','MarkerFaceColor','r');
xline(beta_star,'r--','\beta_0^*','LabelOrientation','horizontal',...
      'LabelVerticalAlignment','bottom');
yline(1,'--k'); grid on;
xlabel('\beta_0'); ylabel('\rho(\beta_0)');
title('Floquet multiplier vs. \beta_0');


function rho = local_rho(b0, s0, gamma, delta, kappa, Fvec, dt)
    Fstar = (Fvec.^kappa) ./ mean(Fvec.^kappa);
    Y = eye(2);
    a1 = 0.5 - sqrt(3)/6;
    a2 = 0.5 + sqrt(3)/6;
    n  = numel(Fstar);
    for j = 1:n-1
        B1 = b0 * Fstar(j);   S1 = s0 * Fstar(j);
        B2 = b0 * Fstar(j+1); S2 = s0 * Fstar(j+1);
        J1 = [-(S1+delta), B1; S1, -gamma];
        J2 = [-(S2+delta), B2; S2, -gamma];
        Y  = expm(dt*(a1*J1 + a2*J2)) * expm(dt*(a2*J1 + a1*J2)) * Y;
    end
    rho = max(abs(eig(Y)));
end