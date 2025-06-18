clear; clc;


dt_hr   = 1;
Tper_hr = 168;
nPer    = Tper_hr/dt_hr;        

Spp     = load('rhythm/preprocessed_data/weibo_spline_pp_98.mat');
pp      = mkpp(Spp.breaks, Spp.coefs);
F_hour  = @(tHr) ppval(pp, mod(tHr/2,84));     
F_vec   = F_hour((0:nPer-1)' * dt_hr); 


POOL_FILE = 'experiment/param_pool.csv';
optsP     = detectImportOptions(POOL_FILE,'Delimiter',',');
P         = readtable(POOL_FILE,optsP);
row       = P(29826,:);

pulse.n     = row.nPulse;
pulse.tc    = str2double(strsplit(row.tc{1},'|'))';
pulse.tau   = str2double(strsplit(row.tau{1},'|'))';
pulse.A     = str2double(strsplit(row.A{1},'|'))';
pulse.theta = row.theta;


N_all  = 1;          
sigma0 = row.sigma0;
kappa  = row.kappa;
delta  = row.delta;
omega  = row.omega;
gamma  = row.gamma;


beta0_base = row.beta0;
nGrid      = 1000;
beta_grid  = logspace(log10(beta0_base*0.01), log10(beta0_base*10), nGrid);


MAX_CYCLES = 500;
MAX_STEPS  = MAX_CYCLES * nPer;


meanI_periods = nan(nGrid, 100);

for i = 1:nGrid
    beta0 = beta_grid(i);
    fprintf('Scanning β₀ = %.5f  (%d/%d)\n', beta0, i, nGrid);
    Y     = run_full_integrate( ...
              beta0, sigma0, kappa, delta, omega, gamma, ...
              F_vec, dt_hr, N_all, MAX_STEPS, pulse);
    I     = Y(:,3);  

    
    lastIdx  = (MAX_STEPS - 100*nPer + 1) : MAX_STEPS;
    I_tail   = I(lastIdx);

    
    I_mat    = reshape(I_tail, nPer, 100);

    
    meanI_periods(i, :) = mean(I_mat, 1);
end


beta_col   = repelem(beta_grid.', 100);    
period_idx = repmat((1:100).', nGrid, 1);  
meanI_col  = reshape(meanI_periods.', [], 1); 

T = table(beta_col, period_idx, meanI_col, ...
          'VariableNames', {'beta0', 'period', 'mean_I'});
writetable(T, 'experiment/bifurcation_scan_periodic_means.csv');
fprintf('◎  experiment/bifurcation_scan_periodic_means.csv\n');


figure;
C = period_idx;   
sz = 8;           
scatter(beta_col, meanI_col, sz, C, 'filled');
set(gca, 'XScale', 'log');
colormap(flipud(gray(100)));   
caxis([1 100]);                
colorbar;
xlabel('\beta_0 (log scale)');
ylabel('Mean I in each of last 100 periods');
title('Bifurcation Scan: Periodic Means of I');
set(gca, 'FontSize', 12);



function Y = run_full_integrate(beta0,sigma0,kappa,delta,omega,gamma, ...
                                F_vec,dt_hr,N_all,MAX_STEPS,pulse)
    nPer  = numel(F_vec);
    Fstar = (F_vec.^kappa) ./ mean(F_vec.^kappa);
    Jpat  = sparse([1 1 0; 1 1 1; 0 1 1]);

    tvec  = (0:MAX_STEPS) * dt_hr;      
    X0    = [0;0;0];                    
    opts  = odeset('JPattern',Jpat,'RelTol',1e-8,'AbsTol',1e-8);

    [~,Y3] = ode15s(@rhs_lin, tvec, X0, opts);
    E = Y3(:,1);  I = Y3(:,2);  R = Y3(:,3);
    S = N_all - (E+I+R);
    Y = [S,E,I,R];

    function dX = rhs_lin(tHr,X)
        E = max(X(1),0);  I = max(X(2),0);  R = max(X(3),0);
        S = max(N_all - (E+I+R),0);

        idx  = mod(round(tHr/dt_hr), nPer) + 1;
        Ffac = Fstar(idx);
        phi  = sum(pulse.A .* exp(-0.5*((tHr-pulse.tc)./pulse.tau).^2));
        ph   = phi .* Ffac;

        dE = beta0*Ffac*S*I + (1-pulse.theta)*ph*S ...
             - sigma0*Ffac*E - delta*E;
        dI = sigma0*Ffac*E +   pulse.theta*ph*S - gamma*I;
        dR = gamma*I + delta*E - omega*R;
        dX = [dE; dI; dR];
    end
end