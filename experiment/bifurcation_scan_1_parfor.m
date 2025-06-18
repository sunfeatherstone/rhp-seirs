
clear; clc;

ckptFile   = 'experiment/bif_scan_ckpt.mat';  
csvFile    = 'experiment/bifurcation_scan_periodic_means_1.csv';
saveStep   = 50;                             
nTail      = 100;                            
dt_hr      = 1;                              
Tper_hr    = 168;                              
nPer       = Tper_hr/dt_hr;                    
MAX_CYCLES = 500;
MAX_STEPS  = MAX_CYCLES * nPer;


Spp    = load('rhythm/preprocessed_data/weibo_spline_pp_98.mat');
pp     = mkpp(Spp.breaks, Spp.coefs);
F_hour = @(tHr) ppval(pp, mod(tHr/2,84));
F_vec  = F_hour((0:nPer-1)'*dt_hr);


POOL_FILE = 'experiment/param_pool.csv';
P         = readtable(POOL_FILE);
row       = P(29826,:);

pulse.n     = row.nPulse;
pulse.tc    = str2double(split(row.tc{1},'|'))';
pulse.tau   = str2double(split(row.tau{1},'|'))';
pulse.A     = str2double(split(row.A{1},'|'))';
pulse.theta = row.theta;

N_all  = 1;
sigma0 = row.sigma0;  kappa = row.kappa; delta = row.delta;
omega  = row.omega;   gamma = row.gamma;
beta0_base = row.beta0;


if isfile(ckptFile)
    load(ckptFile,'beta_grid','meanI_periods','processed');
    fprintf('[Resume]  %d / %d \n', nnz(processed), numel(beta_grid));
else
    nGrid        = 3000;
    beta_grid    = logspace(log10(0.16), log10(0.65), nGrid);

    meanI_periods = nan(nGrid, nTail);   
    processed     = false(1, nGrid);

    save(ckptFile,'beta_grid','meanI_periods','processed','-v7.3');
    fprintf('[Init] β₀=%d\n', nGrid);
end
nGrid = numel(beta_grid);

parpool('local');            

blkSize   = 50;                    
todoIdx   = find(~processed);       
nTodo     = numel(todoIdx);
nBlocks   = ceil(nTodo / blkSize);

for b = 1:nBlocks
    blkRange = todoIdx( (b-1)*blkSize+1 : min(b*blkSize, nTodo) );
    nBlk     = numel(blkRange);

    
    meanBlk  = nan(nBlk, nTail);
    procFlag = false(1, nBlk);

    parfor k = 1:nBlk
        i     = blkRange(k);
        beta0 = beta_grid(i);

        
        

        Y = run_full_integrate(beta0,sigma0,kappa,delta,omega,gamma, ...
                               F_vec,dt_hr,N_all,MAX_STEPS,pulse);

        I  = Y(end-nTail*nPer+1:end,3);
        I_mat          = reshape(I,nPer,nTail);
        meanBlk(k,:)   = mean(I_mat,1);
        procFlag(k)    = true;
    end

    
    meanI_periods(blkRange,:) = meanBlk;
    processed(blkRange)       = procFlag;

    save(ckptFile,'beta_grid','meanI_periods','processed','-v7.3');
    fprintf('[Block %d/%d] save（ %.0f%%）\n', ...
            b, nBlocks, 100*nnz(processed)/numel(processed));
end

delete(gcp('nocreate'));     


if all(processed)
    beta_col   = repelem(beta_grid.', nTail);
    period_idx = repmat((1:nTail).', nGrid, 1);
    meanI_col  = reshape(meanI_periods.',[],1);
    T = table(beta_col,period_idx,meanI_col, ...
              'VariableNames',{'beta0','period','mean_I'});
    writetable(T, csvFile);
    fprintf('◎ %s\n', csvFile);

    delete(ckptFile);
    fprintf('◎ ok\n');
end





function Y = run_full_integrate(beta0,sigma0,kappa,delta,omega,gamma, ...
                                F_vec,dt_hr,N_all,MAX_STEPS,pulse)
    eps = 1e-8;
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
        E = max(X(1),eps);  I = max(X(2),eps);  R = max(X(3),eps);
        
        
        S = max(N_all - (E+I+R),eps);

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