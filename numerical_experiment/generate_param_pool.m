
clear; clc;
rng(20250527, 'twister');  

dt_hr        = 0.5;        
T_period_hr  = 168;        
nStepPerPer  = T_period_hr / dt_hr;


S   = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp  = mkpp(S.breaks, S.coefs);
F_hour  = @(t_hr) ppval(pp, mod(t_hr/2, 84));   % t→F(t)
tGrid   = (0:nStepPerPer-1)' * dt_hr;
F_vec   = F_hour(tGrid);                       


vec2str  = @(v) sprintf('%.4g|', v);       
stripEnd = @(s) s(1:end-(~isempty(s)));     

floquet_rho = @(b0,s0,g,d,k) local_rho(b0,s0,g,d,k,F_vec,dt_hr);

targetN = 1e4;           
rows    = cell(targetN,1);
rowID   = 0;

while rowID < targetN
    rate6   = 10.^unifrnd([-3;-4;-4;-4;-5;-1], [1;0;0;0;-1;1]);
    b0 = rate6(1); s0 = rate6(2); d0 = rate6(3);
    o  = rate6(4); g  = rate6(5); kappa = rate6(6);

    rho = floquet_rho(b0, s0, g, d0, kappa);
    if ~(rho > 1.01 & rho < 1e6);  continue;  end

    nPulse = max(1, min(100, round(10.^unifrnd(0,2))));
    tc   = 30*24*rand(nPulse,1);                      % 0–30
    tau  = 10.^unifrnd(log10(1/6), log10(12), [nPulse,1]);  % 10 min–12 h
    Aamp = 10.^unifrnd(-4, 2, [nPulse,1]);            % 1e-4–1e2
    theta = 1/(1 + exp(-unifrnd(-5,5)));              % logistic

    rowID = rowID + 1;
    rows{rowID} = sprintf([ ...
        '%d,%.4e,%.4e,%.4e,%.4e,%.4e,%.4e,%.5f,%d,%.5f,%.0f,"%s","%s","%s"'], ...
        rowID,b0,s0,kappa,d0,o,g,rho, ...
        nPulse,theta,1,...                     
        stripEnd(vec2str(tc)), ...
        stripEnd(vec2str(tau)), ...
        stripEnd(vec2str(Aamp)));
    
    rowID
end
hdr = ['paramID,beta0,sigma0,kappa,delta,omega,gamma,rho,' ...
       'nPulse,theta,Nall,tc,tau,A'];
fid = fopen('numerical_experiment/param_pool.csv','w');
fprintf(fid, '%s\n', hdr);
fprintf(fid, '%s\n', rows{:});
fclose(fid);
function rho = local_rho(b0,s0,gamma,delta,kappa,Fv,dt)
Fstar = (Fv.^kappa) ./ mean(Fv.^kappa);
Y = eye(2);
for jj = 1:numel(Fv)
    B = b0 * Fstar(jj);
    S = s0 * Fstar(jj);
    J = [-(S+delta) , B ; S , -gamma];
    Y = expm(J*dt) * Y;      
end
rho = max(abs(eig(Y)));
end