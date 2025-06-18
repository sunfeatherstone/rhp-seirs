
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

targetN = 1e5;           
rows    = cell(targetN,1);
rowID   = 0;

while rowID < targetN
    rate6   = 10.^unifrnd([-5;-5;-5;-5;-5;-1], [1;1;1;1;1;1]);
    b0 = rate6(1); s0 = rate6(2); d0 = rate6(3);
    o  = rate6(4); g  = rate6(5); kappa = rate6(6);

    rho = floquet_rho(b0, s0, g, d0, kappa);
    if ~(rho > 1.01);  continue;  end

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
fid = fopen('experiment/param_pool.csv','w');
fprintf(fid, '%s\n', hdr);
fprintf(fid, '%s\n', rows{:});
fclose(fid);
function rho = local_rho(b0,s0,gamma,delta,kappa,Fvec,dt)
Fstar = (Fvec.^kappa) ./ mean(Fvec.^kappa);
Y = eye(2);
n  = numel(Fvec);
a1 = 0.5 - sqrt(3)/6;
a2 = 0.5 + sqrt(3)/6;
for j = 1:n-1         
    B1 = b0 * Fstar(j);     S1 = s0 * Fstar(j);
    B2 = b0 * Fstar(j+1);   S2 = s0 * Fstar(j+1);

    J1 = [-(S1+delta),  B1 ;  S1 , -gamma];
    J2 = [-(S2+delta),  B2 ;  S2 , -gamma];

    Y  = expm(dt*(a1*J1+a2*J2)) * expm(dt*(a2*J1+a1*J2)) * Y;
end
rho = max(abs(eig(Y)));
end