clear functions

function red = shrink_theta(theta_full,stage)
    idx = (stage-1)*3 + (12:14);
    red = theta_full(idx);
end

function full = expand_theta(theta_red,theta_base,stage)
    full = theta_base;
    idx = (stage-1)*3 + (12:14);
    full(idx) = theta_red;
end
function r = resid_stage(th_red, stage, theta_base, ...
    F, T, dt_hr, inc, idx, dispersion_k)
    th_full = expand_theta(th_red, theta_base, stage);
    C_full  = forward(th_full, F, T, dt_hr, true);
    lam_full = forward(th_full, F, T, dt_hr, false);
    r = nb_dev_res( inc(idx) , lam_full(idx) , dispersion_k );
end

function r = nb_dev_res(y, mu, k)
    eps = 1e-9;
    mu = max(mu,eps); k = max(k,eps);
    d2 = 2*( y.*log(max(y,eps)./mu) - (y+k).*log((y+k)./(mu+k)) );
    r  = sign(y-mu) .* sqrt(max(d2,0));
end

function out = forward(theta,F,T,dt_hr,cumFlag,final_test)
    b0     = exp(theta(1)); 

    s0     = exp(theta(2)); 

    N     = exp(theta(3)); E0     = 0;
    I0     = 0; R0     = 0;
    S0      = N;
    thetaV       = 1/(1+exp(-theta(4)));;

    kappa  = exp(theta(5));
    gamma     = exp(theta(6)); delta=exp(theta(7));
    omega    = exp(theta(8));

    T0 = datetime;  
    persistent t_origin day_idx_global
    if isempty(t_origin)
        vars = evalin('base', 'whos');
        t_origin = evalin('base','t0');
        day_idx_global = evalin('base','day_idx');
    end
    t_vec_hr = (0:T-1)' * dt_hr;

    tau = exp(theta(9));
    A1   = exp(theta(10)); 
    A2   = exp(theta(11)); 
    A3   = exp(theta(12)); 
    A4   = exp(theta(13)); 
    A5   = exp(theta(14)); 
    A6   = exp(theta(15)); 
    A7   = exp(theta(16)); 
    A8   = exp(theta(17)); 
    A9   = exp(theta(18)); 
    A10   = exp(theta(19)); 
    A11   = exp(theta(20)); 



    phi = zeros(T,1);
    dist1 = abs(mod(t_vec_hr-(8.25-4)+12,24)-12);
    dist2 = abs(mod(t_vec_hr-(13.25-4)+12,24)-12);
    dist3 = abs(mod(t_vec_hr-(9-4)+12,24)-12);
    dist4 = abs(mod(t_vec_hr-(9-4)+12,24)-12);
    dist5 = abs(mod(t_vec_hr-(11.5-4)+12,24)-12);
    dist6 = abs(mod(t_vec_hr-(9-4)+12,24)-12);
    dist7 = abs(mod(t_vec_hr-(8-4)+12,24)-12);
    dist8 = abs(mod(t_vec_hr-(9-4)+12,24)-12);
    dist9 = abs(mod(t_vec_hr-(15.5-4)+12,24)-12);
    dist10 = abs(mod(t_vec_hr-(9-4)+12,24)-12);
    dist11 = abs(mod(t_vec_hr-(9-4)+12,24)-12);
    idx1 = (day_idx_global == 1);
    idx2 = (day_idx_global == 2);
    idx3 = (day_idx_global == 3);
    idx4 = (day_idx_global == 4);
    idx5 = (day_idx_global == 5);
    idx6 = (day_idx_global == 6);
    idx7 = (day_idx_global == 7);
    idx8 = (day_idx_global == 8);
    idx9 = (day_idx_global == 9);
    idx10 = (day_idx_global == 10);
    idx11 = (day_idx_global == 11);
    phi(idx1) = A1 * exp(-0.5*(dist1(idx1)./tau).^2);
    phi(idx2) = A2 * exp(-0.5*(dist2(idx2)./tau).^2);
    phi(idx3) = A3*exp(-0.5*(dist3(idx3)./tau).^2);
    phi(idx4) = A4*exp(-0.5*(dist4(idx4)./tau).^2);
    phi(idx5) = A5*exp(-0.5*(dist5(idx5)./tau).^2);
    phi(idx6) = A6*exp(-0.5*(dist6(idx6)./tau).^2);
    phi(idx7) = A7*exp(-0.5*(dist7(idx7)./tau).^2);
    phi(idx8) = A8*exp(-0.5*(dist8(idx8)./tau).^2);
    phi(idx9) = A9*exp(-0.5*(dist9(idx9)./tau).^2);
    phi(idx10) = A10*exp(-0.5*(dist10(idx10)./tau).^2);
    phi(idx11) = A11*exp(-0.5*(dist11(idx11)./tau).^2);

    Fp        = F .^ kappa;           
    normFac   = mean(Fp);             
    F_star    = Fp ./ normFac;       
    y = nan(4,T); y(:,1) = [S0;E0;I0;R0];
    for k=1:T-1
        b = b0*F_star(k);
        g = gamma;
        d = delta;
        s = s0*F_star(k);
        ph = phi(k)*F_star(k);
        S=y(1,k); E=y(2,k); I=y(3,k); R=y(4,k);
        dS= -b*S*I/N + omega*R - ph*S;
        dE=  b*S*I/N  + (1-thetaV)*ph*S - s*E - d*E;
        dI=  s*E + thetaV*ph*S - g*I ;
        dR=  g*I + d*E - omega*R;
        y(:,k+1) = y(:,k) + dt_hr*[dS;dE;dI;dR];
    end
    sigma_vec = s0 * F_star;
    lam_phi = thetaV * (phi .* F_star) .* y(1,:)' * dt_hr;
    lam_sigma = sigma_vec .* y(2,:)' * dt_hr;        % E→I 流
    lam = lam_sigma + lam_phi;

    if cumFlag
        out = cumsum(lam);
    else
        out = lam;
    end
end

function nll = nb_nll_cal(y, mu, k)
    
        eps = 1e-9;
        mu = max(mu, eps);
        k  = max(k,  eps);
        ll = gammaln(y + k) - gammaln(k) - gammaln(y + 1) ...
           + k .* log(k ./ (k + mu)) ...
           + y .* log(mu ./ (k + mu));
    
        nll = -sum(ll);
end

clc;  close all;
clear all;
format long g
dispersion_k = 5;
weibo_file = "cascades/Qingdao_out.csv";
matFile = 'best_parameter/Qingdao.mat';
data    = load(matFile, 'theta_refined');
theta_refined = data.theta_refined;
theta_refined
tbl  = readtable(weibo_file);
bucketMin = 30;
S = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp = mkpp(S.breaks, S.coefs);     % true cubic spline on [0,84]
F_base = @(tau) ppval(pp, mod(tau, 84));  
F_hour = @(t_hr) F_base(t_hr/2);
t30 = (0:0.5:167.5)';          % 336 × 1
ts = datetime(tbl.created_at,'InputFormat','yyyy-MM-dd HH:mm:ss');
t0_raw   = min(ts);
tEnd_raw = max(ts);
t0   = datetime(year(t0_raw),month(t0_raw),day(t0_raw), ...
                hour(t0_raw),floor(minute(t0_raw)/bucketMin)*bucketMin,0);
tEnd = datetime(year(tEnd_raw),month(tEnd_raw),day(tEnd_raw), ...
                hour(tEnd_raw),ceil(minute(tEnd_raw)/bucketMin)*bucketMin,0);
dt      = minutes(bucketMin);
tvec    = (t0:dt:tEnd)';
inc     = histcounts(ts, [tvec; tEnd+dt])';
cum_obs = cumsum(inc);
T       = numel(inc);
dt_hr   = minutes(dt)/60;
base_phase = mod(weekday(t0)-2,7)*24 + hour(t0) + minute(t0)/60;
t_rel  = mod(base_phase + (0:T-1)*dt_hr, 168);
F = F_hour(t_rel);
F = F(:);
F = F / mean(F);
t_hours = (0:T-1)' * dt_hr;
t_abs   = t0 + hours(t_hours);
day_idx =  1 + floor(days(t_abs - dateshift(t0,'start','day')));
assignin('base','t0',     t0);
assignin('base','day_idx', day_idx);
lam_pred = forward(theta_refined,F,T,dt_hr,false);
cum_pred = forward(theta_refined,F,T,dt_hr,true);
eps0 = 1e-9;
mape = mean( abs(cum_obs - cum_pred)./max(cum_obs,eps0))*100


nll = nb_nll_cal(inc, lam_pred, dispersion_k)  
Psel     = {'beta0','sigma0','theta','kappa', ...
            'gamma','delta','omega','tau','A0'};   
idxTheta = [1 2 4 5 6 7 8 9];    
idxA     = 10:20;                 
Nsamp    = 2000;
epsRange = 0.20;
dayCut   = 5;
mask5 = day_idx <= dayCut;
sumTop2_5 = zeros(Nsamp,1);
cumTot    = zeros(Nsamp,1);
fprintf('\n=== PRCC (no lock) — sumTop2Peak & cumTotal ====================\n');

Psel      = {'beta0','sigma0','theta','kappa', ...
             'gamma','delta','omega','tau','A0'};  
idxTheta  = [1 2 4 5 6 7 8 9];
idxA      = 10:20;
prccSampleFile = 'sensitivity_analysis/sample/prcc_sample.mat'; 
if isfile(prccSampleFile)
    load(prccSampleFile,'ThetaS','X_u','Nsamp','epsRange');
else
    Nsamp    = 5000;   
    epsRange = 0.20;   
    Psel     = {'beta0','sigma0','theta','kappa', ...
                'gamma','delta','omega','tau','A0'};
    idxTheta = [1 2 4 5 6 7 8 9];
    idxA     = 10:20;

    X_u  = lhsdesign(Nsamp, numel(Psel), 'criterion', 'maximin');
    logL = log(1 - epsRange);
    logH = log(1 + epsRange);

    ThetaS = repmat(theta_refined, Nsamp, 1);
    for j = 1:numel(idxTheta)
        ThetaS(:, idxTheta(j)) = ThetaS(:, idxTheta(j)) + ...
            (logL + (logH-logL).*X_u(:, j));
    end
    ThetaS(:, idxA) = ThetaS(:, idxA) + ...
            (logL + (logH-logL).*X_u(:, end));

    save(prccSampleFile, 'ThetaS','X_u','Nsamp','epsRange', '-v7.3');
end


sumTop2_5  = zeros(Nsamp,1);
cumTot     = zeros(Nsamp,1);

for n = 1:Nsamp
    lam = forward(ThetaS(n,:), F, T, dt_hr, false);
    dailyTop2 = accumarray(day_idx(mask5), lam(mask5), [dayCut,1], ...
                 @(v) sum(maxk(v,2)), 0);
    sumTop2_5(n) = sum(dailyTop2);
    cumTot(n)    = sum(lam);
end
dof  = Nsamp - numel(Psel) - 1;
t95  = tinv(0.975,dof);
PCOR = @(Y) partialcorr([tiedrank(X_u) tiedrank(Y)],'type','Spearman');

[Rmat, Pmat] = PCOR(sumTop2_5);      
R9_peak = Rmat(1:end-1,end).';       
p9_peak = Pmat(1:end-1,end).';
CI9_peak= t95 .* sqrt((1-R9_peak.^2)./dof);

[Rmat, Pmat] = PCOR(cumTot);
R9_tot  = Rmat(1:end-1,end).';
p9_tot  = Pmat(1:end-1,end).';
CI9_tot = t95 .* sqrt((1-R9_tot.^2)./dof);
fprintf('\nPRCC (no lock) — Y = sumTop2Peak\n');
for k = 1:numel(Psel)
    fprintf('%-7s  ρ=%+6.3f  (p=%.2g)\n', Psel{k}, R9_peak(k), p9_peak(k));
end
fprintf('\nPRCC (no lock) — Y = cumTotal, 20天\n');
for k = 1:numel(Psel)
    fprintf('%-7s  ρ=%+6.3f  (p=%.2g)\n', Psel{k}, R9_tot(k),  p9_tot(k));
end
figs = {'PRCC-sumTop2Peak','PRCC-cumTotal'};
figName = 'PRCC-sumTop2Peak';
delete(findobj('Type','figure','Name',figName));
f1=figure('Name',figName,'NumberTitle','off','Color','w',...
       'NextPlot','replacechildren');

bar(R9_peak,'FaceColor',[0.20 0.60 0.90],'EdgeColor','none'); hold on
errorbar(1:numel(R9_peak), R9_peak, CI9_peak,...
         'k.','CapSize',6,'LineWidth',1);
yline(0,'k'); grid on
xticks(1:numel(Psel));  xticklabels(Psel);
ylabel('PRCC');  title({'PRCC';'(sumTop2Peak, first 5 days)'});
exportgraphics(f1,[figs{1} '.pdf'], ...
               'ContentType','vector','BackgroundColor','none')
hold off
figName = 'PRCC-cumTotal';
delete(findobj('Type','figure','Name',figName));

f2=figure('Name',figName,'NumberTitle','off','Color','w', ...
       'NextPlot','replacechildren');

bar(R9_tot,'FaceColor',[0.20 0.60 0.90],'EdgeColor','none'); hold on
errorbar(1:numel(R9_tot), R9_tot, CI9_tot, ...
         'k.','CapSize',6,'LineWidth',1);              % 95 % CI

yline(0,'k'); grid on
xticks(1:numel(Psel));  xticklabels(Psel);
xlim([0.5 numel(Psel)+0.5])
ylabel('PRCC');
title({'PRCC'; '(cumTotal, 20 days)'});
exportgraphics(f2,[figs{2} '.pdf'], ...
               'ContentType','vector','BackgroundColor','none')
hold off