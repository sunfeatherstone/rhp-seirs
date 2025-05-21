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

    Fp        = F .^ kappa;          
    normFac   = mean(Fp);             
    F_star    = Fp ./ normFac;        



    h_loc  = [ 8.25  13.25   9   9  11.5   9   8   9  15.5   9   9 ]; 
    d_idx  = [1,2,3,4,5,6,7,8,9,10,11];
    startHr= 4;
    t_c = h_loc + (d_idx-1)*24 - startHr;        

    
    A_vec = [A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11];
    phi   = zeros(T,1);

    for i = 1:11
        phi = phi + A_vec(i) * exp( -0.5 * ((t_vec_hr - t_c(i))./tau).^2 );
    end

    % 积分
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
    lam_sigma = sigma_vec .* y(2,:)' * dt_hr;       
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
    Nsamp    = 10000;   
    epsRange = 0.5;   
    Psel     = {'beta0','sigma0','theta','kappa', ...
                'gamma','delta','omega','tau','A0'};
    idxTheta = [1 2 4 5 6 7 8 9];
    idxA     = 10:20;

    X_u  = lhsdesign(Nsamp, numel(Psel), 'criterion', 'maximin');
    logL = log(epsRange);
    logH = log(1/epsRange);

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