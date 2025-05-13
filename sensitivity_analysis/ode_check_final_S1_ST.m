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
S = load("rhythm/posts_uidp_interpolation_alignment_spline_pp.mat");
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

day_idx_const = day_idx(:);
fprintf('\n=== 10) Sobol–Jansen Global Sensitivity  (top-2, new names) ===\n');
parNames = {'beta0','sigma0','theta','kappa', ...
            'gamma','delta','omega','tau','A0'};
idxTheta = [1 2 4 5 6 7 8 9];   
idxA     = 10:20;                
P        = numel(parNames);
N   = 10000;
sobA = scramble(sobolset(P,'Skip',1000,'Leap',200),'MatousekAffineOwen');
sobB = scramble(sobolset(P,'Skip',2000,'Leap',200),'MatousekAffineOwen');

A = net(sobA,N);           % N×P
B = net(sobB,N);
epsRange = 0.5;  logL = log(1-epsRange);  logH = log(1+epsRange);

ThetaA = repmat(theta_refined,N,1);
ThetaB = repmat(theta_refined,N,1);

for j = 1:numel(idxTheta)
    ThetaA(:,idxTheta(j)) = ThetaA(:,idxTheta(j)) + (logL+(logH-logL).*A(:,j));
    ThetaB(:,idxTheta(j)) = ThetaB(:,idxTheta(j)) + (logL+(logH-logL).*B(:,j));
end
ThetaA(:,idxA) = ThetaA(:,idxA) + (logL+(logH-logL).*A(:,end));
ThetaB(:,idxA) = ThetaB(:,idxA) + (logL+(logH-logL).*B(:,end));

ThetaMix = cell(P,1);
for i = 1:P
    Tmix = ThetaA;
    if i<=8, Tmix(:,idxTheta(i)) = ThetaB(:,idxTheta(i));
    else,    Tmix(:,idxA)        = ThetaB(:,idxA);
    end
    ThetaMix{i} = Tmix;
end
calc_metric = @(lam) deal( ...
        sum( accumarray(day_idx_const(day_idx_const<=5), ...
                        lam(day_idx_const<=5), [5,1], ...
                        @(v)sum(maxk(v,min(2,numel(v)))), 0) ), ...
        sum(lam) );

cacheFile = fullfile('sensitivity_analysis','sobol_Y_sample.mat');  
if ~isfolder(fileparts(cacheFile)), mkdir(fileparts(cacheFile)); end

if isfile(cacheFile)
    load(cacheFile, ...
        'Y1_A','Y2_A','Y1_B','Y2_B','Y1_M','Y2_M', ... 
        'ThetaA','ThetaB','ThetaMix','meta');          
else
    Y1_A = zeros(N,1);  Y2_A = zeros(N,1);
    Y1_B = zeros(N,1);  Y2_B = zeros(N,1);
    Y1_M = zeros(N,P);  Y2_M = zeros(N,P);

    for n = 1:N
        [Y1_A(n),Y2_A(n)] = calc_metric( forward(ThetaA(n,:),F,T,dt_hr,false) );
        [Y1_B(n),Y2_B(n)] = calc_metric( forward(ThetaB(n,:),F,T,dt_hr,false) );
    end

    for i = 1:P
        for n = 1:N
            [Y1_M(n,i),Y2_M(n,i)] = calc_metric( forward(ThetaMix{i}(n,:),F,T,dt_hr,false) );
        end
    end
    meta = struct('N',N,'P',P,'timestamp',datetime); 
    save(cacheFile, ...
         'Y1_A','Y2_A','Y1_B','Y2_B','Y1_M','Y2_M', ...
         'ThetaA','ThetaB','ThetaMix','meta', ...
         '-v7.3'); 
end


VarY1 = var([Y1_A;Y1_B],1);
VarY2 = var([Y2_A;Y2_B],1);
S1 = zeros(P,2);   ST = zeros(P,2);
S1_conf = zeros(P,2);  ST_conf = zeros(P,2);
for i = 1:P
%——— S1 / ST  ——%
    d1 = 0.5*mean((Y1_B - Y1_M(:,i)).^2);  
    d2 = 0.5*mean((Y2_B - Y2_M(:,i)).^2);
    S1(i,1) = 1 - d1 / VarY1;
    S1(i,2) = 1 - d2 / VarY2;

    DT1 = 0.5*mean((Y1_A - Y1_M(:,i)).^2);
    DT2 = 0.5*mean((Y2_A - Y2_M(:,i)).^2);
    ST(i,1) = DT1 / VarY1;
    ST(i,2) = DT2 / VarY2;
    %  Var( (Y_B - Y_C)^2 ) / (2 VarY)^2 / N
    se_S1_peak = sqrt(var((Y1_B - Y1_M(:,i)).^2,1) / (4*VarY1^2) / N);
    se_S1_cum  = sqrt(var((Y2_B - Y2_M(:,i)).^2,1) / (4*VarY2^2) / N);
    se_ST_peak = sqrt(var((Y1_A - Y1_M(:,i)).^2,1) / (4*VarY1^2) / N);
    se_ST_cum  = sqrt(var((Y2_A - Y2_M(:,i)).^2,1) / (4*VarY2^2) / N);

    z = 1.96;          % 95 %
    S1_conf(i,1) = z * se_S1_peak;
    S1_conf(i,2) = z * se_S1_cum;
    ST_conf(i,1) = z * se_ST_peak;
    ST_conf(i,2) = z * se_ST_cum;

sum_check = S1(i,1) + ST(i,1);
end  
fprintf('Parameter   S1_peak  ST_peak   S1_cum  ST_cum\n');
for i = 1:P
    fprintf('%-8s   %6.3f   %6.3f   %6.3f   %6.3f\n', ...
        parNames{i}, S1(i,1), ST(i,1), S1(i,2), ST(i,2));
end
clr = [0.30 0.65 0.93];

figs = {'Sobol-S1-Peak','Sobol-ST-Peak','Sobol-S1-Cum','Sobol-ST-Cum'};
titles = {['S1 — ',sprintf('H Top-2 (first 5 d)')], ...
          ['ST — ',sprintf('H Top-2 (first 5 d)')], ...
          ['S1 — ',sprintf('C (20 d)')], ...
          ['ST — ',sprintf('C (20 d)')]};
datas = {S1(:,1), ST(:,1), S1(:,2), ST(:,2)};
ci_arrays = {S1_conf(:,1), ST_conf(:,1), S1_conf(:,2), ST_conf(:,2)};

for k = 1:4
    delete(findobj('Type','figure','Name',figs{k}));
    figure('Name',figs{k},'NumberTitle','off','Color','w',...
           'Units','pixels','Position',[100 100 560 420]);

    bar(datas{k},'FaceColor',clr,'EdgeColor','none');
    hold on                                    
    grid on;  yline(0,'k');
    xticks(1:P); xticklabels(parNames);
    ylabel('Sobol index'); title(titles{k});
    set(gca,'Box','off');
    errorbar(1:P, datas{k}, ci_arrays{k}, 'k',...
             'linestyle','none', 'linewidth',1, 'CapSize',4);

    fname = [figs{k},'.pdf'];
    exportgraphics(gcf,fname,'ContentType','vector',...
                   'BackgroundColor','none');
end