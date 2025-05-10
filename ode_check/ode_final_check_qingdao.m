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
clear all
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
F_base = @(tau) ppval(pp, mod(tau, 84));   % tau 单位 = 2 h
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


nll = nb_nll_cal(inc, lam_pred, dispersion_k)   % <-- 计算 NLL
eps0_inc = 1e-9;
mask     = inc > eps0_inc;
eps0_w   = 1e-9;
wmape_inc = sum( abs(inc - lam_pred) ) ...
          / max( sum(inc), eps0_w )  * 100;
fprintf('\nIncrement-level WMAPE = %.4f %%\n', wmape_inc);
winList = 2:6;
eps0_w  = 1e-9;

fprintf('\n—— Rolling-Mean Increment WMAPE ——\n');
for w = winList
    filt     = ones(w,1) / w;          % w 桶均值卷积核
    inc_s    = filter(filt,1,inc);
    pred_s   = filter(filt,1,lam_pred);
    validIdx = (1:T)' >= w;
    num   = sum( abs(inc_s(validIdx) - pred_s(validIdx)) );
    den   = sum( inc_s(validIdx) );
    wmape = num / max(den,eps0_w) * 100;

    fprintf('Increment WMAPE  (window = %d) : %.4f %%\n', w, wmape);
end
figure(1); clf;
bar(t_hours,inc); hold on;
plot(t_hours,lam_pred,'r-');

figure(2); clf;
plot(t_hours,cum_obs,'k-'); hold on;
plot(t_hours,cum_pred,'r--');

figure(3); clf;
plot(t_hours,cum_obs-cum_pred,'b-');
lam_pred = forward(theta_refined,F,T,dt_hr,false);
cum_pred = forward(theta_refined,F,T,dt_hr,true);
eps0 = 1e-9;
mape = mean( abs(cum_obs - cum_pred)./max(cum_obs,eps0))*100;
nll  = nb_nll_cal(inc, lam_pred, dispersion_k);   % <-- 计算 NLL
fig1 = figure(1); clf(fig1,'reset');
ax1  = axes(fig1);
cBar = ax1.ColorOrder(1,:);
bar(ax1, t_hours, inc, 'FaceColor', cBar);
hold(ax1,'on');
plot(ax1, t_hours, lam_pred, 'r-');
xlabel(ax1,'t (hour)');  ylabel(ax1,'increment');
legend(ax1, {'observed','predicted'}, 'Location','best');

drawnow;
exportgraphics(fig1,'fig1_increment.pdf','ContentType','vector');
fig2 = figure(2); clf(fig2,'reset');
ax2  = axes(fig2);
plot(ax2, t_hours, cum_obs, ...
     '*', 'Color', cBar, 'MarkerSize', 2, ...
     'DisplayName', 'observed');
hold(ax2,'on');
plot(ax2, t_hours, cum_pred, 'r-', ...
     'LineWidth', 0.8, 'DisplayName', 'fitted');

xlabel(ax2,'t (hour)');
ylabel(ax2,'cumulative');
legend(ax2,'Location','best');

drawnow;
exportgraphics(fig2,'fig2_cumulative.pdf','ContentType','vector');
fig3 = figure(3); clf(fig3,'reset');
ax3  = axes(fig3);
plot(ax3,t_hours,cum_obs-cum_pred,'b-');
xlabel(ax3,'t (hour)'); ylabel(ax3,'cum\_obs - cum\_pred');
drawnow;
exportgraphics(fig3,'fig3_residuals.pdf','ContentType','vector');
b0    = exp(theta_refined(1));   % baseline β₀
s0    = exp(theta_refined(2));   % baseline σ₀
beta1 = exp(theta_refined(5));
Fp      = F .^ beta1;
normFac = mean(Fp);
F_star  = Fp ./ normFac;
beta_t  = b0 * F_star;
sigma_t = s0 * F_star;
idxPlot = t_hours <= 168;
t_sub   = t_hours(idxPlot);
beta_sub  = beta_t(idxPlot);
sigma_sub = sigma_t(idxPlot);
fig4 = figure(4);  clf(fig4,'reset');
ax4  = axes(fig4);

plot(ax4, t_sub, beta_sub,  'b-', 'LineWidth',0.9, 'DisplayName','\beta(t)'); hold(ax4,'on');
plot(ax4, t_sub, sigma_sub, 'r--','LineWidth',0.9, 'DisplayName','\sigma(t)');

xlabel(ax4,'t (hour)');
ylabel(ax4,'rate');
title(ax4,'\beta(t) and \sigma(t) over first 168 hours');
legend(ax4,'Location','best');
xlim(ax4,[0 168]);

drawnow;
exportgraphics(fig4,'fig4_beta_sigma_168h.pdf','ContentType','vector');
inv_beta_sub  = 1 ./ beta_sub;     % 1/β(t)
inv_sigma_sub = 1 ./ sigma_sub;    % 1/σ(t)
fig5 = figure(5);  clf(fig5,'reset');
ax5  = axes(fig5);

plot(ax5, t_sub, inv_beta_sub,  'b-', 'LineWidth',0.9, ...
     'DisplayName','\beta^{-1}(t)');  hold(ax5,'on');
plot(ax5, t_sub, inv_sigma_sub, 'r--','LineWidth',0.9, ...
     'DisplayName','\sigma^{-1}(t)');

xlabel(ax5,'t (hour)');
ylabel(ax5,'reciprocal rate');
title(ax5,'1/\beta(t) and 1/\sigma(t) over first 168 hours');
legend(ax5,'Location','best');
xlim(ax5,[0 168]);

drawnow;
exportgraphics(fig5,'fig5_inv_beta_sigma_168h.pdf','ContentType','vector');