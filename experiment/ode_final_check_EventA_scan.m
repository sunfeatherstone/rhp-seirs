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

    
    
        epsEI  = 1e-15 * N;                 
        U0     = log([E0; I0; R0] + epsEI); 
        tspan  = (0:T-1) * dt_hr;           

        if final_test
            opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
        else
            opts = odeset('RelTol',1e-3,'AbsTol',1e-3);
        end


        [~,U] = ode15s(@rhs_log, tspan, U0, opts);   

        E = exp(U(:,1)) - epsEI;
        I = exp(U(:,2)) - epsEI;
        R = exp(U(:,3)) - epsEI;
        S = N - (E + I + R);

        y = [S.'; E.'; I.'; R.'];           
    

        
        function dU = rhs_log(t,U)
            E = max(exp(U(1)) - epsEI, 0);
            I = max(exp(U(2)) - epsEI, 0);
            R = max(exp(U(3)) - epsEI, 0);
            S = max(N - (E + I + R),   0);

            idx   = min(floor(t/dt_hr)+1, T); 
            Ffac  = F_star(idx);
            phi_t = phi(idx);

            b = b0 * Ffac;
            s = s0 * Ffac;
            d = delta;
            g = gamma;

            dE =  b*S*I/N + (1-thetaV)*phi_t*S - s*E - d*E;
            dI =  s*E     +   thetaV*phi_t*S   - g*I;
            dR =  g*I + d*E - omega*R;

            dU = [ dE/(E+epsEI);
                dI/(I+epsEI);
                dR/(R+epsEI) ];
        end
    sigma_vec = s0 * F_star;
    lam_phi = thetaV * (phi .* F_star) .* y(1,:)' * dt_hr;
    lam_sigma = sigma_vec .* y(2,:)' * dt_hr;        % E→I 
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
weibo_file = "events/EventA_out.csv";
matFile = 'best_parameter/EventA.mat';
data    = load(matFile, 'theta_refined');
theta_refined = data.theta_refined;
theta_refined
tbl  = readtable(weibo_file);
bucketMin = 30;
S = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp = mkpp(S.breaks, S.coefs);     
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







horizon_days = 20;                           
T20   = horizon_days * 24 / dt_hr;           
t20   = (0:T20-1)' * dt_hr;                 
t_rel20 = mod(base_phase + t20, 168);        
F20   = F_hour(t_rel20);                     
F20   = F20 / mean(F20);                     


startDay  = 11;                              
startIdx  = startDay * 24 / dt_hr + 1;       
tPlot     = t20(startIdx:end);               
Tplot     = numel(tPlot);                    


kList = logspace(0, log10(4), 75);            
nK    = numel(kList);
incMat    = zeros(nK, Tplot);
Rper_vec  = zeros(nK, 1);
Ri_vec    = zeros(nK, 1);


beta_hat  = theta_refined(1);
gamma_hat = theta_refined(6);

for jj = 1:nK
    k   = kList(jj);
    th  = theta_refined;

    th(1) = log(exp(beta_hat) / k);          
    th(6) = log(exp(gamma_hat) * k);         

    
    lam_tmp          = forward(th, F20, T20, dt_hr, false, true);  
    incMat(jj,:)     = lam_tmp(startIdx:end).';    % 只存 12–20 天

    
    Rper_vec(jj) = floquet_multiplier(th, F_hour); 
    Ri_vec(jj)   = sup_Rinst(th, F_hour);          
end


[~, idx_per] = min(abs(Rper_vec - 1));
[~, idx_ri ] = min(abs(Ri_vec   - 1));


figRi = figure('Units','centimeters','Position',[3 3 60 30]); clf;
axRi  = axes(figRi);  hold(axRi,'on');
cmap  = parula(nK);

for jj = 1:nK
    lw = 0.8;  ls = '-';
    if jj == idx_per, lw = 4; end          
    if jj == idx_ri
        lw = 4;
        if jj ~= idx_per, ls = '-'; end      % Ri≈1
    end
    
    plot(axRi, tPlot, incMat(jj,:), ...
         'Color', cmap(jj,:), 'LineWidth', lw, 'LineStyle', ls);
end

xlabel(axRi, 't  (hour)');          % 改标签
ylabel(axRi, 'predicted inc');
title (axRi, 'Event A · β_0/k & γ·k  (k = 1…4, log-uniform) · hour '+string(startDay*24)+'–'+string(horizon_days*24));
grid(axRi,'on');  box(axRi,'on');


xlim(axRi, [min(tPlot)  max(tPlot)]);          
ylim(axRi, [min(incMat(:))  max(incMat(:))]);  


colormap(axRi, parula(nK));                    
caxis(axRi, [log10(kList(1)) log10(kList(end))]);  

cb = colorbar(axRi);                           
cb.Label.String = 'k  (β_0--β_0/k,  γ--γ·k)';

tick_idx   = round(linspace(1, nK, 5));
cb.Ticks   = log10(kList(tick_idx));
cb.TickLabels = arrayfun(@(v)sprintf('%.2f',v), kList(tick_idx), 'uni',0);


cbAx = cb.Axes;             

y_per = log10(kList(idx_per));   
y_ri  = log10(kList(idx_ri));    

kList(idx_per)
kList(idx_ri)

hold(cbAx,'on');
plot(cbAx,[0 1],[y_per y_per],'k-','LineWidth',1.6);   % 实线
plot(cbAx,[0 1],[y_ri  y_ri ],'k-','LineWidth',1.6);  % 虚线
hold(cbAx,'off');

exportgraphics(figRi,'inc_Ri_scan_EventA_tail.pdf','ContentType','vector');



function rho = floquet_multiplier(theta20, F_hour, dt)












if nargin < 3
    dt = 0.05;                 
end


beta0  = exp(theta20(1));
sigma0 = exp(theta20(2));
kappa  = exp(theta20(5));
gamma  = exp(theta20(6));
delta  = exp(theta20(7));

T_tot  = 168;                       
nStep  = round(T_tot / dt);         


t_vec  = (0:nStep) * dt;            
F      = F_hour(t_vec);             
Fp     = F.^kappa;
meanF  = mean(Fp);                  
Fs     = Fp / meanF;                


t_mid  = t_vec(1:end-1) + 0.5*dt;
F_mid  = F_hour(t_mid);
Fs_mid = (F_mid.^kappa) / meanF;    


Phi = eye(2);                       

for k = 1:nStep
    
    b1 = beta0  * Fs(k);            
    s1 = sigma0 * Fs(k);            
    J1 = [-(s1+delta) ,  b1 ;
            s1        , -gamma];

    
    b2 = beta0  * Fs_mid(k);
    s2 = sigma0 * Fs_mid(k);
    J2 = [-(s2+delta) ,  b2 ;
            s2        , -gamma];

    
    b3 = beta0  * Fs(k+1);
    s3 = sigma0 * Fs(k+1);
    J3 = [-(s3+delta) ,  b3 ;
            s3        , -gamma];

    K1 = J1 * Phi;
    K2 = J2 * (Phi + 0.5*dt*K1);
    K3 = J2 * (Phi + 0.5*dt*K2);    
    K4 = J3 * (Phi + dt*K3);

    Phi = Phi + (dt/6)*(K1 + 2*K2 + 2*K3 + K4);
end

rho = max(abs(eig(Phi)));
end


function Ri = sup_Rinst(theta20, F_hour, dt)




if nargin < 3,  dt = 0.05;  end


beta0  = exp(theta20(1));
sigma0 = exp(theta20(2));
kappa  = exp(theta20(5));
gamma  = exp(theta20(6));
delta  = exp(theta20(7));

T_tot = 168;
nStep = round(T_tot/dt);
t_hr  = (0:nStep-1) * dt;      
F     = F_hour(t_hr);
Fstar = (F.^kappa) / mean(F.^kappa);

sigma_t = sigma0 * Fstar;
beta_t  = beta0  * Fstar;

R_inst  = (beta_t .* sigma_t) ./ ((sigma_t + delta) * gamma);
Ri      = max(R_inst);
end