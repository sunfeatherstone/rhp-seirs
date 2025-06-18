clear;  close all;  clc;


matEvent = 'best_parameter/EventA.mat';    
load(matEvent,'theta_refined');    

splineFile = 'rhythm/preprocessed_data/weibo_spline_pp_98.mat';
S = load(splineFile);                        
pp = mkpp(S.breaks, S.coefs);                
F_base  = @(tau) ppval(pp, mod(tau,84));     
F_hour  = @(h) F_base(h/2);                  


b0      = exp(theta_refined(1));
s0      = exp(theta_refined(2));
Npop    = exp(theta_refined(3));
thetaV  = 1/(1+exp(-theta_refined(4)));
kappa   = exp(theta_refined(5));
gamma   = exp(theta_refined(6));
deltaE  = exp(theta_refined(7));
omega   = exp(theta_refined(8));
tauG    = exp(theta_refined(9));             
A1      = exp(theta_refined(10));            


dt_ref = 0.25;                               
t_ref  = 0:dt_ref:168-dt_ref;                
Fp     = F_hour(t_ref).^kappa;
meanFp = mean(Fp);                           
driver = @(t) (F_hour(t).^kappa) / meanFp;   


dt_hr  = 0.5;                                
dt_hr_1  = 0.01; 
t0_vec = 0:dt_hr:168;                        
nT0    = numel(t0_vec);

cum_lin  = zeros(1,nT0);                     
cum_full = zeros(1,nT0);                     

optsODE = odeset('RelTol',1e-6,'AbsTol',1e-9,...
                 'MaxStep',dt_hr/2);

for kk = 1:nT0
    t0 = t0_vec(kk);                         
    
    pulse = @(t) A1 * exp(-0.5 * ((t - t0)/tauG).^2) .* driver(t);

    
    t_span = t0:dt_hr_1:t0+24-dt_hr_1;           
    nt     = numel(t_span);
    
    
    odeLin = @(t,Y) rhs_lin(t,Y,pulse,driver,...
                            s0,deltaE,gamma,thetaV,Npop);

    Y0   = [0; 0; 0];                               
    tEnd = t0 + 24;
    sol  = ode15s(odeLin,[t0 tEnd],Y0,...
                odeset('RelTol',1e-6,'AbsTol',1e-9));

    cum_lin(kk) = sol.y(3,end);                     

    
    odeFun = @(t,X) rhs_full(t,X,pulse,driver,...
                  b0,s0,gamma,deltaE,omega,thetaV,Npop);
    X0     = [Npop;0;0;0];                   
    sol    = ode15s(odeFun,[t0 t0+24],X0,optsODE);
    tt = (t0:dt_hr:t0+24);                   
    YY = deval(sol,tt);
    S = YY(1,:);   E = YY(2,:);   I = YY(3,:);
    drv_tt  = driver(tt);
    sig_tt  = s0 * drv_tt;
    lam_phi = thetaV * pulse(tt) .* S;       
    lam_sig = sig_tt .* E;
    cum_full(kk) = trapz(tt, lam_phi + lam_sig);
end


mu_lin  = mean(cum_lin);          
mu_full = mean(cum_full);

cum_lin_n  = cum_lin  ./ mu_lin;  
cum_full_n = cum_full ./ mu_full;


fig = figure('Units','centimeters','Position',[3 3 34 17]); clf;
plot(t0_vec, cum_lin_n ,'b-' , 'LineWidth',1.4); hold on;
plot(t0_vec, cum_full_n,'r-','LineWidth',1.4);


for d = 0:7
    x = d*24;
    line([x x],[0 max(cum_lin_n)],'Color',[.7 .7 .7],'LineStyle','--');
end
xticks(0:24:168);
xticklabels({'Mon','Tue','Wed','Thu','Fri','Sat','Sun','Mon'});
ylabel('normalised 24h cumulative increment');
legend({'DFE-linear / mean','full RhP-SEIRS / mean'},'Location','northwest');
legend({'DFE‑linear approx.','full RhP‑SEIRS'},'Location','northwest');
title('Event A · 24h cumulative increment vs. pulse centre');
grid on;  box on;
ylim([min(min(cum_lin_n),min(cum_full_n))   max(max(cum_lin_n),max(cum_full_n))]); 
xlim([0-0.001   168+0.001]); 

exportgraphics(fig,'pulse_centre_scan.pdf','ContentType','vector');


function dX = rhs_full(t,X,pulse,driver,...
                       b0,s0,gamma,deltaE,omega,thetaV,N)
    S = X(1);  E = X(2);  I = X(3);  R = X(4);

    drv = driver(t);
    beta = b0 * drv;
    sigma= s0 * drv;
    P    = pulse(t);

    dS = -beta*S*I/N  + omega*R - P*S;
    dE =  beta*S*I/N  + (1-thetaV)*P*S - sigma*E - deltaE*E;
    dI =  sigma*E      + thetaV*P*S   - gamma*I;
    dR =  gamma*I + deltaE*E - omega*R;
    dX = [dS; dE; dI; dR];
end


function dY = rhs_lin(t,Y,pulse,driver,...
                      s0,deltaE,gamma,thetaV,N)
    E = Y(1);  I = Y(2);

    drv   = driver(t);
    sig   = s0 * drv;
    P     = pulse(t);

    lam_phi = thetaV * P * N;           
    lam_sig = sig * E;

    dE = (1-thetaV)*P*N - (sig + deltaE)*E;
    dI = lam_phi + lam_sig - gamma*I;
    dC = lam_phi + lam_sig;             

    dY = [dE; dI; dC];
end