




dt_hr   = 0.5;                               
dt_day  = dt_hr/24;                          
Tmax_step = 1e6;                             


kGrid_tmp = logspace(0, log10(4), 100);
Rper_tmp  = arrayfun(@(k) floquet_multiplier( ...
                    replace_beta_gamma(theta_refined,k), F_hour), kGrid_tmp);
idx0   = find(Rper_tmp<1, 1, 'first');
if isempty(idx0)
    error('k up to 4 still gives Rper>1');
end
k_dag  = kGrid_tmp(idx0)


kList = logspace(log10(k_dag), log10(4), 30);   
nK    = numel(kList);
tDecDay = NaN(1,nK);                            
Rper_vec = zeros(1,nK);


for jj = 1:nK
    k      = kList(jj);
    th_k   = replace_beta_gamma(theta_refined, k);
    Rper_vec(jj) = floquet_multiplier(th_k, F_hour);

    
    tDecDay(jj) = decay_time_I1(th_k, F_hour, dt_hr, Tmax_step) * dt_day;
end


k_lo = 1;  k_hi = 4;
for it = 1:30
    k_mid = sqrt(k_lo*k_hi);
    if sup_Rinst(replace_beta_gamma(theta_refined,k_mid),F_hour) > 1
        k_lo = k_mid;
    else
        k_hi = k_mid;
    end
end
k_Ri1 = k_hi;


figDT = figure('Units','centimeters','Position',[3 3 16 10]); clf;
loglog(kList, tDecDay, '-o', 'LineWidth',1.6,'MarkerSize',4); hold on;
xline(k_Ri1,'--r','LineWidth',1.5);

xlabel('k  (log scale)');
ylabel('time to I<10  (day, log)');     % y 轴也说明是 log
title('Event A · decay time vs. k      β{\_0}\toβ{\_0}/k ,  γ\toγ·k');

grid on; box on;                       
legend({'decay time','Ri = 1'},'Location','northwest');

exportgraphics(figDT,'decay_time_vs_k_loglog.pdf','ContentType','vector');






function th2 = replace_beta_gamma(th, k)
    th2      = th;
    th2(1)   = log(exp(th(1))/k);   
    th2(6)   = log(exp(th(6))*k);   
end


function step_decay = decay_time_I1(theta, F_hour, dt_hr, Nstep_max)
    
    T_sim = Nstep_max;                 
    t_hr  = (0:T_sim-1)' * dt_hr;  
    F     = F_hour(mod(t_hr,168));     
    F     = F / mean(F);

    
    
    b0     = exp(theta(1));
    s0     = exp(theta(2));
    Npop   = exp(theta(3));
    thetaV = 1/(1+exp(-theta(4)));
    kappa  = exp(theta(5));
    gamma  = exp(theta(6));
    delta  = exp(theta(7));
    omega  = exp(theta(8));
    tau    = exp(theta(9));
    A_vec  = exp(theta(10:20));

    
    h_loc  = [ 8.25  13.25   9   9  11.5   9   8   9  15.5   9   9 ];
    d_idx  = 1:11;
    startHr= 4;
    t_c    = h_loc + (d_idx-1)*24 - startHr;
    phi    = zeros(T_sim,1);
    for i = 1:11
        phi = phi + A_vec(i) * exp( -0.5*((t_hr - t_c(i))./tau).^2 );
    end
    F_star = (F.^kappa) / mean(F.^kappa);

    
    epsEI  = 1e-15 * Npop;
    U0     = log([0;0;0]+epsEI);             
    opts   = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxStep',dt_hr);

    
    step_decay = NaN;
    U = U0;
    for kk = 1:Nstep_max
        t_span = [(kk-1)*dt_hr, kk*dt_hr];
        [~,UU] = ode15s(@rhs_log, t_span, U, opts);
        U = UU(end,:)';
        I_now = exp(U(2)) - epsEI;
        if I_now < 10
            step_decay = kk;
            break
        end
    end
    if isnan(step_decay)
        error("...");
    end

    
    function dU = rhs_log(t,U)
        E = max(exp(U(1)) - epsEI, 0);
        I = max(exp(U(2)) - epsEI, 0);
        R = max(exp(U(3)) - epsEI, 0);
        S = max(Npop - (E + I + R), 0);

        idx   = min(floor(t/dt_hr)+1, T_sim);
        Ffac  = F_star(idx);
        phi_t = phi(idx);

        b = b0 * Ffac;
        s = s0 * Ffac;

        dE =  b*S*I/Npop + (1-thetaV)*phi_t*S - s*E - delta*E;
        dI =  s*E        +  thetaV*phi_t*S   - gamma*I;
        dR =  gamma*I + delta*E - omega*R;

        dU = [ dE/(E+epsEI);
               dI/(I+epsEI);
               dR/(R+epsEI) ];
    end
end




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