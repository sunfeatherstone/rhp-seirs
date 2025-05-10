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
function [theta_refined, t_hours, inc, cum_obs, lam_pred, cum_pred] = stagewise_fit(lb,ub)
    clc;  close all;
    format long g
    dispersion_k = 5;
    weibo_file = "cascades/Qingdao_out.csv";
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
    resfun_full = @(th) ...
        nb_dev_res( inc , forward(th,F,T,dt_hr,false) , dispersion_k );

    psoOpts = optimoptions('particleswarm', 'SwarmSize',60,'MaxIterations',100,'Display','iter');
    obj_pso = @(th) sum(resfun_full(th).^2);
    theta0_all    = particleswarm(obj_pso, length(lb), lb, ub, psoOpts);
    lsqOpts = optimoptions('lsqnonlin', ...
    'Algorithm','trust-region-reflective', ...  % TRR 是默认，显式写出来方便改
    'StepTolerance',1e-36, ...
    'OptimalityTolerance',1e-36, ...
    'FunctionTolerance',1e-36, ...
    'ScaleProblem','jacobian', ...
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',1e9, ...
    'Display','iter');
    theta_refined = lsqnonlin(resfun_full, theta0_all, lb, ub, lsqOpts);
    lam_pred = forward(theta_refined,F,T,dt_hr,false);
    cum_pred = forward(theta_refined,F,T,dt_hr,true);
    if any(lam_pred < -1e-2)
        idx = find(lam_pred < 0, 1);
        error("lamPred:Negative", ...
            "lam\\_pred contains negative value at index %d (%.3g)", ...
            idx, lam_pred(idx));
    end

    eps0 = 1e-9;
    mape = mean( abs(cum_obs - cum_pred)./max(cum_obs,eps0))*100;
    nll = nb_nll_cal(inc, lam_pred, dispersion_k);   % <-- 计算 NLL

    fn = sprintf('CMA-PSO-LSQ_pipiline/MAPE_%05.2f_NLL_%08.1f', mape, nll);
    save([fn '.mat'], 'theta_refined','mape');

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
clear; clc; close all;
lb_11  = [-7; -7; log(1e4); -4; log(0.1); log(1/(7*24));log(1/(7*24)); log(1/(30*24)); log(0.1)];
ub_11  = [ 2;  2; log(1e6);  4; log(5);   log(1); log(1);        log(1/12);      log(12)];
lb_p   = -7*ones(11,1);     % A1…A11
ub_p   =  3*ones(11,1);
D      = numel(lb_p);

dispersion_k = 5;
maxRestart        = 10000;
budgetPerRestart  = 40*D*2;
lambda            = 4 + floor(3*log(D));
stateF            = 'cma_state_local_new.mat';
if isfile(stateF)
    load(stateF,"bestx","bestFit","restartID");
else
    bestFit   = inf;
    bestx     = lb_p' + rand(1,D).*(ub_p'-lb_p');
    restartID = 1;
end
while restartID<=maxRestart
    fprintf('\n===== Restart %d / %d =====\n',restartID,maxRestart);
    shrink = 0.5^(restartID-1);
    m0     = bestx + shrink*randn(1,D).*(ub_p'-lb_p');
    m0     = max(min(m0,ub_p'),lb_p');      % clip
    sigma0 = shrink * 0.3*mean(ub_p-lb_p);
    opts = struct('lambda',lambda,...
                  'maxEval',budgetPerRestart,...
                  'bounds',[lb_p';ub_p'],...
                  'verb_disp',10);

    [xmin,fmin] = simple_cmaes(@(x)evaluate_stage(x,lb_11,ub_11,dispersion_k),...
                               m0',sigma0,opts);

    fprintf('  Restart %d best NLL = %.4g\n',restartID,fmin);

    if fmin<bestFit
        bestFit=fmin; bestx=xmin';
        fprintf('  >>> New global best! NLL = %.4g <<<\n',bestFit);
    end

    save(stateF,"bestFit","bestx","restartID");
    restartID = restartID+1;
end

fprintf('\n========= DONE =========\n');
function nll=evaluate_stage(A_log,lb11,ub11,disp_k)
    row=A_log(:); lb=[lb11;row]; ub=[ub11;row];
    try
        [theta_refined, t_hours, inc, cum_obs, lam_pred, cum_pred] = stagewise_fit(lb,ub);
        nll=nb_nll_cal(inc,lam_pred,disp_k);
        if ~isfinite(nll), nll=inf; 
        
        end
    catch
        nll=inf;
    end
end




function [xbest, fbest] = simple_cmaes(obj, x0, sigma, opts)
    
    D = numel(x0);
    lambda = get_opt(opts,'lambda',4+floor(3*log(D)));
    mu     = get_opt(opts,'mu',floor(lambda/2));
    w      = log(mu+0.5) - log(1:mu)';
    w      = w / sum(w);                  % normalize
    mueff  = 1/sum(w.^2);
    
    maxEval = get_opt(opts,'maxEval',1e4);
    bound   = get_opt(opts,'bounds',[-inf(1,D);inf(1,D)]);
    verb    = get_opt(opts,'verb_disp',20);
    
    cc   = 4/(D+4);             % cumulation for C
    cs   = (mueff+2)/(D+mueff+3);
    c1   = 2/((D+1.3)^2+mueff);
    cmu  = min(1-c1, 2*(mueff-2+1/mueff)/((D+2)^2+mueff));
    damps= 1 + 2*max(0,sqrt((mueff-1)/(D+1))-1) + cs;
    
    pc = zeros(D,1); ps = zeros(D,1);
    B = eye(D); Ddiag = ones(D,1); C = eye(D);
    eigeneval = 0; chiN = D^0.5*(1-1/(4*D)+1/(21*D^2));
    
    mean = x0(:);
    xbest = mean; fbest = obj(reflect(mean,bound));
    
    counteval = 1;
    
    while counteval < maxEval
        arz  = randn(D,lambda);
        ary  = B*(Ddiag .* arz);           % D.*N(0,1)
        arx  = reflect(mean + sigma*ary, bound);
        fitness = zeros(1,lambda);
        for k=1:lambda
            fitness(k) = obj(arx(:,k));
        end
        counteval = counteval + lambda;
        [fmin,idx] = min(fitness);
        if fmin < fbest
            fbest=fmin; xbest=arx(:,idx);
        end
        [~,ord] = sort(fitness);
        xsel = arx(:,ord(1:mu));
        zsel = arz(:,ord(1:mu));
        old_mean = mean;
        mean = xsel * w;
        ps = (1-cs)*ps + sqrt(cs*(2-cs)*mueff)/(sigma) * (B*zsel* w);
        hsig = norm(ps)/sqrt(1-(1-cs)^(2*counteval/lambda))/chiN < 1.4+2/(D+1);
        pc = (1-cc)*pc + hsig*sqrt(cc*(2-cc)*mueff)/(sigma) * (B*Ddiag.*zsel* w);
        artmp = (1/sigma) * (xsel - old_mean) ;
        C = (1-c1-cmu)*C + c1*(pc*pc'+ (1-hsig)*cc*(2-cc)*C) ...
            + cmu * artmp * diag(w) * artmp';
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        if counteval - eigeneval > lambda/(c1+cmu)/D/10
            eigeneval = counteval;
            [B,tmp] = eig((C+C')/2);
            Ddiag  = sqrt(max(diag(tmp),1e-32));
        end
        if verb && mod(counteval,verb*lambda)<lambda
            fprintf('eval=%d  fbest=%.3g  step=%.3g\n',...
                     counteval,fbest,sigma);
        end
    end
    end % ======= END simple_cmaes =======
    function val=get_opt(s,name,default)
        if isfield(s,name), val=s.(name); else, val=default; end
    end
    function X=reflect(X,b)
        lb=b(1,:); ub=b(2,:);
        if all(isinf(lb) & isinf(ub)), return; end
        X = max(min(X,ub'),lb');
    end





