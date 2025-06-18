
clear; clc;

POOL_FILE   = 'experiment/param_pool.csv';
RESULT_FILE = 'experiment/sampling_results.csv';
CHUNK       = 10;               
dt_hr       = 0.5;
T_period_hr = 168;
burnin      = 5;
chkCycles   = 5;
tolPer      = 1e-4;
maxSteps    = 1e6;
k_check     = 3;                 
N_all       = 1;                

S   = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp  = mkpp(S.breaks,S.coefs);
F_hour  = @(t_hr) ppval(pp, mod(t_hr/2, 84));
nStepPerPer  = T_period_hr / dt_hr;
F_vec = F_hour((0:nStepPerPer-1)'*dt_hr);

opts = detectImportOptions(POOL_FILE, 'Delimiter', ',');
P    = readtable(POOL_FILE, opts);            
nTot = height(P);

doneMask = false(nTot,1);
if isfile(RESULT_FILE)
    opts2 = detectImportOptions(RESULT_FILE, 'Delimiter', ',');
    R     = readtable(RESULT_FILE, opts2);
    doneMask(R.paramID) = true;                
    fprintf('✔  %d \n', nnz(doneMask));
else
    hdr = ['paramID,flag,steps'];
    fid = fopen(RESULT_FILE,'w'); fprintf(fid,'%s',hdr); fclose(fid);
end

fid = fopen(RESULT_FILE, 'a');          
cleanupObj = onCleanup(@() fclose(fid));

buf = {};  bufCnt = 0;                 

for idx = find(~doneMask)'              
    b0   = P.beta0(idx);  s0  = P.sigma0(idx);
    kappa= P.kappa(idx);  d0  = P.delta(idx);
    omega= P.omega(idx);  g   = P.gamma(idx);

    nP = P.nPulse(idx);
    tc   = str2double(strsplit(P.tc{idx},  '|'))';
    tau  = str2double(strsplit(P.tau{idx}, '|'))';
    Aamp = str2double(strsplit(P.A{idx},  '|'))';
    theta= P.theta(idx);

    pulse = struct('n', nP, 'tc', tc, 'tau', tau, ...
                   'A', Aamp, 'theta', theta);

    [flag, nSteps] = simulate( ...
        b0, s0, kappa, d0, omega, g, F_vec, dt_hr, T_period_hr, N_all, ...
        burnin, chkCycles, tolPer, maxSteps, pulse, k_check);

    bufCnt = bufCnt + 1;
    buf{bufCnt} = sprintf('\n%d,%s,%d', ...
            idx, flag, nSteps);  %#ok<*SAGROW>

    fprintf('[%4d/%4d]  %-9s  steps=%d\n', ...
            idx, nTot, flag, nSteps);

    if bufCnt >= CHUNK
        fprintf(fid, '%s', buf{:});
        buf = {}; bufCnt = 0;
        fseek(fid,0,'cof'); 
    end
end

if bufCnt > 0
    fprintf(fid, '%s', buf{:});
    fseek(fid,0,'cof'); 
end




function [flag,nSteps] = simulate( ...
        b0,s0,k,d0,o,g,Fvec,dt,Tper,N_all, ...
        burnin,chk,tolPer,maxSteps,pulse,k_check)
aCycles  = 10;           

Jpat = sparse([1 1 0;
               1 1 1;
               0 1 1]);

nPer  = numel(Fvec);
Fstar = (Fvec.^k) ./ mean(Fvec.^k);
Tmax  = maxSteps * dt;

warmSteps = (burnin + chk)*nPer;
tspan1    = 0:dt:warmSteps*dt;
X0 = [0; 0; 0]; 
opts = odeset('RelTol',1e-5,...  
              'AbsTol',1e-5,...
              'MaxStep',dt);
opts.JPattern = Jpat;

[tWarm,YWarm] = ode15s(@rhs_lin,tspan1,X0,opts);
state0 = YWarm(end,:)';   

bufSize = chk;
meanBuf = zeros(4, bufSize);                
for i = 1:chk
    idx = size(YWarm,1) - chk*nPer + (i-1)*nPer + (1:nPer);
    meanEIR = mean(YWarm(idx,:),1)';         % 3×1
    Smean   = N_all - sum(meanEIR);          % 1×1
    meanBuf(:,i) = [Smean; meanEIR];        
end
prevMean = meanBuf(:,end);              
bufPtr   = 1;                  


cyclesHit   = 0;
prevMean    = meanBuf(:,end);
cyclesDone  = burnin;          
stepGlobal  = warmSteps;       


while stepGlobal < maxSteps
    blockCycles = min(aCycles, (maxSteps-stepGlobal)/nPer);
    if blockCycles < 1, break; end
    stepsBlk = blockCycles * nPer;
    tspanBlk = (stepGlobal+1)*dt : dt : (stepGlobal+stepsBlk)*dt;

    recentMean = mean(meanBuf,2);           % 4×1
    AbsTolVec3  = max(1e-6, 1e-5*recentMean(2:4));
    opts.AbsTol = AbsTolVec3;
    opts.RelTol = 1e-6;

    blkStart = tic;                       
    [~,Yblk]  = ode15s(@rhs_lin,tspanBlk,state0,opts);
    state0 = Yblk(end,:)';
    blkTime  = toc(blkStart);             
    for c = 1:blockCycles
        cyclesDone = cyclesDone + 1;
        rng = (c-1)*nPer + (1:nPer);
        if size(Yblk,1) < rng(end)
        flag='odefail'; periodLen=NaN; nSteps=stepGlobal+size(Yblk,1); return;
        end
        meanWeekEIR = mean(Yblk(rng,:),1)';
        Smean       = N_all - sum(meanWeekEIR);  
        meanWeek    = [Smean; meanWeekEIR];     

        % fprintf('[cycle %3d]  <S,E,I,R>/N = %.6g %.6g %.6g %.6g\n', ...
        %     cyclesDone, meanWeek(1)/N_all, meanWeek(2)/N_all, ...
        %     meanWeek(3)/N_all, meanWeek(4)/N_all);

        if cyclesDone > burnin
            tolVec   = max(tolPer, meanWeek .* 10.^(-k_check));
            meanDiff = abs(meanWeek - prevMean);

            if nnz(meanDiff < tolVec) >= 3     
                cyclesHit = cyclesHit + 1;
            else
                cyclesHit = 0;
            end
            prevMean = meanWeek;

            if cyclesHit >= chk
                flag      = 'periodic';
                nSteps    = stepGlobal + c*nPer;
                return;
            end
        end

        meanBuf(:,bufPtr) = meanWeek;
        bufPtr = mod(bufPtr,bufSize) + 1;
    end

    % fprintf('[%s]  %2d cycles (%5d steps)  %.3f s  ⇢  %.2e %.2e %.2e\n', ...
    %     datestr(now,'HH:MM:SS'), blockCycles, stepsBlk, blkTime, AbsTolVec3);

    state0      = Yblk(end,:)';
    stepGlobal  = stepGlobal + stepsBlk;
end

% X0        = [0;0;0];
epsEI     = 1e-15 * N_all;
X0_log    = log([1e-12;1e-12;1e-12] + epsEI);

% [tWarm,YWarm] = ode15s(@rhs,tspan1,X0,opts);
[tWarm,UWarm] = ode15s(@rhs_log,tspan1,X0_log,opts);
YWarm = exp(UWarm) - epsEI;
state0_log = UWarm(end,:)';       

bufSize = chk;
meanBuf = zeros(4, bufSize);                
for i = 1:chk
    idx = size(YWarm,1) - chk*nPer + (i-1)*nPer + (1:nPer);
    meanEIR = mean(YWarm(idx,:),1)';         % 3×1
    Smean   = N_all - sum(meanEIR);          % 1×1
    meanBuf(:,i) = [Smean; meanEIR];        
end
prevMean = meanBuf(:,end);                
bufPtr   = 1;                  

cyclesHit   = 0;
prevMean    = meanBuf(:,end);
cyclesDone  = burnin;           
stepGlobal  = warmSteps;        

hitCnt     = zeros(4,1);   
cyclesDone = burnin;         
stepGlobal = warmSteps;

while stepGlobal < maxSteps

    blockCycles = min(aCycles, (maxSteps-stepGlobal)/nPer);
    if blockCycles < 1, break; end
    stepsBlk = blockCycles * nPer;
    tspanBlk = (stepGlobal+1)*dt : dt : (stepGlobal+stepsBlk)*dt;

    opts.AbsTol = 1e-5;
    opts.RelTol = 1e-5;

    blkStart = tic;                   
    [~,Ublk]  = ode15s(@rhs_log,tspanBlk,state0_log,opts);
    Yblk      = exp(Ublk) - epsEI;
    state0_log = Ublk(end,:)';
    blkTime  = toc(blkStart);       

    for c = 1:blockCycles
        cyclesDone = cyclesDone + 1;
        rng = (c-1)*nPer + (1:nPer);
        if size(Yblk,1) < rng(end)
        flag='odefail'; periodLen=NaN; nSteps=stepGlobal+size(Yblk,1); return;
        end
        meanWeekEIR = mean(Yblk(rng,:),1)';
        Smean       = N_all - sum(meanWeekEIR);  
        meanWeek    = [Smean; meanWeekEIR];       


        % fprintf('[cycle %3d]  <S,E,I,R>/N = %.6g %.6g %.6g %.6g\n', ...
        %     cyclesDone, meanWeek(1)/N_all, meanWeek(2)/N_all, ...
        %     meanWeek(3)/N_all, meanWeek(4)/N_all);


        if cyclesDone > burnin
            tolVec   = max(tolPer, meanWeek .* 10.^(-k_check));
            meanDiff = abs(meanWeek - prevMean);

            if nnz(meanDiff < tolVec) >= 3     
                cyclesHit = cyclesHit + 1;
            else
                cyclesHit = 0;
            end
            prevMean = meanWeek;

            if cyclesHit >= chk
                flag      = 'periodic-after-further-check';
                periodLen = Tper;
                nSteps    = stepGlobal + c*nPer;
                return;
            end
        end
        

        meanBuf(:,bufPtr) = meanWeek;
        bufPtr = mod(bufPtr,bufSize) + 1;
    end

    % fprintf('[%s]  %2d cycles (%5d steps)  %.3f s  ⇢  %.2e %.2e %.2e\n', ...
    %     datestr(now,'HH:MM:SS'), blockCycles, stepsBlk, blkTime, AbsTolVec3);

    state0      = Yblk(end,:)';
    stepGlobal  = stepGlobal + stepsBlk;
end


flag      = 'unsettled';
nSteps    = stepGlobal;








function dX = rhs_lin(tHr,X)
    eps=1e-8;
    E = max(X(1),eps);  I = max(X(2),eps);  R = max(X(3),eps);
        % E = max(X(1),0);  I = max(X(2),0);  R = max(X(3),0);
        % E = X(1);  I = X(2);  R = X(3);
    S = max(N_all - (E + I + R), eps);

    idx  = mod(round(tHr/dt), nPer) + 1;
    Ffac = Fstar(idx);
    phi  = sum(pulse.A .* exp(-0.5*((tHr-pulse.tc)./pulse.tau).^2));
    ph   = phi * Ffac;
    dE =  b0*Ffac*S*I + (1-pulse.theta)*ph*S - s0*Ffac*E - d0*E;
    dI =  s0*Ffac*E   +   pulse.theta*ph*S   - g*I;
    dR =  g*I + d0*E  - o*R;

    dX = [dE; dI; dR];
end

function dU = rhs_log(tHr,U)
    eps=1e-8;

    E = max(exp(U(1)) + epsEI, eps);
    I = max(exp(U(2)) + epsEI, eps);
    R = max(exp(U(3)) + epsEI, eps);
    S = max(N_all - (E + I + R), eps);

    idx  = mod(round(tHr/dt), nPer) + 1;
    Ffac = Fstar(idx);
    phi  = sum(pulse.A .* exp(-0.5*((tHr-pulse.tc)./pulse.tau).^2));
    ph   = phi * Ffac;


    dE =  b0*Ffac*S*I + (1-pulse.theta)*ph*S - s0*Ffac*E - d0*E;
    dI =  s0*Ffac*E   +   pulse.theta*ph*S   - g*I;
    dR =  g*I + d0*E  - o*R;

 
    dU = [(dE)/(E+epsEI);
          (dI)/(I+epsEI);
          (dR)/(R+epsEI)];
end

end
