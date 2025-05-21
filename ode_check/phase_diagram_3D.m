% ===================  R0 / rho  +  (E/N,I/N) scatter  ====================
% 2025-05-20  –  Yushi Sun
% -------------------------------------------------------------------------

clear; close all; clc;
S  = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp = mkpp(S.breaks, S.coefs);
F_base = @(tau) ppval(pp, mod(tau, 84));
F_hour = @(t_hr) F_base(t_hr/2);

paramFiles = {
    'best_parameter/Qingdao.mat'
    'best_parameter/Kashgar.mat'
    'best_parameter/Chengdu.mat' };
cityNames = {'Qingdao','Kashgar','Chengdu'};
colRGB    = [ 0.04 0.28 0.62 ;   % 
              0.10 0.45 0.80 ;   % 
              0.55 0.75 0.92 ];  % 


nCase   = numel(paramFiles);
R0_vec  = zeros(1,nCase);
rho_vec = zeros(1,nCase);
allE = cell(1,nCase);
allI = cell(1,nCase);
allR = cell(1,nCase);
Tweek   = 168;      % h
dt_h    = 0.5;      % h
Tsteps  = 100000;    % 15k ×0.5h ≈ 312.5h


fwdFcn = {@forward_qingdao , @forward_kashgar , @forward_chengdu};

for c = 1:nCase
    data = load(paramFiles{c}, 'theta_refined');
    th   = data.theta_refined(:);           

    % ---- 1) R0, rho ---------------------------------------------------
    b0   = exp(th(1));  s0 = exp(th(2));
    g    = exp(th(6));  d  = exp(th(7));
    R0_vec(c) = b0 * s0 / ( g * (s0+d) );   

    rho_vec(c)= floquet_multiplier(th,F_hour); 


[E,I,R,Npop] = fwdFcn{c}(th,Tsteps,dt_h);  
allE{c} = E(end-960:end)./Npop;
allI{c} = I(end-960:end)./Npop;
allR{c} = R(end-960:end)./Npop;
end


fig1 = figure('Units','centimeters','Position',[2 2 10 7]);
bh   = bar([R0_vec.' rho_vec.'],'grouped');
set(bh(1),'FaceColor',[.6 .6 .6]);
set(bh(2),'FaceColor',[.2 .2 .8]);
set(gca,'XTickLabel',cityNames,'FontSize',9);
ylabel('value'); legend({'$\mathcal R_0$','$\rho(F)$'},'Interpreter','latex');
title('Basic reproduction index vs. Floquet multiplier');
exportgraphics(fig1,'R0_rho_bar.pdf','ContentType','vector','BackgroundColor','none');

% for c = 1:nCase
%     fig = figure('Units','centimeters', ...
%                  'Position',[2 2 11 11]);        
%     scatter(allE{c}, allI{c}, 6, 'filled', ...
%             'MarkerFaceColor', colRGB(c,:), ...
%             'MarkerEdgeColor', 'none');
%     xlabel('E / N');  ylabel('I / N');
%     axis equal;  grid on;  box on;
%     title(['Normalised limit cycle — ' cityNames{c}]);
    
%     pdfName = ['EI_limit_cycle_' cityNames{c} '.pdf'];
%     exportgraphics(fig, pdfName, ...
%                    'ContentType','vector', 'BackgroundColor','none');
%     close(fig);
% end

%% ---------- 图 2-bis : 3-D normalised limit cycles----------
for c = 1:nCase
    fig3 = figure('Units','centimeters','Position',[2 2 12 10]);  %

plot3(allE{c}, allI{c}, allR{c}, '-', ...
      'Color', [0.10 0.45 0.80], 'LineWidth', 1.2);

xlabel('E / N'); ylabel('I / N'); zlabel('R / N');
grid on; box on; axis equal tight; view(45,25);


fixedH = 10;                      
ax = gca;
range = [diff(ax.XLim) diff(ax.YLim) diff(ax.ZLim)];
aspect = max(range) / min(range);
newW   = fixedH * aspect;
set(gcf,'Units','centimeters',...
        'Position',[2 2 newW fixedH]);   


ax.XRuler.Exponent = 0;  ax.YRuler.Exponent = 0; ax.ZRuler.Exponent = 0;
ax.XTickLabel = arrayfun(@(v) sprintf('%.3f',v), ax.XTick,'uni',0);
ax.YTickLabel = arrayfun(@(v) sprintf('%.3f',v), ax.YTick,'uni',0);
ax.ZTickLabel = arrayfun(@(v) sprintf('%.3f',v), ax.ZTick,'uni',0);

title(['Normalised limit cycle — ' cityNames{c}]);

pdfName = ['EIR_limit_cycle_' cityNames{c} '.pdf'];
exportgraphics(gcf, pdfName, ...
               'ContentType','vector','BackgroundColor','none');


end

function [E,I,R,Npop] = integrate_RhP(theta, pulseSpec, Tsteps, dt_h, F_hour)
    % pulseSpec = struct('h_loc',[],'d_idx',[],'tau',[],'tau1',[]);
    % theta 

    S  = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp = mkpp(S.breaks, S.coefs);
F_base = @(tau) ppval(pp, mod(tau, 84));
F_hour = @(t_hr) F_base(t_hr/2);


    b0 = exp(theta(1));   s0 = exp(theta(2));
    Npop = exp(theta(3));
    thetaV = 1/(1+exp(-theta(4)));
    kappa  = exp(theta(5));
    gamma  = exp(theta(6));
    delta  = exp(theta(7));
    omega  = exp(theta(8));


    t_hr = (0:Tsteps-1)' * dt_h;
    F    = F_hour(mod(t_hr,168));          % raw rhythm
    F    = F / mean(F);                    % mean-1
    Fk   = F .^ kappa;                     % raise
    Fk   = Fk / mean(Fk);                  

    % ---------- φ(t) ---------------
    tau  = pulseSpec.tau;
    tau1 = pulseSpec.tau1;
    h_loc = pulseSpec.h_loc;
    d_idx = pulseSpec.d_idx;         
    A_vec = exp(theta(10:10+numel(h_loc)-1));

    t_c = h_loc + (d_idx-1)*24 - 4;     % centre
    phi = zeros(Tsteps,1);
    for i = 1:numel(h_loc)
        if i==1 && ~isempty(tau1)
        w = tau1;
    else
        w = tau;
    end
        phi = phi + A_vec(i)*exp(-0.5*((t_hr-t_c(i))/w).^2);
    end

    S = zeros(Tsteps,1); E = S; I = S; R = S;
    S(1)=Npop;  
    for k=1:Tsteps-1
        b  = b0 * Fk(k);
        s  = s0 * Fk(k);
        ph = phi(k) * Fk(k);

        dS = -b*S(k)*I(k)/Npop + omega*R(k) - ph*S(k);
        dE =  b*S(k)*I(k)/Npop + (1-thetaV)*ph*S(k) - s*E(k) - delta*E(k);
        dI =  s*E(k) + thetaV*ph*S(k) - gamma*I(k);
        dR =  gamma*I(k) + delta*E(k) - omega*R(k);

        S(k+1)=S(k)+dt_h*dS;
        E(k+1)=E(k)+dt_h*dE;
        I(k+1)=I(k)+dt_h*dI;
        R(k+1)=R(k)+dt_h*dR;
    end
end

% ========= 1) Qingdao ====================================================
function [E,I,R,Npop] = forward_qingdao(theta, Tsteps, dt_h)
    pulseSpec.h_loc = [ 8.25  13.25  9  9  11.5  9  8  9  15.5  9  9 ];
    pulseSpec.d_idx = 1:11;
    pulseSpec.tau   = exp(theta(9));
    pulseSpec.tau1  = [];            
    [E,I,R,Npop] = integrate_RhP(theta, pulseSpec, Tsteps, dt_h);
end

% ========= 2) Kashgar ====================================================
function [E,I,R,Npop] = forward_kashgar(theta, Tsteps, dt_h)
    timeShift = theta(17);           
    pulseSpec.h_loc = [20.5+timeShift  9 14 20 9 9];
    pulseSpec.d_idx = [1 2 2 2 3 4];
    pulseSpec.tau   = exp(theta(9));
    pulseSpec.tau1  = exp(theta(16));
    [E,I,R,Npop] = integrate_RhP(theta, pulseSpec, Tsteps, dt_h);
end

% ========= 3) Chengdu ====================================================
function [E,I,R,Npop] = forward_chengdu(theta, Tsteps, dt_h)
    timeShift = theta(18);
    pulseSpec.h_loc = [12+timeShift 11 10 21.25 10 12 10.5];
    pulseSpec.d_idx = [1 2 3 3 4 8 10];
    pulseSpec.tau   = exp(theta(9));
    pulseSpec.tau1  = exp(theta(17));
    [E,I,R,Npop] = integrate_RhP(theta, pulseSpec, Tsteps, dt_h);
end


%====================== floquet_multiplier.m ================
function rho = floquet_multiplier(theta20,F_hour)
    % -----  unpack ( beta0 sigma0 kappa gamma delta)
    beta0 = exp(theta20(1));
    sigma0= exp(theta20(2));
    kappa = exp(theta20(5));
    gamma = exp(theta20(6));
    delta = exp(theta20(7));

    T_tot = 168;   dt = 0.5;  nStep = T_tot/dt;
    t_hr  = (0:nStep-1).' * dt;
    F     = F_hour(t_hr);         
    Fp    = F.^kappa;  F_star = Fp / mean(Fp);

    %y' = J(t)*y   (2×2)
    y = eye(2);
    for kk = 1:nStep
        b = beta0 * F_star(kk);
        s = sigma0* F_star(kk);
        J = [- (s+delta) ,  b ;
               s         , -gamma];
        y = y + dt * J * y;         % Euler–Maruyama 
    end
    rho = max(abs(eig(y)));
end


function score = hsic_rff_1d(X, Z, D)

% X: n×1, Z: n×1
% D: RFF

n  = numel(X);
Zc = Z - mean(Z);
sig = median(pdist(X));
W   = randn(1,D) / sig;
b   = 2*pi*rand(1,D);
Phi = sqrt(2/D) * cos( X*W + b );
Phi = Phi - mean(Phi,1);

c      = Phi' * Zc;                 % D×1
hsicXZ = sum(c.^2)       / n^2;
M      = Phi' * Phi;               % D×D
hsicXX = sum(M(:).^2)    / n^2;
hsicZZ = (Zc'*Zc)^2      / n^2;

score  = hsicXZ / sqrt(hsicXX*hsicZZ + eps);
end