%====================== target_hsic.m ======================
% 手搓目标导向 HSIC：RhP–SEIRS 阈值 (rho>1) 为目标
rng(0);                           % 可重复


clc; close all; clearvars -except forward floquet_multiplier    % 保留函数
format long g

% ---------- 公共：节律样条 ----------
S  = load("rhythm/preprocessed_data/weibo_spline_pp_98.mat");
pp = mkpp(S.breaks, S.coefs);
F_base = @(tau) ppval(pp, mod(tau, 84));
F_hour = @(t_hr) F_base(t_hr/2);

% ---------- 参数 ----------
paramFiles = { ...
    'best_parameter/Qingdao.mat', ...
    'best_parameter/Kashgar.mat', ...
    'best_parameter/Chengdu.mat' };

nSample = 10000;     % 1e4 足够稳定；如需更快可调 5000
D_feat  = 200;       % RFF 维度

% ---------- 主循环 ----------
for f = 1:numel(paramFiles)
    data = load(paramFiles{f}, 'theta_refined');
    [ tbl, raw ]  = run_target_hsic(data.theta_refined, F_hour, nSample, D_feat);
    
    fprintf('\n=== %s : Target-oriented HSIC (ρ>1) ===\n', ...
            erase(paramFiles{f}, 'best_parameter/'));
    for r = 1:size(tbl,1)
        fprintf('%8s :  %.4f\n', tbl{r,1}, tbl{r,2});
    end
    score_stack(f,:) = raw; 

    theta_hat=data.theta_refined;

%========== (r_beta, r_gamma) 等高线 ρ≈1 ==========
beta_hat  = theta_hat(1);
gamma_hat = theta_hat(6);

rB = logspace(log10(0.5), log10(2), 101);
rG = logspace(log10(0.5), log10(2), 101);
[RB, RG] = meshgrid(rB, rG);
rho_grid = zeros(size(RB));

for ii = 1:numel(RB)
    th_tmp      = theta_hat;
    th_tmp(1)   = log(RB(ii)) + (beta_hat);
    th_tmp(6)   = log(RG(ii)) + (gamma_hat);
    rho_grid(ii)= floquet_multiplier(th_tmp, F_hour);
end

% ---------- 画图（对数坐标 + 科学计数色条） ----------
fig2 = figure('Visible','off',...
              'Units','centimeters','Position',[2 2 14 12]);   % 略放大
ax  = axes(fig2);

% ---------- 1) 先把 rho 转成对数 -------------
Lgrid = log10(rho_grid);          % ★ 新增：对数网格

% ------------- 伪色图改用 Lgrid --------------
contourf(ax, RB, RG, Lgrid, 24, 'LineColor','none'); hold(ax,'on');
% ------------- 黑线仍然画 rho = 1 -------------
contour (ax, RB, RG, rho_grid, [1 1], 'k', 'LineWidth',2.0);  % 不动

% 3) 当前点 (1,1)
plot    (ax, 1, 1, 'r^', 'MarkerFaceColor','r','MarkerSize',8);

% -------- 坐标轴改成对数并设置 tick --------
% ===== 建议：把 tick 值单独列出来 =====
ticks = [0.5 : 0.1 : 2];          % 0.5, 0.6, 0.7, …, 2.0
set(ax, 'XScale','log', 'YScale','log', ...
        'XLim',[0.5 2], 'YLim',[0.5 2], ...
        'XTick',ticks, 'YTick',ticks, ...
        'FontSize',9,'Box','on', ...
        'GridAlpha',0.25,'LineWidth',0.7);

% （可选）漂亮地打印 tick 标签，防止 MATLAB 自动用指数记号
set(ax,'XTickLabel',arrayfun(@(v)sprintf('%.1f',v),ticks,'uni',0));
set(ax,'YTickLabel',arrayfun(@(v)sprintf('%.1f',v),ticks,'uni',0));

grid(ax,'on');

% （可选）给刻度文字转个角度，防止挤在一起
% xtickangle(ax,45);   % ytickangle 同理

xlabel(ax,'$r_{\beta}=\beta_0/\hat{\beta}_0$',...
           'Interpreter','latex','FontSize',11);
ylabel(ax,'$r_{\gamma}=\gamma/\hat{\gamma}$',...
           'Interpreter','latex','FontSize',11);

[~, city] = fileparts(paramFiles{f});
title(ax, sprintf('%s: $\\rho(r_{\\beta},r_{\\gamma})$', city),...
          'Interpreter','latex');

cb = colorbar(ax);                 % ← 已有
cb.Title.String      = 'lg\rho';         % 把说明文字写在 Title 里
cb.Title.Interpreter = 'tex';            % 保持原来的解释器
cb.Ruler.TickLabelFormat = '%.1f'; % ★ 统一保留 1 位小数
% cb.Label.FontSize    = 9;
% cb.Label.Rotation    = 0;
% cb.Label.VerticalAlignment   = 'bottom';

% ---------- 导出 ----------
pdfName = sprintf('rho_contour_%s.pdf', city);
exportgraphics(fig2, pdfName,...
               'ContentType','vector','BackgroundColor','none');
close(fig2);
end


%% --------- 画柱状图并导出 PDF ---------
paramNames = {'\beta_0','\sigma_0','\kappa','\gamma','\delta'};
cityNames  = {'Qingdao','Kashgar','Chengdu'};

fig = figure('Units','centimeters','Position',[2 2 24 10]);

% 1) 画柱并拿到句柄
bh = bar(paramNames, score_stack.', 'grouped');   % bh(1:3) 是三条系列

% 2) 统一配色：深蓝 / 中蓝 / 浅蓝
set(bh(1),'FaceColor',[0.04 0.28 0.62]);   % 深蓝
set(bh(2),'FaceColor',[0.10 0.45 0.80]);   % 中蓝
set(bh(3),'FaceColor',[0.55 0.75 0.92]);   % 浅蓝

ylabel('Target-oriented HSIC score','FontSize',11)
legend(cityNames,'Location','northeast','FontSize',9)
ylim([0 0.42]); grid on; box on;

exportgraphics(fig,'HSIC_grouped_bar.pdf', ...
               'ContentType','vector','BackgroundColor','none');








function [HSIC_table, raw_scores] = run_target_hsic(theta_refined, F_hour, nSample, D_feat)
% === 目标导向 RFF–HSIC，返回 6×2 表：{name, score} ===
    rng(0);                                % 固定随机性
    names = {'beta0','sigma0','kappa','gamma','delta'};
    idx_use = [1 2 5 6 7];               % 六个可控率
    theta_hat = theta_refined(:);
    p_hat  = theta_hat(idx_use);
    k      = numel(idx_use);

    % --- 对数均匀采样: r = exp( log0.5 + U*(log2-log0.5) )
    S_unit = lhsdesign(nSample,k,'criterion','maximin','iterations',50);
    log_lo = log(0.5);   log_hi = log(2);
    Theta_s = p_hat.' + (log_lo + S_unit*(log_hi-log_lo));  % n×k

    % ---------- 计算 rho(Θ) ----------
    rho_vec = zeros(nSample,1);
    for ii = 1:nSample
        th_tmp          = theta_hat;
        th_tmp(idx_use) = Theta_s(ii,:);
        rho_vec(ii)     = floquet_multiplier(th_tmp, F_hour);
    end
    Z      = double(rho_vec > 1);          % 目标二分类
    Zc     = Z - mean(Z);                  % 居中
    Zc2    = sum(Zc.^2);
    
    % ---------- 低内存 RFF–HSIC ----------
    S_HSIC = zeros(1,k);
    for j = 1:k
        Xj = Theta_s(:,j);
        sig = median(pdist(Xj));
        W   = randn(1,D_feat) / sig;
        b   = 2*pi*rand(1,D_feat);
        Phi = sqrt(2/D_feat) * cos( Xj*W + b );
        Phi = Phi - mean(Phi,1);
        c   = Phi' * Zc;
        hsicXZ = sum(c.^2) / nSample^2;
        M     = Phi' * Phi;
        hsicXX= sum(M(:).^2) / nSample^2;
        hsicZZ= (Zc2^2) / nSample^2;
        S_HSIC(j) = hsicXZ / sqrt(hsicXX*hsicZZ + eps);
    end
    
    % 排序并返回 cell 表
    raw_scores = S_HSIC; 
    [score_sorted, ord] = sort(S_HSIC,'descend');
    HSIC_table = [ names(ord).' , num2cell(score_sorted(:)) ];
end












%====================== floquet_multiplier.m ================
function rho = floquet_multiplier(theta20,F_hour)
% 返主 Floquet 乘子 (2×2 线化块) — 给 HSIC 调用
    % -----  unpack (只用 beta0 sigma0 kappa gamma delta)
    beta0 = exp(theta20(1));
    sigma0= exp(theta20(2));
    kappa = exp(theta20(5));
    gamma = exp(theta20(6));
    delta = exp(theta20(7));

    % 周期 T=168 h; 细分步长同前
    T_tot = 168;   dt = 0.5;  nStep = T_tot/dt;
    t_hr  = (0:nStep-1).' * dt;
    F     = F_hour(t_hr);          % 直接复用节律样条
    Fp    = F.^kappa;  F_star = Fp / mean(Fp);

    % 积分变分方程 y' = J(t)*y   (2×2)
    y = eye(2);
    for kk = 1:nStep
        b = beta0 * F_star(kk);
        s = sigma0* F_star(kk);
        J = [- (s+delta) ,  b ;
               s         , -gamma];
        y = y + dt * J * y;         % Euler–Maruyama 已够用
    end
    rho = max(abs(eig(y)));
end


function score = hsic_rff_1d(X, Z, D)
% 低内存 RFF–HSIC (适用于一维 X；Z 任意向量或 0/1)
% X: n×1, Z: n×1
% D: RFF 维数

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

