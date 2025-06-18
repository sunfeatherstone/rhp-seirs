clear; clc;

T = readtable('experiment/A0_bifurcation_scan_periodic_means.csv');

A0     = T.A0;   
period_idx= T.period;   
meanI     = T.mean_I;  

figure;
sz = 8; 
scatter(A0, meanI, sz, period_idx, 'filled');
set(gca, 'XScale','log');        
colormap(flipud(gray(100)));   
caxis([1 100]);                
hcb = colorbar;               
hcb.Label.String = 'Period index';
xlabel('\A_0 (log scale)');
ylabel('Mean I per period');
title('Bifurcation Scan: Periodic Means of I');
grid on;
set(gca,'FontSize',12);