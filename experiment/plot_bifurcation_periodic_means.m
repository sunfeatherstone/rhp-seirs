
clear; clc;


T = readtable('experiment/bifurcation_scan_periodic_means.csv');

beta0     = T.beta0;          
period_idx= T.period;         
meanI     = T.mean_I;         


uniq_beta = unique(beta0,'stable');        
num_beta  = numel(uniq_beta);

cval = nan(num_beta,1);                      

for k = 1:num_beta
    mask   = beta0 == uniq_beta(k);  
    yvals  = meanI(mask);
    D      = abs(yvals - yvals.');     
    D(D==0)= inf;   
    mNND   = mean(min(D));                   

    cval(k)= log10(mNND + eps);              
end


[~,loc] = ismember(beta0,uniq_beta);         
c       = cval(loc);                         


fig = figure( ...
    'Units','centimeters', ...  
    'Position',[0 0 40 20]);



ptSize = 2;
scatter(beta0, meanI, ptSize, c, 'filled', ...
        'MarkerFaceAlpha',0.85,'MarkerEdgeAlpha',0.85);

set(gca,'XScale','log');


nCol   = 256;
cmGray = flipud(gray(nCol));
skip   = round(nCol * 0.3);      
colormap(cmGray(skip+1:end,:));   


q = prctile(cval,[2 98]);
caxis(q);

hcb = colorbar;
hcb.Label.String  = 'log_{10}(mNND) per \beta_0 stripe';
hcb.TickDirection = 'out';

xlabel('\beta_0 (log scale)');
ylabel('Mean I per period');
title('Bifurcation scan coloured by log_{10} mNND (grayscale)');
grid on;box on;
set(gca,'FontSize',12);
ylim([0-0.001   max(meanI)+0.001]); 

ax = gca;
ti = ax.TightInset;
set(ax,'LooseInset',[ti(1) ti(2) 0.04 ti(4)]);   

betaStar = 0.0686;                     
xline(betaStar,'k--','\beta_0 = 0.0686', ...
      'LabelHorizontalAlignment','left', ...
      'LabelVerticalAlignment','middle');


fig = gcf;                                  
fig.Units          = 'centimeters';


set(fig,'PaperUnits','centimeters');
set(fig,'PaperPosition',fig.Position);     
set(fig,'PaperSize',fig.Position(3:4));    
set(fig,'PaperPositionMode','manual');     

print(fig,'-dpdf','-opengl','-r300','-loose','bif_scan.pdf');






