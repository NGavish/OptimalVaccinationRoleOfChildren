function main
% This function prepares SI Figure "Effect of vaccine efficacy on critical reproduction numbers."

close all
%% Parameters
fname='dataEfficacy'

% Define vaccine efficacy scenarios
betaVac=0.1; effVac=0.5;
betaVacB=0.05;betaVacC=0.15;
recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;

%% Select Country
countryIdx=1;
countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"};
countryListNames={"Belgium", "USA", "India", "Spain", "Zimbabwe", "Brazil", "China", "South Africa", "Poland"};
subplot(3,3,[1 2 4 5])
%% Run loop
s10=ones(1,9)*betaVac;s10(1)=1;
s20=s10;s20(2)=1;

% For each country, compute the assortivity of children (Rcritical) for
% different vaccine efficacy levels
M=numel(countryList);
for ix=1:M
    %% Load country data
    country=countryList{ix};
    countryData=load(join(['../countryData/',country,'_data.mat'],''));
    C=countryData.contactMatrix;
    Ni=countryData.agDist';

    s10=ones(1,9)*betaVac;s10(1)=1;
    s20=s10;s20(2)=1;

    A=diag(s10.*Ni)*C*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10(ix)=1./max(d(:));

    A=diag(s20.*Ni)*C*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20(ix)=1./max(d(:));

    prctOfChildren20(ix)=100*sum(countryData.agDist(1:2));
    prctOfChildren10(ix)=100*sum(countryData.agDist(1));

    %% 95% vaccine efficacy
    s10=ones(1,9)*betaVacB;s10(1)=1;
    s20=s10;s20(2)=1;
    C2=C;
    A=diag(s10.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10_B(ix)=1./max(d(:));

    A=diag(s20.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20_B(ix)=1./max(d(:));

    %% 85% vaccine efficacy
    s10=ones(1,9)*betaVacC;s10(1)=1;
    s20=s10;s20(2)=1;
    C2=C;
    A=diag(s10.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10_C(ix)=1./max(d(:));

    A=diag(s20.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20_C(ix)=1./max(d(:));

end
defineColors;
assortivityOfChildren20List{1}=assortivityOfChildren20;
assortivityOfChildren20List{2}=assortivityOfChildren20_B;
assortivityOfChildren20List{3}=assortivityOfChildren20_C;

assortivityOfChildren10List{1}=assortivityOfChildren10;
assortivityOfChildren10List{2}=assortivityOfChildren10_B;
assortivityOfChildren10List{3}=assortivityOfChildren10_C;

% Present data
t=tiledlayout(2,1, 'TileSpacing','compact','Padding','none')
for ix=2:3
    nexttile;
    assortivityOfChildren20=assortivityOfChildren20List{ix}
    assortivityOfChildren10=assortivityOfChildren10List{ix}

    % Plot R^20_critical data
    plot([1 1]*mean(assortivityOfChildren20),[15 55],'b:','linewidth',1);hold on;
    hline=plot([1 1]*mean(assortivityOfChildren20List{1}),[15 55],'-.','Color',0.7*[1 1 1],'linewidth',1.2);
    hold on;plot([1 1]*mean(assortivityOfChildren10),[15 55],'r:','linewidth',1);
    plot([1 1]*mean(assortivityOfChildren10List{1}),[15 55],'-.','Color',0.7*[1 1 1],'linewidth',1.2);
    hndl(1,ix)=scatter(assortivityOfChildren20,prctOfChildren20,25,'s','filled','MarkerFaceColor',blue,'MarkerEdgeColor',blue);hold on;


    h=text(0.03+assortivityOfChildren20,prctOfChildren20,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',blue);



    ylabel('Size of age group 0-19')
    ytickformat( 'percentage');
    yyaxis right;

    % Adjust text labels to avoid overlaps
    h(2).VerticalAlignment='top';
    h(4).VerticalAlignment='top';h(9).VerticalAlignment='bottom';

    % Plot R^10_critical data
    hndl(2,ix)=scatter(assortivityOfChildren10,prctOfChildren10,25,'o','filled','MarkerFaceColor',orange,'MarkerEdgeColor',orange);hold on;

    h=text(0.03+assortivityOfChildren10,prctOfChildren10,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',orange);
  
    ylabel('Size of age group 0-9')
    box on;
    xlabel('R_{\rm critical}');
    ytickformat('percentage');
    yyaxis left;hold on;

    % Adjust text labels to avoid overlaps
    h(1).VerticalAlignment='bottom';h(1).HorizontalAlignment='right';h(1).Position=h(1).Position+[0 0.2 0];
    h(2).VerticalAlignment='bottom';h(2).HorizontalAlignment='center';h(2).Position=h(2).Position+[0 0.2 0];
    h(7).VerticalAlignment='middle';h(7).HorizontalAlignment='right';h(7).Position=h(7).Position+[-0.05 0 0];
    h(9).VerticalAlignment='middle';h(9).HorizontalAlignment='right';h(9).Position=h(9).Position+[-0.05 0.2 0];
    h(8).VerticalAlignment='bottom';h(8).HorizontalAlignment='right';h(8).Position=h(8).Position+[0 0.2 0];

    ax = gca;
    ax.YAxis(2).Color = orange;ax.YAxis(1).Color = blue;
    ytickformat('percentage');grid on
    set(gca,'xtick',[1:0.25:4]);

    axHndl(ix)=gca;
end

% Add title and subplot reference text
set(gcf,'Position',[323 243 655 554]);
linkaxes([axHndl(2) axHndl(3)],'xy');
lg=legend(axHndl(2),[hndl(1,3),hndl(2,3) hline],'R_{\rm critical}^{20+}','R_{\rm critical}^{10+}','90% vaccine efficacy','fontsize',10,'Orientation','horizontal');
lg.Location = 'southoutside';

for ix=[4 7 9]
    h(ix).VerticalAlignment='top';h(ix).HorizontalAlignment='left';h(ix).Position=h(ix).Position+[-1e-2 0 0];
end
text(axHndl(2),0.05,0.9,'A','units','normalized','fontsize',13)
text(axHndl(3),0.05,0.9,'B','units','normalized','fontsize',13)
title(axHndl(2),'95% vaccine efficacy');
title(axHndl(3),'85% vaccine efficacy');

shg

printGraph('../graphs/Rthreshold_variousCountries_vaccineEfficacy');
return
