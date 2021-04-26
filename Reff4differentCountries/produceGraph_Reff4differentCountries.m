function produceGraph_Reff4differentCountries
close all

%% Parameters
betaVac=0.2; effVac=0.75;
recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;

%% Select Country
countryIdx=1;
countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"};
countryListNames={"Belgium", "USA", "India", "Spain", "Zimbabwe", "Brazil", "China", "South Africa", "Poland"};
subplot(3,3,[1 2 4 5])
%% Run loop
s10=ones(1,9)*betaVac;s10(1)=1;
s20=s10;s20(2)=1;
susceptibilityProfile=ones(1,9);susceptibilityProfile(1:2)=1.5;

M=numel(countryList);
for ix=1:M
    %% Load country data
    country=countryList{ix};
    countryData=load(join(['../countryData/',country,'_data.mat'],''));
    C=countryData.contactMatrix;
    Ni=countryData.agDist';
    A=diag(s10.*Ni)*C*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10(ix)=1./max(d(:));
    
    A=diag(s20.*Ni)*C*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20(ix)=1./max(d(:));
    
    prctOfChildren20(ix)=100*sum(countryData.agDist(1:2));
    prctOfChildren10(ix)=100*sum(countryData.agDist(1));
    
    C2=diag(susceptibilityProfile)*C;M2=diag(Ni)*C2*diag(1./Ni);[V,d]=eig(M2);C2 =C2/max(d(:));
    
    A=diag(s10.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10_B(ix)=1./max(d(:));
    
    A=diag(s20.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20_B(ix)=1./max(d(:));
    
    
end
defineColors;

t=tiledlayout(2,1, 'TileSpacing','compact','Padding','none')
for ix=1:2
    nexttile(ix);
    hndl(1,ix)=scatter(assortivityOfChildren20,prctOfChildren20,25,'s','filled','MarkerFaceColor',blue,'MarkerEdgeColor',blue);hold on;
    
    
    h=text(0.03+assortivityOfChildren20,prctOfChildren20,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',blue);
    h(2).VerticalAlignment='top';%h(9).VerticalAlignment='bottom';
    h(4).VerticalAlignment='top';h(9).VerticalAlignment='bottom';
    ylabel('Size of age group 0-19')
    ytickformat( 'percentage');
    yyaxis right;
    hndl(2,ix)=scatter(assortivityOfChildren10,prctOfChildren10,25,'o','filled','MarkerFaceColor',orange,'MarkerEdgeColor',orange);hold on;
    
    h=text(0.03+assortivityOfChildren10,prctOfChildren10,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',orange);
    h(2).VerticalAlignment='bottom';h(2).HorizontalAlignment='center';h(2).Position=h(2).Position+[0 0.2 0];
    h(7).VerticalAlignment='top';h(7).HorizontalAlignment='center';h(7).Position=h(7).Position-[0 0.2 0];
    h(9).VerticalAlignment='bottom';h(9).HorizontalAlignment='center';h(9).Position=h(9).Position+[0 0.2 0];
    h(8).VerticalAlignment='bottom';h(8).HorizontalAlignment='center';h(8).Position=h(8).Position+[0 0.2 0];
    ylabel('Size of age group 0-9')
    box on;
    xlabel('R_{\rm critical}');
    ytickformat('percentage');
    
    yyaxis left;hold on;plot([1 1]*mean(assortivityOfChildren20),[15 55],'b:','linewidth',0.75);
    hold on;plot([1 1]*mean(assortivityOfChildren10),[15 55],'r:','linewidth',0.75);
    ax = gca;
    ax.YAxis(2).Color = orange;ax.YAxis(1).Color = blue;
    ytickformat('percentage');grid on
    set(gca,'xtick',[1:0.25:4]);
    
    
    assortivityOfChildren20=assortivityOfChildren20_B;
    assortivityOfChildren10=assortivityOfChildren10_B;
    
    axHndl(ix)=gca;
end
set(gcf,'Position',[323 243 655 554]);
linkaxes(axHndl,'xy');
lg=legend(axHndl(1),[hndl(1,1),hndl(2,1)],'R_{\rm critical}^{20+}','R_{\rm critical}^{10+}','fontsize',10,'Orientation','horizontal');
lg.Location = 'southoutside';

for ix=[4 7 9]
    h(ix).VerticalAlignment='top';h(ix).HorizontalAlignment='left';h(ix).Position=h(ix).Position+[-1e-2 0 0];
end

text(axHndl(1),0.05,0.9,'A','units','normalized','fontsize',13)
text(axHndl(2),0.05,0.9,'B','units','normalized','fontsize',13)
title(axHndl(1),'Estimated age dependent susceptibly profile');
title(axHndl(2),'Increased susceptibly of age group 0-19');
shg
printGraph('../graphs/Rthreshold_variousCountries');
return

