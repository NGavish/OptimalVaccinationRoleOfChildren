function main
close all
%% Parameters
fname='dataRecovered'

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
% susceptibilityProfile=ones(1,9);%susceptibilityProfile(1:2)=1.5;
    M=numel(countryList);
    for ix=1:M
        %% Load country data
        country=countryList{ix};
        countryData=load(join(['../countryData/',country,'_data.mat'],''));
        C=countryData.contactMatrix;
        Ni=countryData.agDist';
        
        A=diag(Ni)*C*diag(1./Ni);[V,d]=eig(A);[dummy,idx]=max(max(d));
        r10=0.1*V(:,idx)'/sum(V(:,idx)'.*Ni);
        r20=2*r10;
        
        % None recovered
        s10=ones(1,9)*betaVac;s10(1)=1;
        s20=s10;s20(2)=1;
        
        A=diag(s10.*Ni)*C*diag(1./Ni);[V,d]=eig(A);
        assortivityOfChildren10(ix)=1./max(d(:));
        
        A=diag(s20.*Ni)*C*diag(1./Ni);[V,d]=eig(A);
        assortivityOfChildren20(ix)=1./max(d(:));
        
        prctOfChildren20(ix)=100*sum(countryData.agDist(1:2));
        prctOfChildren10(ix)=100*sum(countryData.agDist(1));
        
        % 10% recovered
        s10=(1-r10)*betaVac;s10(1)=1-r10(1);
        s20=s10;s20(2)=1-r10(2);
        C2=C;
        A=diag(s10.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
        assortivityOfChildren10_B(ix)=1./max(d(:));
        
        A=diag(s20.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
        assortivityOfChildren20_B(ix)=1./max(d(:));
        
        % 20% recovered
        s10=(1-r20)*betaVac;s10(1)=1-r20(1);
        s20=s10;s20(2)=1-r20(2);
        %C2=diag(susceptibilityProfile)*C;M2=diag(Ni)*C2*diag(1./Ni);[V,d]=eig(M2);C2 =C2/max(d(:));
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
    t=tiledlayout(2,1, 'TileSpacing','compact','Padding','none')
    for ix=2:3
    nexttile;
    assortivityOfChildren20=assortivityOfChildren20List{ix}
    assortivityOfChildren10=assortivityOfChildren10List{ix}
    
    plot([1 1]*mean(assortivityOfChildren20),[15 55],'b:','linewidth',1);hold on;
    hline=plot([1 1]*mean(assortivityOfChildren20List{1}),[15 55],'-.','Color',0.7*[1 1 1],'linewidth',1.2);
    hold on;plot([1 1]*mean(assortivityOfChildren10),[15 55],'r:','linewidth',1);
    plot([1 1]*mean(assortivityOfChildren10List{1}),[15 55],'-.','Color',0.7*[1 1 1],'linewidth',1.2);
    hndl(1,ix)=scatter(assortivityOfChildren20,prctOfChildren20,25,'s','filled','MarkerFaceColor',blue,'MarkerEdgeColor',blue);hold on;
        
    h=text(0.03+assortivityOfChildren20,prctOfChildren20,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',blue);

    h(2).VerticalAlignment='top';%h(9).VerticalAlignment='bottom';
    h(4).VerticalAlignment='top';h(9).VerticalAlignment='bottom';
    ylabel('Size of age group 0-19')
    ytickformat( 'percentage');
    yyaxis right;
    hndl(2,ix)=scatter(assortivityOfChildren10,prctOfChildren10,25,'o','filled','MarkerFaceColor',orange,'MarkerEdgeColor',orange);hold on;
ytickformat( 'percentage');
    h=text(0.03+assortivityOfChildren10,prctOfChildren10,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',orange);
        h(1).VerticalAlignment='bottom';h(1).HorizontalAlignment='right';h(1).Position=h(1).Position+[0 0.2 0];

    h(2).VerticalAlignment='bottom';h(2).HorizontalAlignment='center';h(2).Position=h(2).Position+[0 0.2 0];
    h(7).VerticalAlignment='top';h(7).HorizontalAlignment='right';h(7).Position=h(7).Position-[0 0.2 0];
    h(9).VerticalAlignment='bottom';h(9).HorizontalAlignment='center';h(9).Position=h(9).Position+[0 0.2 0];
    h(8).VerticalAlignment='bottom';h(8).HorizontalAlignment='center';h(8).Position=h(8).Position+[0 0.2 0];
    ylabel('Size of age group 0-9')
    box on;%legend('R_{\rm critical}^{20+}','R_{\rm critical}^{10+}','location','best');
    xlabel('R_{\rm critical}');
    
    yyaxis left;hold on;
    ax = gca;
    ax.YAxis(2).Color = orange;ax.YAxis(1).Color = blue; 
    ytickformat('percentage');grid on
    xlim([1 5]);set(gca,'xtick',[1:0.25:5]);

    
    axHndl(ix)=gca;
    end
    set(gcf,'Position',[323 243 655 554]);
linkaxes([axHndl(2) axHndl(3)],'xy');
lg=legend(axHndl(2),[hndl(1,3),hndl(2,3) hline],'R_{\rm critical}^{20+}','R_{\rm critical}^{10+}','Reference (none recovered)','fontsize',10,'Orientation','horizontal');
 lg.Location = 'southoutside';

for ix=[4 7 9]
h(ix).VerticalAlignment='top';h(ix).HorizontalAlignment='left';h(ix).Position=h(ix).Position+[-1e-2 0 0];
end
text(axHndl(2),0.05,0.9,'A','units','normalized','fontsize',13)
text(axHndl(3),0.05,0.9,'B','units','normalized','fontsize',13)
title(axHndl(2),'10% recovered or with prior immunity');
title(axHndl(3),'20% recovered or with prior immunity');

shg

      printGraph('../graphs/Rthreshold_variousCountries_recovered');
    return

