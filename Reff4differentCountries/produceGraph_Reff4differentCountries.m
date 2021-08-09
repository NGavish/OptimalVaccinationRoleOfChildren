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
susceptibilityProfile=ones(1,9);susceptibilityProfile(2)=2;
susceptibilityProfileC=ones(1,9);susceptibilityProfileC(1:2)=2;

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

    %% 1.5 1.5 1 ... 1
    C2=diag(susceptibilityProfile)*C;M2=diag(Ni)*C2*diag(1./Ni);[V,d]=eig(M2);C2 =C2/max(d(:));

    A=diag(s10.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10_B(ix)=1./max(d(:));

    A=diag(s20.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20_B(ix)=1./max(d(:));

    %% 1 2 1 ... 1
    C2=diag(susceptibilityProfileC)*C;M2=diag(Ni)*C2*diag(1./Ni);[V,d]=eig(M2);C2 =C2/max(d(:));

    A=diag(s10.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren10_C(ix)=1./max(d(:));

    A=diag(s20.*Ni)*C2*diag(1./Ni);[V,d]=eig(A);
    assortivityOfChildren20_C(ix)=1./max(d(:));

end
defineColors;

t=tiledlayout(3,1, 'TileSpacing','compact','Padding','none')
for ix=1:3
    nexttile(ix);
    hndl(1,ix)=scatter(assortivityOfChildren20,prctOfChildren20,25,'s','filled','MarkerFaceColor',blue,'MarkerEdgeColor',blue);hold on;


    h20{ix}=text(0.03+assortivityOfChildren20,prctOfChildren20,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',blue);
%     h(2).VerticalAlignment='top';%h(9).VerticalAlignment='bottom';
%     h(4).VerticalAlignment='top';h(9).VerticalAlignment='bottom';
    ylabel('Size of age group 0-19')
    ytickformat( 'percentage');
    yyaxis right;
    hndl(2,ix)=scatter(assortivityOfChildren10,prctOfChildren10,25,'o','filled','MarkerFaceColor',orange,'MarkerEdgeColor',orange);hold on;

    h10{ix}=text(0.03+assortivityOfChildren10,prctOfChildren10,countryListNames,'HorizontalAlignment','left','VerticalAlignment','middle','color',orange);
%     h(2).VerticalAlignment='bottom';h(2).HorizontalAlignment='center';h(2).Position=h(2).Position+[0 0.2 0];
%     h(7).VerticalAlignment='top';h(7).HorizontalAlignment='center';h(7).Position=h(7).Position-[0 0.2 0];
%     h(9).VerticalAlignment='bottom';h(9).HorizontalAlignment='center';h(9).Position=h(9).Position+[0 0.2 0];
%     h(8).VerticalAlignment='bottom';h(8).HorizontalAlignment='center';h(8).Position=h(8).Position+[0 0.2 0];
    ylabel('Size of age group 0-9')
    box on;
    xlabel('R_{\rm critical}');
    ytickformat('percentage');

    yyaxis left;hold on;plot([1 1]*mean(assortivityOfChildren20),[-5 155],'b:','linewidth',0.75);ylim([0 60])
    hold on;plot([1 1]*mean(assortivityOfChildren10),[-5 155],'r:','linewidth',0.75);
    ax = gca;
    ax.YAxis(2).Color = orange;ax.YAxis(1).Color = blue;
    ytickformat('percentage');grid on
    set(gca,'xtick',[1:0.25:5]);

    if ix==1
        assortivityOfChildren20=assortivityOfChildren20_B;
        assortivityOfChildren10=assortivityOfChildren10_B;
    else
        assortivityOfChildren20=assortivityOfChildren20_C;
        assortivityOfChildren10=assortivityOfChildren10_C;
    end

    axHndl(ix)=gca;
end
set(gcf,'Position',[323 96 812 701]);
linkaxes(axHndl,'xy');
lg=legend(axHndl(3),[hndl(1,1),hndl(2,1)],'R_{\rm critical}^{20+}','R_{\rm critical}^{10+}','fontsize',10,'Orientation','horizontal');
lg.Location = 'southoutside';



text(axHndl(1),0.03,0.9,'A','units','normalized','fontsize',13)
text(axHndl(2),0.03,0.9,'B','units','normalized','fontsize',13)
text(axHndl(3),0.03,0.9,'C','units','normalized','fontsize',13)
title(axHndl(1),'Children (age group 0-19) are half as susceptible as adults');
title(axHndl(2),'Adolescents (ages 10-19) are equally susceptible as adults');
title(axHndl(3),'Children (ages 0-19) are equally susceptible as adults');

% for ix=[1 7 9]
%      h20{1}(ix).VerticalAlignment='top';h20{1}(ix).HorizontalAlignment='left';h20{1}(ix).Position=h20{1}(ix).Position+[-1e-2 0 0];
% end

% Adjust labels - A
h10{1}(7).HorizontalAlignment='center';h10{1}(7).VerticalAlignment='bottom';h10{1}(7).Position=h10{1}(7).Position+[0 0.5 0];
h20{1}(8).HorizontalAlignment='center';h20{1}(8).VerticalAlignment='bottom';
h20{1}(9).HorizontalAlignment='center';h20{1}(9).VerticalAlignment='top';

% Adjust labels - B
h10{2}(7).HorizontalAlignment='center';h10{2}(7).VerticalAlignment='top';%h10{2}(7).Position=h10{1}(7).Position+[0 -0.5 0];
h10{2}(9).HorizontalAlignment='center';h10{2}(9).VerticalAlignment='top';%h10{2}(7).Position=h10{1}(7).Position+[0 -0.5 0];
h20{2}(8).HorizontalAlignment='center';h20{2}(8).VerticalAlignment='bottom';h20{2}(8).Position=h20{2}(8).Position+[0 0.6 0];
ix=4;h20{2}(ix).HorizontalAlignment='left';h20{2}(ix).VerticalAlignment='top';h20{2}(ix).Position=h20{2}(ix).Position+[-0.01 0 0];
ix=9;h20{2}(ix).HorizontalAlignment='center';h20{2}(ix).VerticalAlignment='top';h20{2}(ix).Position=h20{2}(ix).Position+[-0.04 0 0];

% Adjust labels - C
kx=3;
ix=1;h10{kx}(ix).HorizontalAlignment='left';h10{kx}(ix).VerticalAlignment='top';
ix=4;h20{kx}(ix).HorizontalAlignment='left';h20{kx}(ix).VerticalAlignment='top';h20{kx}(ix).Position=h20{kx}(ix).Position+[-0.01 0 0];
ix=8;h20{kx}(ix).HorizontalAlignment='left';h20{kx}(ix).VerticalAlignment='bottom';h20{kx}(ix).Position=h20{kx}(ix).Position+[-0.01 0 0];
ix=9;h20{kx}(ix).HorizontalAlignment='center';h20{kx}(ix).VerticalAlignment='top';h20{kx}(ix).Position=h20{kx}(ix).Position+[-0.04 0 0];

%h20{1}(9).HorizontalAlignment='center';h20{1}(9).VerticalAlignment='top';

h10{2}
shg
printGraph('../graphs/Rthreshold_variousCountries');
return

