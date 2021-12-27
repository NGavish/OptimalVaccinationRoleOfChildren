function produceAllGraphs(collectData)
addpath('../Codes/AllCountries/AuxilaryFunctions/');
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
defineColors
global R0
if collectData
    display('Computing data for Pareto front (this is not a quick computation)')
    runComputationOfParetoFronts;
end

close all;
defineColors;

Graphs=[
    true, % all-or-none
    true, % hesistancy
    true, % efficacy
    true, % coverage
    true, % comparison to 20+
    true] % recovered
Graphs(1:6)=false;Graphs(1)=true;
R0=6;produceParetoGraph(Graphs,R0) % Produce Figure 6 & all pareto front graphs for  SI 
Graphs(1:6)=false;Graphs(5)=true;
R0=4;produceParetoGraph(Graphs,R0) % Produce Figure 5
R0=8;produceParetoGraph(Graphs,R0) % Produce Figure 7

function produceParetoGraph(whichGraphsToProduce,R0)
defineColors
ix=1;
%% Default scenario
susceptibilityFactor=1;maxPrct=100;VcPrct=73.2;betaVac=0.1;nuVac=1;effVac=1-0.05/betaVac;
presentationOptions.presentUniform=false;

%% All-or-none
if whichGraphsToProduce(ix)

    % Produce SI graph - Pareto front with all-or-none vaccine
    betaVacAllorNone=0;nuVacAllorNone=0.9;

    t=tiledlayout(3,1, 'TileSpacing','compact','Padding','none')

    % Pareto front graph
    axVec{1}=nexttile(1);

    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'All')
    set(h,'Color',blue,'LineStyle','--');g=gca;xlimSave=g.XLim;

    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'above10')
    set(h,'Color',red,'LineStyle','-.');

    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'above20')

    xlim(xlimSave);ylim([-0.25 3.5])

    g=gca;
    lg=legend([g.Children(11) g.Children(8) g.Children(10) g.Children(9) g.Children(5)],'Random allocations - all ages','Random allocations - ages 20 and older','Pareto front - all ages','Pareto front - ages 10 and older','Pareto front - ages 20 and older','location','southwest','autoupdate','off','fontsize',12);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));

    % Allocations along Pareto front for ages 20+
    axVec{2}=axes('Position',[0.40 0.400952380952381 0.59 0.2]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'above20')
    title('Vaccine allocation restricted to ages 20 and above');box on
    grid on
    g2=gca;
    delete(g2.Children(10));delete(g2.Children(10));

    % Allocations along Pareto front for all ages
    axVec{3}=axes('Position',[0.08 0.1 0.87 0.2]);

    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac)
    g3=gca;
    delete(g3.Children(10));delete(g3.Children(10));

    grid on
    title('No age restriction on vaccine allocation')

    set(gcf,'Position',[520 234 610 563]);
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    text(axVec{2},148,20,Ages{3},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},148,70,Ages{4},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{2},148,110,Ages{5},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{2},148,145,Ages{6},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},154,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},154,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},154,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');

    text(axVec{3},136,30,Ages{2},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},136,70,Ages{3},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},136,110,Ages{4},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{3},136,150,Ages{5},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{3},155,90,Ages{6},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},155,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},155,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');

    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    lg2 = legend(axVec{2},Ages,'fontsize',10,'Orientation','vertical','location','westoutside');%legend(nexttile(2), [line1,line2,line3,line4]);
    
    set(gcf,'Position',[520 199 736 598]);

    printGraph('../graphs/ParetoFront_AllorNone');
    close all;
end
ix=ix+1;

%% Hesistancy
if whichGraphsToProduce(ix)
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')

    % Pareto front
    axVec{1}=nexttile([1,2]);

    % Default
    defineColors;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');

    % 10% Hesistancy
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,90,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',red,'LineStyle','-');

    % 20% Hesistancy
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,80,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',blue,'LineStyle','-.');

    xlim([163 197]);

    lg=legend('No vaccine hesistancy','Maximum 90% vaccination in age group','Maximum 80% vaccination in age group','location','best');
    %     set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));

    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,90,VcPrct,betaVac,nuVac,effVac)
    title('Maximum 90% vaccination in age group');

    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,80,VcPrct,betaVac,nuVac,effVac)
    title('Maximum 80% vaccination in age group');

    % Link axes

    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    lg = legend(Ages,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    %     lg.Location = 'southoutside';\
    lg.Layout.Tile='south'

    set(gcf,'Position',[520 392 633 405]);
    shg;
    printGraph('../graphs/ParetoFront_Hesistancy');
    close all;
end
ix=ix+1;

%% Efficacy
if whichGraphsToProduce(ix)
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')

    % Pareto front
    axVec{1}=nexttile([1,2]);

    % Default
    defineColors;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');

    % 95% efficacy
    betaVacVE95=0.05;effVacVE95=1-0.05/betaVacVE95;  % 95% efficacy
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacVE95,nuVac,effVacVE95);hold on;
    set(h,'Color',red,'LineStyle','-');

    % 85% efficacy
    betaVacVE85=0.15;effVacVE85=1-0.05/betaVacVE85;  % 85% efficacy
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacVE85,nuVac,effVacVE85)
    set(h,'Color',blue,'LineStyle','-.');

    xlim([130 227]);ylim([0 4])

    lg1=legend('90% Vaccine Efficacy (as in Figure 5)','95% Vaccine Efficacy','85% Vaccine Efficacy','location','northeast');
    set(lg1.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));


    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,100,VcPrct,betaVacVE95,nuVac,effVacVE95)
    title('95% Vaccine Efficacy');

    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,100,VcPrct,betaVacVE85,nuVac,effVacVE85)
    title('85% Vaccine Efficacy');

    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    lg = legend(Ages,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    %     lg.Location = 'southoutside';\
    lg.Layout.Tile='south'

    set(gcf,'Position',[520 392 676 405]);
    shg
    printGraph('../graphs/ParetoFront_Efficacy');
    close all;
end

ix=ix+1;

%% Vaccination coverage
if whichGraphsToProduce(ix)
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')

    % Pareto front
    axVec{1}=nexttile([1,2]);

    % Default
    defineColors;
    presentationOptions.presentUniform=false;presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');

    % 80% coverage
    VcPrct80=80;VcPrct80=86;VacPrct80=93.2;VacPrct80=100;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct+5,betaVac,nuVac,effVac)
    set(h,'Color',red,'LineStyle','-');

    % 66% coverage
    VcPrct66=66;VcPrct66=59;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-5,betaVac,nuVac,effVac)
    set(h,'Color',blue,'LineStyle','-.');

    xlim([140 215]);ylim([0 4])
    g=gca;
    legend([g.Children(1),g.Children(3),g.Children(2)],'52% vaccine coverage','55% vaccine coverage  (as in Figure 6)','58% vaccine coverage','location','northeast');

    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-5,betaVac,nuVac,effVac)
    title('52% vaccine coverage (68.2% of adult population)');

    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct+5,betaVac,nuVac,effVac)
    title('58% vaccine coverage (78.2% of adult population)');

    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    lg = legend(Ages,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    %     lg.Location = 'southoutside';\
    lg.Layout.Tile='south'
    linkaxes([ axVec{2} axVec{3} ],'y');
    set(gcf,'Position',[520 392 676 405]);
    shg;
    printGraph('../graphs/ParetoFront_coverage');
    close all;
end

ix=ix+1;
defineColors;
%% Comparison to 20+
if whichGraphsToProduce(ix)
    betaVac=0.1;
    t=tiledlayout(3,1, 'TileSpacing','compact','Padding','none')

    % Pareto front
    axVec{1}=nexttile(1);

    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All')
    set(h,'Color',blue,'LineStyle','--');

    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above10')
    set(h,'Color',red,'LineStyle','-.');

    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above20')


    xlim([163.5 204])

    g=gca;
    lg=legend([g.Children(11) g.Children(8) g.Children(10) g.Children(9) g.Children(5)],'Random allocations - all ages','Random allocations - ages 20 and older','Pareto front - all ages','Pareto front - ages 10 and older','Pareto front - ages 20 and older','location','southwest','autoupdate','off','fontsize',12);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));

    axVec{2}=axes('Position',[0.53 0.400952380952381 0.456696428571429 0.2]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above20')
    title('Vaccine allocation restricted to ages 20 and above')
    grid on
    g2=gca;
    delete(g2.Children(10));delete(g2.Children(10));
    axVec{3}=axes('Position',[0.08 0.1 0.89 0.2]);

    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    g3=gca;
    delete(g3.Children(10));delete(g3.Children(10));

    grid on
    title('No age restriction on vaccine allocation')

    xlim(axVec{1},[162 202]);
        Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};

    switch R0
        case 8
            xlim(axVec{1},[180 230])
            ylim(axVec{1},[0.25 4])

            text(axVec{2},215,30,Ages{3},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},215,70,Ages{4},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{2},215,110,Ages{5},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{2},215,145,Ages{6},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},222,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},222,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},222,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');

            text(axVec{3},194.5,20,Ages{1},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},194.5,60,Ages{2},'fontsize',12,'color','w','fontweight','bold');
             text(axVec{3},206,90,Ages{6},'fontsize',12,'color','w','fontweight','bold');
             text(axVec{3},206,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
             text(axVec{3},206,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
             text(axVec{3},206,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');

             set(axVec{3},'Position',[0.32 0.1 0.43 0.2])
            set(axVec{2},'Position',[0.7 0.400952380952381 0.185 0.2])
         case 4
            axis(axVec{1},[95 170 0 3])


            text(axVec{2},130,25,Ages{3},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},130,65,Ages{4},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{2},130,110,Ages{5},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{2},130,145,Ages{6},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},155,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},155,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},155,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');

            text(axVec{3},100,20,Ages{2},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},100,60,Ages{3},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},100,110,Ages{4},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},100,150,Ages{5},'fontsize',12,'color','w','fontweight','bold');

             text(axVec{3},150,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
             text(axVec{3},150,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
             text(axVec{3},150,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');

            set(axVec{3},'Position',[0.08 0.1 0.76 0.2])
            set(axVec{2},'Position',[0.48 0.400952380952381 0.4 0.2])
        case 6
            text(axVec{2},183,30,Ages{3},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},183,70,Ages{4},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{2},183,110,Ages{5},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{2},183,145,Ages{6},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},197,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},197,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{2},197,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');

            text(axVec{3},167,20,Ages{2},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},167,70,Ages{3},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},167,112,Ages{4},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{3},167,160,Ages{5},'fontsize',12,'color','k','fontweight','bold');
            text(axVec{3},194,90,Ages{6},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},194,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
            text(axVec{3},194,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');


             set(axVec{3},'Position',[0.1 0.1 0.78 0.2])
            set(axVec{2},'Position',[0.52 0.400952380952381 0.42 0.2])
    end

    set(gcf,'Position',[520 199 736 598]);
    printGraph(['../graphs/ParetoFront_R0=',num2str(R0)]);
    close all;
end

ix=ix+1;

%% Recovered population
if whichGraphsToProduce(ix)
    recoveredprct10=10
    recoveredprct20=20
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')

    % Pareto front
    axVec{1}=nexttile([1,2]);

    % Default
    defineColors;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');

    % 10% recovered
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct10,betaVac,nuVac,effVac,'All',10)
    set(h,'Color',red,'LineStyle','-');

    % 20% recovered
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'All',20)
    set(h,'Color',blue,'LineStyle','-.');

    xlim([158 220]);

    legend('None recovered','10% recovered','20% recovered','location','southwest');

    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct10,betaVac,nuVac,effVac,'All',10)
    title('10% recovered');

    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'All',20)
    title('20% recovered');

    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    lg = legend(Ages,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    lg.Layout.Tile='south'
    linkaxes([ axVec{2} axVec{3} ],'y');
    set(gcf,'Position',[520 392 676 405]);
    shg;
    printGraph('../graphs/ParetoFront_recovered');
    close all;


    %%
    t=tiledlayout(3,1, 'TileSpacing','compact','Padding','none')

    % Pareto front
    axVec{1}=nexttile(1);
%     VcPrct80=80;VcPrct66=50;
%     presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
%     h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All');hold on;
%     set(h,'Color',0.7*[1 1 1],'LineStyle','--');
%     presentationOptions.presentUniform=false;
%     presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
%     h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'Above20')
%     set(h,'Color',0.7*[1 1 1],'LineStyle','-')

    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'Above20',20)
    g=gca;xlimSave=g.XLim;

    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'Above10',20)
    set(h,'Color',yellow,'LineStyle','-.');
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'All',20)
    set(h,'Color',purple,'LineStyle','--');
   

    g=gca;
    g.Children
        lg=legend([g.Children(4) g.Children(3) g.Children(2) g.Children(1)],'Random allocations','Pareto front - ages 20 and older eligible (20% recovered)','Pareto front - ages 10 and older eligible (20% recovered)','Pareto front - all eligible (20% recovered)','location','southwest','autoupdate','off','fontsize',12);

%     lg=legend([g.Children(4) g.Children(3) g.Children(5) g.Children(2) g.Children(1) g.Children(6)],'Random allocations','Pareto front - ages 20 and older eligible (20% recovered)','Pareto front - ages 20 and older eligible (none recovered)','Pareto front - ages 10 and older eligible (20% recovered)','Pareto front - all eligible (20% recovered)','Pareto front - all eligible (none recovered)','location','southwest','autoupdate','off','fontsize',12);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));


    axVec{2}=axes('Position',[0.717857142857143 0.400952380952381 0.27 0.18]);
    presentationOptions.presentScatterplot=false ;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'Above20',20)
    title('Eligibilty - ages 20 and above')
    grid on

    axVec{3}=axes('Position',[0.4527857142857143 0.1 0.5351 0.18]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,'All',20)
    grid on
    title('Vaccine eligibilty - all ages')
    xlim(axVec{1},[160 225]);ylim(axVec{1},[-0.5,3.5]);set(axVec{1},'ytick',[0:1:3])
    % Manually add age group labels

    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};

    lg = legend(axVec{2},Ages,'fontsize',10,'Orientation','vertical');%legend(nexttile(2), [line1,line2,line3,line4]);
    lg.Location = 'eastoutside';
    set(lg,...
        'Position',[0.086735993349466 0.390763760632677 0.112687813021703 0.2],...
        'FontSize',10);
%     ylim(axVec{1},[-0.5,0.8]);
    set(gcf,'Position',[520 234 599 563]);shg
    printGraph('../graphs/ParetoFront_recovered20');
    close all;
end
return