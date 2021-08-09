function produceAllGraphs(collectData)
if collectData
    display('Computing data for Pareto front (this is not a quick computation)')
    computeParetoFront;
end

close all;
defineColors;
produceGraph=[false, % default
    false, % all-or-none
    false, % hesistancy]
    false, % efficacy
    false, % susceptibility
    false, % coverage
    false, % comparison to 20+
    false] % recovered 
produceGraph(:)=true;
ix=1;
%% Default scenario
susceptibilityFactor=1;maxPrct=100;VcPrct=73.2;betaVac=0.2;nuVac=1;effVac=1-0.05/betaVac;
presentationOptions.presentUniform=false;
if produceGraph(ix)
    t=tiledlayout(2,1, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile(1);
    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    
%     set(axVec{1},'Position',[0.0930232558139535 0.762767710049423 0.317013463892289 0.219110378912685]);
    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=true;
    presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above20')
    

    % Distributions along Pareto front
    axVec{2}=nexttile(2);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    xlim([75 163])
    % Link axes
    linkaxes([axVec{1} axVec{2} ],'x');%ylabel(axVec{1},'Non-vaccinated Infected','fontsize',11);
    
    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    text(axVec{2},120,15,Ages{2},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},120,40,Ages{3},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},120,75,Ages{4},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{2},120,115,Ages{5},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{2},153,85,Ages{6},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},153,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},153,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},153,174,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');
    set(gcf,'Position',[520 190 817 607])
    return
    printGraph('../graphs/ParetoFront');
    close all;
end
ix=ix+1;
%% All-or-none
if produceGraph(ix)
    betaVacAllorNone=0;nuVacAllorNone=0.8;
    t=tiledlayout(3,1, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile(1);
    betaVac,nuVac
    
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');
    
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'Above20')
    set(h,'Color',0.7*[1 1 1],'LineStyle','-')
    
    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'Above20')
    g=gca;xlimSave=g.XLim;
    

    
    %         presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    %         h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'Above10')
    %         set(h,'Color',yellow,'LineStyle','-.');
    %         xlim(xlimSave);
    
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac)
    set(h,'Color',purple,'LineStyle','--');
    xlim([45 135]);
    
    g=gca;
            g.Children
    lg=legend([g.Children(3) g.Children(2) g.Children(4) g.Children(1) g.Children(5)],'Random allocations','Pareto front - ages 20 and older (All-or-none)','Pareto front - ages 20 and older (leaky)','Pareto front - all ages (all-or-none)','Pareto front - all ages (leaky)','location','southwest','autoupdate','off','fontsize',12);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));

    %         legend(g.Children([1 2 3],'Random allocations','Homogenous allocation','Pareto front','Allocation minimizing infections','Allocation minimizing mortality','Section of Pareto front - all ages eligible','location','southwest','autoupdate','off','fontsize',12)
    %         legend([g.Children(9) g.Children(6) g.Children(2) g.Children(1)],'Random allocations - ages 20 and older eligible','Pareto front - ages 20 and older eligible','Pareto front - ages 10 and older eligible','Pareto front - all ages eligible','location','southwest','autoupdate','off','fontsize',12)
%     legend([g.Children(9) g.Children(6) g.Children(1)],'Random allocations - ages 20 and older eligible','Pareto front - ages 20 and older eligible','Pareto front - all ages eligible','location','southwest','autoupdate','off','fontsize',12)
    % Distributions along Pareto front
    %axVec{2}=nexttile(3);
    
    %         axVec{2}=axes('Position',[0.517857142857143 0.330952380952381 0.456696428571429 0.1375]);
    axVec{2}=axes('Position',[0.447857142857143 0.400952380952381 0.305 0.2]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,'Above20')
    title('Vaccine restricted to ages 20 and above')
    grid on
    %axVec{3}=nexttile(4);
    
    axVec{3}=axes('Position',[0.094 0.1 0.6589 0.2]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac)
    grid on
    title('No age restriction for vaccine allocation')
    % Link axes
    %         linkaxes([axVec{1} axVec{2} axVec{3} ],'x');%ylabel(axVec{1},'Non-vaccinated Infected','fontsize',11);
    set(gcf,'Position',[520 234 599 563]);
    % Manually add age group labels
        xlim(axVec{1},[55 163]);

    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
            lg = legend(axVec{2},Ages,'fontsize',10,'Orientation','vertical');%legend(nexttile(2), [line1,line2,line3,line4]);
        lg.Location = 'eastoutside';
    set(lg,...
    'Position',[0.086735993349466 0.390763760632677 0.112687813021703 0.214031971580817],...
    'FontSize',10);
    ylim(axVec{1},[-0.6 2]);set(axVec{1},'ytick',0:0.5:2)
    shg;
    printGraph('../graphs/ParetoFront_AllorNone');
    close all;
end
ix=ix+1;

%% Hesistancy
if produceGraph(ix)
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
    
    xlim([74 161]);
    
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
if produceGraph(ix)
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile([1,2]);
    
    % Default
    defineColors;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');
    
    % 90% efficacy
    betaVacVE85=0.15;effVacVE85=1-0.05/betaVacVE85;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacVE85,nuVac,effVacVE85)
    set(h,'Color',red,'LineStyle','-');
    
    % 70% efficacy
    betaVacVE75=0.25;effVacVE75=1-0.05/betaVacVE75;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVacVE75,nuVac,effVacVE75)
    set(h,'Color',blue,'LineStyle','-.');
    
    xlim([40 180]);ylim([0 2.1])
    
    legend('80% Vaccine Efficacy (as in Figure 5)','85% Vaccine Efficacy','75% Vaccine Efficacy','location','northeast');
    
    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,100,VcPrct,betaVacVE85,nuVac,effVacVE85)
    title('85% Vaccine Efficacy');
    
    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,100,VcPrct,betaVacVE75,nuVac,effVacVE75)
    title('75% Vaccine Efficacy');
    
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

%% Susceptiblity
if produceGraph(ix)
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile([1,2]);
    
    % Default
    defineColors;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');
    
    % 1.5 susceptibility factor
    susceptibilityFactor15=1.5;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor15,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',red,'LineStyle','-');
    
    % 2 susceptibility factor
    susceptibilityFactor20=2;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor20,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',blue,'LineStyle','-.');
    
    xlim([65 170]);
    
    legend('Susceptiblity profile as in Figure 4','Increased susceptiblity by factor of 1.5','Increased susceptiblity by factor of 2','location','northeast');
    
    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor15,maxPrct,VcPrct,betaVac,nuVac,effVac)
    title('Increased susceptiblity by factor of 1.5');
    
    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor20,maxPrct,VcPrct,betaVac,nuVac,effVac)
    title('Increased susceptiblity by factor of 2');
    
    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    lg = legend(Ages,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    %     lg.Location = 'southoutside';\
    lg.Layout.Tile='south'
    
    set(gcf,'Position',[520 392 676 405]);
    shg;
  printGraph('../graphs/ParetoFront_susceptibility');
  close all;
end

ix=ix+1;

%% Vaccination coverage
if produceGraph(ix)
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile([1,2]);
    
    % Default
    defineColors;
    presentationOptions.presentUniform=false;presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');
    
    % 80% coverage
    VcPrct80=80;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct80,betaVac,nuVac,effVac)
    set(h,'Color',red,'LineStyle','-');
    
    % 66% coverage
    VcPrct66=66;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac)
    set(h,'Color',blue,'LineStyle','-.');
    
    xlim([35 185]);ylim([0 2.3])
    g=gca;
    legend([g.Children(1),g.Children(3),g.Children(2)],'50% vaccine coverage','55% vaccine coverage  (as in Figure 5)','60% vaccine coverage','location','northeast');
    
    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac)
    title('50% vaccine coverage');
    
    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct80,betaVac,nuVac,effVac)
    title('60% vaccine coverage');
    
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
if produceGraph(ix)
    t=tiledlayout(3,1, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile(1);
    
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above10')
    set(h,'Color',yellow,'LineStyle','-.');
    
    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    set(h,'Color',purple,'LineStyle','--');
%     xlim(xlimSave);
    betaVac=0.2;
    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=true;
    presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above20')
    g=gca;xlimSave=g.XLim;
    
    
    xlim(xlimSave);
    
    
    
    g=gca;
    %         legend(g.Children([1 2 3],'Random allocations','Homogenous allocation','Pareto front','Allocation minimizing infections','Allocation minimizing mortality','Section of Pareto front - all ages eligible','location','southwest','autoupdate','off','fontsize',12)
    %          legend([g.Children(2) g.Children(3) g.Children(3) g.Children(3)],'Random allocations - all eligible','Pareto front - all eligible','Random allocations - ages 20 and older eligible','Pareto front - ages 10 and older eligible','Pareto front - all ages eligible','location','southwest','autoupdate','off','fontsize',12)
    g.Children
    %lg=legend([g.Children(4) g.Children(2) g.Children(5) g.Children(3) g.Children(1)],'Random allocations - all eligible','Random allocations - ages 20 and older eligible','Pareto front - all eligible','Pareto front - ages 10 and older eligible','Pareto front - ages 20 and older eligible','location','southwest','autoupdate','off','fontsize',12);
    lg=legend([g.Children(10) g.Children(8) g.Children(11) g.Children(9) g.Children(5)],'Random allocations - all ages','Random allocations - ages 20 and older','Pareto front - all ages','Pareto front - ages 10 and older','Pareto front - ages 20 and older','location','southwest','autoupdate','off','fontsize',12);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));

    axVec{2}=axes('Position',[0.53 0.400952380952381 0.456696428571429 0.2]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above20')
    title('Vaccine allocation restricted to ages 20 and above')
    grid on
    
    axVec{3}=axes('Position',[0.08 0.1 0.89 0.2]);
    
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac)
    grid on
    title('No age restriction on vaccine allocation')
    xlim(axVec{1},[75 163]);
    set(gcf,'Position',[520 234 610 563]);
    % Manually add age group labels
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    text(axVec{2},120,30,Ages{3},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},120,70,Ages{4},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{2},120,110,Ages{5},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{2},120,145,Ages{6},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},154.5,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},154.5,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{2},154.5,174.2,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');
    
    text(axVec{3},80,15,Ages{2},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},80,55,Ages{3},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},80,100,Ages{4},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{3},80,145,Ages{5},'fontsize',12,'color','k','fontweight','bold');
    text(axVec{3},154.5,90,Ages{6},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},154.5,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
    text(axVec{3},154.5,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
%     set(gcf,'Position',[520 164 517 633])
    set(gcf,'Position',[520 199 736 598]);return
  printGraph('../graphs/ParetoFront');
  close all;
end

ix=ix+1;

%% Recovered population
if produceGraph(ix)
    VcPrct66=50;
    t=tiledlayout(2,2, 'TileSpacing','compact','Padding','none')
    
    % Pareto front
    axVec{1}=nexttile([1,2]);
    
    % Default
    defineColors;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac)
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');
    
    % 10% recovered
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'All',10)
    set(h,'Color',red,'LineStyle','-');
    
    % 20% recovered
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'All',20)
    set(h,'Color',blue,'LineStyle','-.');
    
    xlim([80 220]);
    
    legend('None recovered','10% recovered','20% recovered','location','northwest');
    
    axVec{2}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'All',10)
    title('10% recovered');
    
    axVec{3}=nexttile;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'All',20)
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
    VcPrct80=80;VcPrct66=50;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct80,betaVac,nuVac,effVac,'All');hold on;
    set(h,'Color',0.7*[1 1 1],'LineStyle','--');
     presentationOptions.presentUniform=false;
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct80,betaVac,nuVac,effVac,'Above20')
    set(h,'Color',0.7*[1 1 1],'LineStyle','-')
    
    presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'Above20',20)
    g=gca;xlimSave=g.XLim;
    
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'Above10',20)
    set(h,'Color',yellow,'LineStyle','-.');
presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'All',20)
    set(h,'Color',purple,'LineStyle','--');
    axis([25 125 -0.66 0.85]);set(gca,'ytick',[0:0.2:0.8])
    
    g=gca;
    g.Children
    lg=legend([g.Children(4) g.Children(3) g.Children(5) g.Children(2) g.Children(1) g.Children(6)],'Random allocations','Pareto front - ages 20 and older eligible (20% recovered)','Pareto front - ages 20 and older eligible (none recovered)','Pareto front - ages 10 and older eligible (20% recovered)','Pareto front - all eligible (20% recovered)','Pareto front - all eligible (none recovered)','location','southwest','autoupdate','off','fontsize',12);
    set(lg.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
    
    
    axVec{2}=axes('Position',[0.717857142857143 0.400952380952381 0.27 0.2]);
    presentationOptions.presentScatterplot=false ;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'Above20',20)
    title('Eligibilty - ages 20 and above')
    grid on
    
    axVec{3}=axes('Position',[0.4527857142857143 0.1 0.5351 0.2]);
    presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
    h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct66,betaVac,nuVac,effVac,'All',20)
    grid on
    title('Vaccine eligibilty - all ages')
    xlim(axVec{1},[35 185]);
    % Manually add age group labels
    
    Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};

        lg = legend(axVec{2},Ages,'fontsize',10,'Orientation','vertical');%legend(nexttile(2), [line1,line2,line3,line4]);
        lg.Location = 'eastoutside';
    set(lg,...
    'Position',[0.086735993349466 0.390763760632677 0.112687813021703 0.214031971580817],...
    'FontSize',10);
ylim(axVec{1},[-0.5,0.8]);
     set(gcf,'Position',[520 234 599 563]);
printGraph('../graphs/ParetoFront_recovered20');
close all;
end