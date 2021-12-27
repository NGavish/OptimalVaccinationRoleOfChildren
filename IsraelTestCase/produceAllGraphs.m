function produceAllGraphs(collectData)
if collectData
    computeParetoFront;
end

%% Produce Pareto plot
susceptibilityFactor=1;VcPrct=80;maxPrct=100;betaVac=0.1;nuVac=1;effVac=0.5;
t=tiledlayout(2,1, 'TileSpacing','compact','Padding','none')
recoveredPrct=0.1386;
% Pareto front
axVec{1}=nexttile(1);
presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=true;
presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All',recoveredPrct,'dataParetoFront_Israel')

presentationOptions.presentScatterplot=true;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=true;
presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All',recoveredPrct,'dataParetoFront_Israel')

presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=true;presentationOptions.presentDistribution=false;presentationOptions.presentUniform=false;
h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All',0,'dataParetoFront_Israel')

delete(axVec{1}.Children(5));delete(axVec{1}.Children(5));delete(axVec{1}.Children(5));delete(axVec{1}.Children(5));
delete(axVec{1}.Children(7));delete(axVec{1}.Children(7));

legend([axVec{1}.Children(4), axVec{1}.Children(5) ],'Random allocations','Pareto front (all eligible)','location','best')

hold on;plot(1.7357*[1 1],[4 8.7],'-.','LineWidth',1.5,'Color',0.8*[1 1 1])
ylim([4 8.75])

% Distributions along Pareto front
axVec{2}=nexttile(2);
presentationOptions.presentScatterplot=false;presentationOptions.presentParetoGraph=false;presentationOptions.presentDistribution=true;
h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All',0,'dataParetoFront_Israel')


textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
legend(axVec{2},textAges,'AutoUpdate','off')

hold on;plot(1.7357*[1 1],[0 3e2],'-.','LineWidth',1.5,'Color',0.8*[1 1 1])

set(gcf,'Position',[520 190 817 607])
printGraph('../graphs/ParetoIsrael')
return