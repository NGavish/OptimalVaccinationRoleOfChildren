function presentData(collectData)

addpath('../Codes/AllCountries/AuxilaryFunctions/');
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')

defineColors

if collectData
    runComputationOfParetoFronts;
end
close all

textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
tiledlayout(3,3,'TileSpacing','compact');

for jx=1:40
    data=load(['./data/rnd',num2str(jx),'_dataParetoFront_susFactor=100_maxPrct100_VcPrct=732_betaVac=10_nuVac=10_recoveredPrct=0_All_R=6']);
    idx=find(data.overallInfected>100);
    a(jx)=data.overallInfected(idx(1));
    b(jx)=data.overallInfected(idx(end));
for kx=1:9
    nexttile(kx)
    for ix=idx
        dist(ix)=data.distribution{ix}(kx);
    end
    h=scatter(data.overallInfected(idx)/1e6,100*dist(idx),5,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue);hold on;
    alpha(h,0.1)
end
end
[min(a) max(a)]

data=load('dataParetoFront_susFactor=100_maxPrct100_VcPrct=732_betaVac=10_nuVac=10_recoveredPrct=0_All_R=6.mat');
idx=find(data.overallInfected>100);

for kx=1:9
    ax=nexttile(kx)

    for ix=idx
        dist(ix)=data.distribution{ix}(kx);
    end
    plot(data.overallInfected(idx)/1e6,100*dist(idx),'r','linewidth',2);
    title(textAges{kx});box on;grid on;
        xtickformat('%gM');xlabel('Overall infected');ytickformat('percentage');ylabel('% vaccinated')

end

h=scatter(160,100,5,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue);hold on;

lg=legend(ax,[ax.Children(2) ax.Children(1)],'Optimal allocations','Optimal allocations corresponding to perturbed contact matrices','Orientation','horizontal')
lg.Layout.Tile='south'
set(gcf,'Position',[520 281 751 516])
shg

printGraph('../graphs/CijSensitivity')