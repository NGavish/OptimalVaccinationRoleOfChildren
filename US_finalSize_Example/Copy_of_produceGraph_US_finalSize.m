function produceGraph_US_finalSize(allOrNone)
if nargin<1
    allOrNone=false;
end

% Relative susceptibility of vaccined individuals & efficacy of vaccine in preventing disease
R0=10;betaVac=0.1;nuVac=1;effVac=0.5;VcPrct=73.2;susceptibilityFactor=[1 1]; fname='US_finalSize_Example';
if allOrNone
    betaVac=0;nuVac=0.9;fname=[fname,'_allOrNone'];
end

maxPrct=100; % Maximal vaccination per age group
vaccineRangeVec={'All','above20','above20'};
%countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"};
country="USA";
countryData=load(join(['../countryData/',country,'_data.mat'],''));

% Prepare country specific data
Cij=countryData.contactMatrix;
IFR=countryData.IFR';
N=countryData.N;
Ni=N*countryData.agDist;
adultAgesVec{1}=1:9;adultAgesVec{2}=3:9;adultAgesVec{3}=3:9;
Nadult=Ni(adultAgesVec{2});
Vc=sum(Nadult)*VcPrct/100;

susceptibilityProfile=ones(1,9);susceptibilityProfile(1:2)=susceptibilityFactor;
Cij=diag(susceptibilityProfile)*Cij;

% Compute next generation matrix & scale contact matrix
Mij=diag(Ni)*Cij*diag(1./Ni);
[V,d]=eig(Mij);
Cij=Cij/d(1);
infectedDistribution=V(:,1)/sum(V(:,1)); % Extract age distribution of infectvies from largest e.v. of ngm

% Initialize susceptible, recovered, vaccinated and active infectives
recoveredprctVec=[0 0 20];infected_nv_prctVec=[0 0 0.25];infected_v_prctVec=[0 0 0.25];
figure('Renderer', 'painters', 'Position', [10 10 900 600])
tiledlayout(3,1, 'TileSpacing','compact','Padding','none');

for ix=1:3
    vaccineRange=vaccineRangeVec{ix};
    adultAges=adultAgesVec{ix};
    Nadult=Ni(adultAges);
    
    recoveredprct=recoveredprctVec(ix);infected_nv_prct=infected_nv_prctVec(ix);infected_v_prct=infected_v_prctVec(ix);
    
    r=recoveredprct/100*N*infectedDistribution./Ni;
    infected0_nv=infected_nv_prct/100*N*infectedDistribution./Ni;
    infected0_v=0*infected0_nv;
    infected0_v(adultAges)=N*infected_v_prct/100*infectedDistribution(adultAges)/sum(infectedDistribution(adultAges))./Nadult;
    v=zeros(size(Ni));
    s=1-r-infected0_v-infected0_nv-v;
    
    % Compute uniform allocation of vaccines
    vaccinesLeftToDistibute=Vc-sum(v);
    upperBound=max(s(adultAges)-(1-maxPrct/100),0);
    uniformAllocation=computeUniformAllocation(s,adultAges,Nadult,vaccinesLeftToDistibute,upperBound);
    
    axVec{ix}=nexttile(ix);
      
    [result,overallInfected_nv,overallFatality,data]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    presentData(data,Ni,overallFatality)
    overallInfected_nv
end
linkaxes([axVec{1} axVec{2} axVec{3}],'xy');
xlim(axVec{1},[0.5 9.5])
lg = legend(axVec{3},'Vaccinated susceptible','Vaccinated removed','Pre-existing immunity','Non-vaccinated removed','location','northeast','fontsize',10);%,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
axVec{1}.Children(13)
lg = legend([axVec{1}.Children(15) axVec{1}.Children(14) axVec{1}.Children(12)],'Vaccinated susceptible','Vaccinated removed','Non-vaccinated removed','location','northeast','fontsize',10);%,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
lg = legend([axVec{2}.Children(15) axVec{1}.Children(14) axVec{1}.Children(12)],'Vaccinated susceptible','Vaccinated removed','Non-vaccinated removed','location','northeast','fontsize',10);%,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
ylim(axVec{3},[0 50.1])
text(axVec{1},0.02,0.92,'A','units','normalized','fontsize',11)
text(axVec{2},0.02,0.92,'B','units','normalized','fontsize',11)
text(axVec{3},0.02,0.92,'C','units','normalized','fontsize',11)
set(gcf,'Position',[10 83 756 600]);
printGraph(['../graphs/',fname])
return
%%
function uniformAllocation=computeUniformAllocation(s,adultAges,Nadult,vaccinesLeftToDistibute,upperBound)
uniformAllocation=[s(adultAges).*Nadult/sum(s(adultAges).*Nadult)*vaccinesLeftToDistibute]./Nadult;

idxVec=ones(1,numel(adultAges));
idx=find(uniformAllocation>upperBound);

while numel(idx)>0
    vacExtra=(uniformAllocation(idx)-upperBound(idx))'*Nadult(idx);
    uniformAllocation(idx)=upperBound(idx)-1e-15;
    idxVec(idx)=false;
    n_idx=find(idxVec);
    uniformAllocation(n_idx)=uniformAllocation(n_idx)+(s(n_idx).*Nadult(n_idx)/sum(s(n_idx).*Nadult(n_idx))*vacExtra)./Nadult(n_idx);
    idx=find(uniformAllocation>upperBound);
end
return

function presentData(data,N,overallFatality)
sus_v=data(:,1);
infected_v=data(:,2);
r=data(:,3);
infected_nv=data(:,4);
sus_nv=data(:,5);
green=[0.4660 0.6740 0.1880];red=[0.8500 0.3250 0.0980];
darkgreen=[0 0.5 0];
blue=[0 0.4470 0.7410];
h=bar(data/1e6,'stacked');
h(1).FaceColor=darkgreen;h(2).FaceColor=green;h(3).FaceColor=blue;h(3).FaceColor=blue;h(4).FaceColor=red;h(5).FaceColor='w';
hT=[];
ix=4;
yCoord=N-sus_nv;
idx=find(infected_nv<0.5e7);yCoord(idx)=N(idx);
noidx=find(infected_nv>0.5e7)

text(h(ix).XData(noidx)+h(ix).XOffset,yCoord(noidx)/1e6,num2str(100*infected_nv(noidx)./N(noidx),'%2.1f%%'), ...
    'VerticalAlignment','top','horizontalalign','center','fontsize',11,'color','w','fontweight','bold');

text(h(ix).XData(idx)+h(ix).XOffset,yCoord(idx)/1e6,num2str(100*infected_nv(idx)./N(idx),'%2.1f%%'), ...
    'VerticalAlignment','bottom','horizontalalign','center','fontsize',11,'color','k','fontweight','bold');

ax = gca;ax.YRuler.Exponent = 0;ytickformat('%,2.2f M')
agDistLabels = {"0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"};
xlabel('Age groups')

set(gca,'xtick',1:numel(agDistLabels),'xticklabel',agDistLabels);
ax = gca;ax.YRuler.Exponent = 0;

ylim([0 1.09*max(N)]/1e6);ytickformat('%,4.0g M')
ylabel('Number of people');
jf=java.text.DecimalFormat;
set(gca,'TickLength',[0 0])

return
