function h=presentParetoFrontData(presentationOptions,susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,vaccineRange,recoveredprct,fnameParm)
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
R0=6;presentScatterplot=true;

useconvhull=false;
if nargin<8
    vaccineRange='All';
end

if nargin==10
    fname=fnameParm;R0=6;useconvhull=true;
end
data=load('./data/dataParetoFront_Israel_All_R=6.mat');

% data=load('above12_R=15.mat');
% country="USA";
% data=load(['./data/',fname,'_',vaccineRange,'_R=',num2str(R0)]);
data.vaccinesLeftToDistibute/sum(data.Ni)
data.vaccinesLeftToDistibute/sum(data.Ni(3:9))

defineColors;
% alreadyRecovered=0*recoveredprct*sum(data.Ni)/100;
% 
% data.sample_overallInfected=data.sample_overallInfected+alreadyRecovered;
% data.sample_overallInfected1=data.sample_overallInfected1+alreadyRecovered;
% data.sample_overallInfected2=data.sample_overallInfected2+alreadyRecovered;
% data.overallInfected=data.overallInfected+alreadyRecovered;
% data.overallInfected_uniform=data.overallInfected_uniform+alreadyRecovered;
% data.infectionVec=data.infectionVec+alreadyRecovered;

if useconvhull
    [k,av] = convhull([data.infectionVec' data.overallMortality']);
    k(end)=[];
    [dummy,I1]=min(data.infectionVec(k(:)));
    [dummy,I2]=max(data.infectionVec(k(:)));
    k=k(I1:I2);
    
    [dummy,I1]=min(data.overallMortality(k(:)));
    [dummy,I2]=max(data.overallMortality(k(:)));
    k=k(I2:I1);
    data.overallMortality=data.overallMortality(k(:));
    data.infectionVec=data.infectionVec(k(:));
end
switch vaccineRange
    case 'All'
        randColor=turquoise;alphaValue=0.25;
        vaccineRangeTxt='All elibigle';
        markerType='s'
    case 'Above20'
        randColor=yellow;alphaValue=0.25;
        vaccineRangeTxt='Eligibility above 20';
        markerType='o'
            case 'above16'
        randColor=yellow;alphaValue=0.25;
        vaccineRangeTxt='Eligibility above 16';
        markerType='o'
end
%% Present Data
if presentationOptions.presentScatterplot
    scatter1=scatter(data.sample_overallInfected/1e6,data.sample_overallMortality/1e3,8,'filled',markerType,'MarkerFaceColor',randColor,'MarkerEdgeColor',randColor)
    alpha(scatter1,alphaValue)
     hold on;
end

if presentationOptions.presentUniform
        hold on;plot(data.overallInfected_uniform/1e6,data.overallMortality_uniform/1e3,'d','color',blue,'markersize',5,'linewidth',2)
    text(data.overallInfected_uniform/1e6,data.overallMortality_uniform/1e3,{'Homogenous',['(',vaccineRangeTxt,')']},'HorizontalAlignment','center','VerticalAlignment','bottom','color',blue,'fontsize',12,'FontWeight','bold')
% plot(data.overallInfected_uniform/1e6,data.overallMortality_uniform/1e6,'d','color',blue,'markersize',5,'linewidth',2)
%      text(data.overallInfected_uniform/1e6,data.overallMortality_uniform/1e6,'Homogenous','HorizontalAlignment','left','VerticalAlignment','bottom','color',blue,'fontsize',12,'FontWeight','bold')

end
overallMortality=data.overallMortality/1e3;
overallInfected=data.infectionVec/1e6;


if presentationOptions.presentParetoGraph
    
    h=plot(overallInfected,overallMortality,'k','linewidth',1.5);hold on;
    
    
     
    if presentationOptions.presentUniform
           scatter(overallInfected(1),overallMortality(1),35,'s','filled','MarkerFaceColor',green,'MarkerEdgeColor',green);hold on;
    
    scatter(overallInfected( end),overallMortality(end),35,'^','filled','MarkerFaceColor',orange,'MarkerEdgeColor',orange);hold on;

        text(overallInfected(1),overallMortality(1),'Allocation minimizing infections ','HorizontalAlignment','left','VerticalAlignment','bottom','color',green,'fontsize',12,'FontWeight','bold')
        text(overallInfected(end),overallMortality(end),'Allocation minimizing mortality','HorizontalAlignment','right','VerticalAlignment','top','color',orange,'fontsize',12,'FontWeight','bold')
    end
    xlabel('Overall Infected');ylabel('Overall Mortality');
    xlim([min(overallInfected)*0.999,max(overallInfected)*1.001]);
    box on;    ax = gca;ax.YRuler.Exponent = 0;ytickformat('%,7.0f');ax.XRuler.Exponent = 0;xtickformat('%,7.0f')
    ax.YRuler.Exponent = 0;ytickformat('%g K');
    ax.XRuler.Exponent = 0;xtickformat('%g M');
    legend('Random allocations','Homogenous allocation','Pareto front','Allocation minimizing infections','Allocation minimizing mortality','location','southwest','autoupdate','off','fontsize',12)
end

if presentationOptions.presentDistribution
    %M=numel(k)-1;
    if ~useconvhull
        k=1:data.M;
    end
    for ix=1:numel(k)
        vacDist=0*data.r;vacDist(data.adultAges)=data.distribution{k(ix)};
        
        allocationData(ix,:)=vacDist.*data.Ni;%[sum(distribution{k(ix)}(ageAbove60)) sum(distribution{k(ix)}(age40to60)) sum(distribution{k(ix)}(age20to40))];
    end
    
    h=area(overallInfected,allocationData/1e3);%axis([0 1 0 vaccinesLeftToDistibute*1.05/1e6]);
    for ix=1:9
        h(ix).FaceColor=lineColors(ix,:);
    end
    
    xlabel('w');title(['Allocation of vaccines among age groups']);
    % legend('Ages 60+','Ages 40-59','Ages 20-39','location','best');
    
    ylabel('Number of vaccines');ax.YRuler.Exponent = 0;ytickformat('%g K');%ytickformat('%,11.0f');
    xlabel('Overall Infected');ax.XRuler.Exponent = 0;xtickformat('%g M');
    box on;    %ax = gca;ax.XRuler.Exponent = 0;xtickformat('%,7.0f')
    grid on
    axis tight;
end

return