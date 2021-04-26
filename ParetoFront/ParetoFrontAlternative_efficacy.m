function ParetoFront(vaccineRange)
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
addpath('../Codes/AllCountries/AuxilaryFunctions/');
close all;
vaccineRange='All';
%% Parameters
collectData=false;

R0=3;betaVac=0.2; effVac=0.75;fname='dataParetoFront732';
R0=3;betaVac=0.1; effVac=0.5;fname='dataParetoFront732_betaVac01';
R0=3;betaVac=0.5; effVac=0.9;fname='dataParetoFront732_betaVac05';
%recoveredprct=8;infected_nv_prct=0.0045;infected_v_prct=0.0045;
recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;
% The percent of 20+ population that have vaccines avaliable
VcPrct=73.2; maxPrct=100; % Maximal vaccination per age group
%
% recoveredprct=8;
% infected_nv_prct=0;
% infected_v_prct=0;
% VcPrct=80; % The percent of 20+ population that have vaccines avaliable
% betaVac=0.2; % Relative suscebtibility of vaccined individuals
% effVac=0.8; % Efficacy of vaccine in preventing disease
% maxPrct=95; % Maximal vaccination per age group
% R0=3;
%% Select Country
%countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"};
%ix=2;
%country=countryList{ix};
country="USA";
%countryList{ix};
%countryData=load(join(['../countryData/',country,'_data.mat'],''));
vaccineRangeVec={'above20','above16','above12','above6','All'};
% vaccineRange=vaccineRangeVec{1};
%% Prepare data
[uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,0.99*VcPrct,betaVac,effVac,maxPrct,true);

%% Define optimization problem
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',1e4,'Display','none','Algorithm','sqp');%,'Display','iter');
problem.solver = 'fmincon';
problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
problem.A=[];problem.B=[];
problem.lb=0*upperBound;
problem.ub=upperBound;

defineColors;
%% Run loop
if collectData
    
    [dummy,overallInfected_uniform,overallMortality_uniform,data_uniform]=computeFinalSize(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    
    Msample=5000;M=800;
    %% First solve for infection minimizing allocation
    w=0;
    [result,overallInfected_spreaders,overallMortality_spreaders,data_spreaders]=computeFinalSize(spreadersAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    x=fmincon(@(x)computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
    [dummy,overallInfected(1),overallMortality(1),data{1}]=computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    distribution{1}=x;
    parfor ix=1:Msample
        y=randomizePoint(distribution{1},Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.15);
        sampleAllocation{ix}=y;
        [dummy,sample_overallInfected(ix),sample_overallMortality(ix),dummy]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    end
    
    [minResult,minIdx]=min(sample_overallInfected);
    if minResult<overallInfected(1)
        [dummy,overallInfected(1),overallMortality(1),data{1}]=computeFinalSize(sampleAllocation{minIdx},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
        distribution{1}=sampleAllocation{minIdx};
        display('Improved at infection minimizing end');
    end
    
    %% Next solve for mortality minimizing allocation
    w=1;
    [result,overallInfected_riskFocused,overallMortality_riskFocused,data_riskFocused]=computeFinalSize(riskFocusedAllocation,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    x=fmincon(@(x)computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
    [dummy,overallInfected(M),overallMortality(M),data{M}]=computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    distribution{M}=x;
    parfor ix=1:Msample
        y=randomizePoint(distribution{M},Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
        sampleAllocation{ix}=y;
        [dummy,sample_overallInfected(ix),sample_overallMortality(ix),dummy]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    end
    
    [minResult,minIdx]=min(sample_overallMortality);
    if minResult<overallMortality(M)
        [dummy,overallInfected(M),overallMortality(M),data{M}]=computeFinalSize(sampleAllocation{minIdx},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
        distribution{M}=sampleAllocation{minIdx};
        display('Improved at mortality minimizing end');
    end
    
    minInfected=overallInfected(1)
    maxInfected=overallInfected(M)
    
    infectionVec=linspace(minInfected,maxInfected,M);
    
    %% For each number of infections, minimize mortality
    for ix=2:M-1
        infectedConstraint=infectionVec(ix);
        x=fmincon(@(x)computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1),distribution{ix-1},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac),problem.options);
        [dummy,overallInfected(ix),overallMortality(ix),dummy]=computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
        distribution{ix}=x;
    end
    
    % Sweep back
    x=distribution{M};
    for ix=M-1:-1:2
        infectedConstraint=infectionVec(ix);
        x=fmincon(@(x)computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1),distribution{ix+1},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac),problem.options);
        [dummy,overallInfected_aux(ix),overallMortality_aux(ix),dummy]=computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
        if overallMortality_aux(ix)<overallMortality(ix)
            distribution{ix}=x;
            overallMortality(ix)=overallMortality_aux(ix);
            display(['Managed to improve ix=',num2str(ix)]);
        end
    end
    % Sample points
    parfor ix=1:Msample
        y=randomizePoint(distribution{1+mod(ix,M)},Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.25);
        sampleAllocation1{ix}=y;
        [dummy,sample_overallInfected1(ix),sample_overallMortality1(ix),dummy]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
        
        y=randomizePoint(uniformAllocation,Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
        sampleAllocation2{ix}=y;
        [dummy,sample_overallInfected2(ix),sample_overallMortality2(ix),dummy]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
    end
    sample_overallInfected=[sample_overallInfected1 sample_overallInfected2];
    sample_overallMortality=[sample_overallMortality1 sample_overallMortality2];
    
    save([fname,'_',vaccineRange,'_R=',num2str(R0)]);
end
load([fname,'_',vaccineRange,'_R=',num2str(R0)]);

%% Present Data
close all;
t=tiledlayout(2,1, 'TileSpacing','compact','Padding','none')
defineColors;
%green=[0.4660, 0.6740, 0.1880];orange=[0.8500, 0.3250, 0.0980];blue=[0, 0.4470, 0.7410];yellow=[0.9290, 0.6940, 0.1250];red=[0.6350, 0.0780, 0.1840];
axVec{1}=nexttile(1);
%scatter1=scatter([sample1_overallInfected_nv sample2_overallInfected_nv],[sample1_overallMortality sample2_overallMortality],8,'MarkerFaceColor',[1, 0.8, 0.3],'MarkerEdgeColor',[1, 0.8, 0.3])
scatter1=scatter(sample_overallInfected/1e6,sample_overallMortality/1e6,8,'MarkerFaceColor',yellow,'MarkerEdgeColor',yellow)
alpha(scatter1,0.25)
hold on;plot(overallInfected_uniform/1e6,overallMortality_uniform/1e6,'d','color',blue,'markersize',5,'linewidth',2)
text(overallInfected_uniform/1e6,overallMortality_uniform/1e6,'Homogenous','HorizontalAlignment','left','VerticalAlignment','bottom','color',blue,'fontsize',12,'FontWeight','bold')

%[k,av] = convhull([overallInfected_nv' overallMortality']);
%plot(overallInfected_nv(k(1:end-1)),overallMortality(k(1:end-1)),'b','linewidth',1.2);
overallMortality=overallMortality/1e6;
overallInfected=overallInfected/1e6;
plot(overallInfected,overallMortality,'k','linewidth',1.5);hold on;

scatter(overallInfected(1),overallMortality(1),35,'s','filled','MarkerFaceColor',green,'MarkerEdgeColor',green);hold on;

% plot(overallInfected_nv(1),overallMortality(1),'s','color',green,'linewidth',1.2,'markersize',14);
text(overallInfected(1),overallMortality(1),'Allocation minimizing infections ','HorizontalAlignment','left','VerticalAlignment','bottom','color',green,'fontsize',12,'FontWeight','bold')
% plot(overallInfected_nv( end),overallMortality(end),'*','color',orange,'linewidth',1.2,'markersize',14);
scatter(overallInfected( end),overallMortality(end),35,'^','filled','MarkerFaceColor',orange,'MarkerEdgeColor',orange);hold on;

text(overallInfected(end),overallMortality(end),'Allocation minimizing mortality','HorizontalAlignment','right','VerticalAlignment','top','color',orange,'fontsize',12,'FontWeight','bold')

xlabel('Overall Infected');ylabel('Overall Mortality');
xlim([min(overallInfected)*0.999,max([overallInfected overallInfected_riskFocused/1e6])*1.001]);
box on;    ax = gca;ax.YRuler.Exponent = 0;ytickformat('%,7.0f');ax.XRuler.Exponent = 0;xtickformat('%,7.0f')
ax.YRuler.Exponent = 0;ytickformat('%g M');
ax.XRuler.Exponent = 0;xtickformat('%g M');
% hold on;plot([overallInfected_nv_spreaders overallInfected_nv_risk],[overallMortality_spreaders overallMortality_risk],'+','color',red,'markersize',5,'linewidth',2)

% text(overallInfected_nv_spreaders,overallMortality_spreaders,'Spreaders','HorizontalAlignment','left','VerticalAlignment','bottom','color',red,'fontsize',12,'FontWeight','bold')
% text(overallInfected_nv_risk,overallMortality_risk,'High-risk   ','HorizontalAlignment','center','VerticalAlignment','bottom','color',red,'fontsize',12,'FontWeight','bold')
legend('Random allocations','Homogenous allocation','Pareto front','Allocation minimizing infections','Allocation minimizing mortality','location','north','autoupdate','off','fontsize',12)
%axis([55.8e6 64.5e6 9e4 6e5])


axVec{2}=nexttile(2);
%M=numel(k)-1;
k=1:M;
for ix=1:M
    vacDist=0*r;vacDist(adultAges)=distribution{k(ix)};
    
    allocationData(ix,:)=vacDist.*Ni;%[sum(distribution{k(ix)}(ageAbove60)) sum(distribution{k(ix)}(age40to60)) sum(distribution{k(ix)}(age20to40))];
end
% allocationData(M+1,:)=riskFocusedAllocation.*Nadult;%[sum(riskFocusedAllocation(ageAbove60)) sum(riskFocusedAllocation(age40to60)) sum(riskFocusedAllocation(age20to40))]
%load('colorbrewer.mat');
%load('Glasbey_bw.mat');
%load('Glasbey_Category10.mat');
%C = parula(9);%colorbrewer.qual.Set1{9};
%C(5,:)=yellow*255;
%C(9,:)=purple*255;
h=area(overallInfected,allocationData/1e6);%axis([0 1 0 vaccinesLeftToDistibute*1.05/1e6]);
for ix=1:9
    h(ix).FaceColor=lineColors(ix,:);
end

% for ix=1:7
%     currY=currY+dataMat(a,ix);
%     if dataMat(a,ix)>sum(dataMat(a,:))/50
%         text(data.Rvec(a),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
%     end
% end

%area(overallMortality,allocationData/1e6);%axis([0 1 0 vaccinesLeftToDistibute*1.05/1e6]);
xlabel('w');title(['Allocation of vaccines among age groups']);
% legend('Ages 60+','Ages 40-59','Ages 20-39','location','best');

ylabel('Number of vaccines');ax.YRuler.Exponent = 0;ytickformat('%g M');%ytickformat('%,11.0f');
xlabel('Overall Infected');ax.XRuler.Exponent = 0;xtickformat('%g M');
box on;    %ax = gca;ax.XRuler.Exponent = 0;xtickformat('%,7.0f')
grid on
%
%printGraph(['Allocation_',vaccineRange])
linkaxes([axVec{1} axVec{2} ],'x');%ylabel(axVec{1},'Non-vaccinated Infected','fontsize',11);
%xlim([min(overallInfected_nv)*0.999,max([overallInfected_nv overallInfected_nv_risk/1e6])*1.001]);
xlim([min(overallInfected)*0.999,max([overallInfected 0*overallInfected_riskFocused/1e6])*1.001]);
ylim([0 vaccinesLeftToDistibute/1e6]);
%printGraph([fname,'_',vaccineRange,'_R=',num2str(R0)]);
shg

Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
%  legend('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+')
text(axVec{2},120,15,Ages{2},'fontsize',12,'color','w','fontweight','bold');
text(axVec{2},120,45,Ages{3},'fontsize',12,'color','w','fontweight','bold');
text(axVec{2},120,82,Ages{4},'fontsize',12,'color','k','fontweight','bold');
text(axVec{2},120,125,Ages{5},'fontsize',12,'color','k','fontweight','bold');
text(axVec{2},156,125,Ages{7},'fontsize',12,'color','w','fontweight','bold');
text(axVec{2},156,155,Ages{8},'fontsize',12,'color','w','fontweight','bold');
text(axVec{2},156,91,Ages{6},'fontsize',12,'color','w','fontweight','bold');
text(axVec{2},156,174,['  ',Ages{9}],'fontsize',12,'color','w','fontweight','bold');
% annotation('textarrow',[0.955357142857143 0.985714285714285],...
%    [0.171428571428571 0.133333333333334],'Color',[1 1 1]);
set(gcf,'Position',[520 190 817 607])
% printGraph('ParetoFront')
return

%function x=randomizePoint(x,N,Nadult,upperBound,vaccinesLeftToDistibute,delta)
x=max(x+delta*(rand(size(x))-0.5),0);
x=vaccinesLeftToDistibute/sum(x.*Nadult)*x;

idxVec=ones(size(x));
idx=find(x>upperBound);

while numel(idx)>0
    vacExtra=(x(idx)-upperBound(idx))'*Nadult(idx);
    x(idx)=upperBound(idx)-1e-15;
    idxVec(idx)=false;
    n_idx=find(idxVec);
    x(n_idx)=x(n_idx)+(upperBound(n_idx).*Nadult(n_idx)/sum(upperBound(n_idx).*Nadult(n_idx))*vacExtra)./Nadult(n_idx);
    idx=find(x>upperBound);
end


function [c,ceq] = mycon(x,infectedContraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac)

[dummy,overallInfected,dummy,dummy]=computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);

c = [];     % Compute nonlinear inequalities at x.
ceq = infectedContraint-overallInfected;   % Compute nonlinear equalities at x.
return