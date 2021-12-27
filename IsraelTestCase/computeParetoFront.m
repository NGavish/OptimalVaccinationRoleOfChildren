function computeParetoFront
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
addpath('../Codes/AllCountries/AuxilaryFunctions/');

%% Select Country
country="ISR";
vaccineRange='All';
susceptibilityFactor=1;
maxPrct=93;
betaVac=0.1;
nuVac=1;
effVac=0.5;
R0=6;

%% Prepare data
%% Parameters
infected_nv_prct=0;infected_v_prct=0;
fname=['./data/dataParetoFront_Israel'];

%% Select Country
country="ISR";
vaccineRange='All';
recoveredprct=0;VcPrct=0;
%% Prepare data
[uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,0.99*VcPrct,betaVac,effVac,maxPrct,true);
r=([238560 274038 227435 185811 152459 105156  68902  36380  24539]-[6 8 38 73 1609 466 1241 2050 4070])'./Ni;
v=[0.0 659426.2 861606.1 873261.4 838685.1 677547.9 618353.9 443091.9 242818.2];
v=v'./Ni;
v=[0 60 66.11 69.62 75.86 80.26 85.42 87.99 80]'/100;
s=max(0,1-r-v-infected0_v-infected0_nv);
vaccinesLeftToDistibute=3e5;

%% Uniform
s=max(0,1-r-v-infected0_v-infected0_nv);

%% Define optimization problem
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',5e4,'ConstraintTolerance',1e-6,'StepTolerance',1e-10,'Display','none');%,'Algorithm','sqp','Display','none');%,'Display','iter');     
problem.solver = 'fmincon';
problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
problem.A=[];problem.B=[];
problem.lb=0*upperBound;
problem.ub=upperBound;

uniformAllocation =[s(adultAges).*Nadult/sum(s(adultAges).*Nadult)*problem.Beq]./Nadult;
spreadersAllocation=uniformAllocation;
riskFocusedAllocation=uniformAllocation;
defineColors;
%% Run loop

[dummy,overallInfected_uniform,overallMortality_uniform]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

Msample=5000;M=400;
%% First solve for infection minimizing allocation
w=0;
[result,overallInfected_spreaders,overallMortality_spreaders,data_spreaders]=computeFinalSize_generalized(spreadersAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
[x,fval,exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
[dummy,overallInfected(1),overallMortality(1),data{1}]=computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
distribution{1}=x;
parfor ix=1:Msample
    y=randomizePoint(distribution{1},Ni,Nadult,upperBound,vaccinesLeftToDistibute,(ix>1).*(0.05+0.5*mod(ix,3).^2));
    sampleAllocation{ix}=y;
    [dummy,sample_overallInfected(ix),sample_overallMortality(ix),dummy]=computeFinalSize_generalized(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
end

[minResult,minIdx]=min(sample_overallInfected);
if minResult<overallInfected(1)
    [dummy,overallInfected(1),overallMortality(1),data{1}]=computeFinalSize_generalized(sampleAllocation{minIdx},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    distribution{1}=sampleAllocation{minIdx};
    display('Improved at infection minimizing end');
end
infectedConstraint=overallInfected(1);
parfor kx=1:100
    y=randomizePoint(distribution{1},Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.25);
    [aux,overallMortality_sample(kx),exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac),problem.options);
    if exitflag<=0
        overallMortality_sample(kx)=Inf;
    end
    sampleAllocation1{kx}=aux;
end

[minResult,minIdx]=min(overallMortality_sample(1:100));
if overallMortality(1)>minResult
    distribution{1}=sampleAllocation1{minIdx};
    [dummy,overallInfected(1),overallMortality(1),dummy]=computeFinalSize_generalized(distribution{1},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    display(['Improved infection minimizing allocation']);
end

%% Next solve for mortality minimizing allocation
w=1;
[result,overallInfected_riskFocused,overallMortality_riskFocused,data_riskFocused]=computeFinalSize_generalized(riskFocusedAllocation,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

x=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
[dummy,overallInfected(M),overallMortality(M),data{M}]=computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
distribution{M}=x;
parfor ix=1:Msample
    y=randomizePoint(distribution{M},Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
    sampleAllocation{ix}=y;
    [dummy,sample_overallInfected(ix),sample_overallMortality(ix),dummy]=computeFinalSize_generalized(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
end

[minResult,minIdx]=min(sample_overallMortality);
if minResult<overallMortality(M)
    [dummy,overallInfected(M),overallMortality(M),data{M}]=computeFinalSize_generalized(sampleAllocation{minIdx},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    distribution{M}=sampleAllocation{minIdx};
    display('Improved at mortality minimizing end');
end

minInfected=overallInfected(1)
maxInfected=overallInfected(M)

infectionVec=linspace(minInfected,maxInfected,M);
w=1;
%% For each number of infections, minimize mortality
for ix=2:M-1
    infectedConstraint=infectionVec(ix);
    [x,overallMortality(ix),exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),distribution{ix-1},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac),problem.options);
    %         [dummy,overallInfected(ix),overallMortality(ix),dummy]=computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    distribution{ix}=x;
    if exitflag<=0
        distribution{ix}=distribution{ix-1};
        overallMortality(ix)=NaN;
    end
    parfor kx=1:50
        y=randomizePoint(x,Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.25);
        [aux,overallMortality_sample(kx),exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac),problem.options);
        if exitflag<=0
            overallMortality_sample(kx)=Inf;
        end
        sampleAllocation1{kx}=aux;
    end
    
    [minResult,minIdx]=min(overallMortality_sample);
    if isnan(overallMortality(ix)) || overallMortality(ix)>minResult
        distribution{ix}=sampleAllocation1{minIdx};
        [dummy,overallInfected(ix),overallMortality(ix),dummy]=computeFinalSize_generalized(distribution{ix},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
        display(['Improved at ',num2str(overallInfected(ix))]);
    end
end

% Sweep back
x=distribution{M};
for ix=M-1:-1:2
    infectedConstraint=infectionVec(ix);
    [x,overallMortality_aux(ix),exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),distribution{ix+1},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac),problem.options);
    if exitflag>0 && overallMortality_aux(ix)<overallMortality(ix)
        distribution{ix}=x;
        overallMortality(ix)=overallMortality_aux(ix);
        display(['* Managed to improve ix=',num2str(ix)]);
    end
end
% Sample points
parfor ix=1:Msample
    y=randomizePoint(distribution{1+mod(ix,M)},Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.05+0.5*mod(ix,3).^2);
    sampleAllocation1{ix}=y;
    [dummy,sample_overallInfected1(ix),sample_overallMortality1(ix),dummy]=computeFinalSize_generalized(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    
    y=randomizePoint(uniformAllocation,Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
    sampleAllocation2{ix}=y;
    [dummy,sample_overallInfected2(ix),sample_overallMortality2(ix),dummy]=computeFinalSize_generalized(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
end
sample_overallInfected=[sample_overallInfected1 sample_overallInfected2];
sample_overallMortality=[sample_overallMortality1 sample_overallMortality2];

save([fname,'_',vaccineRange,'_R=',num2str(R0)]);
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


function [c,ceq] = mycon(x,infectedContraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac)
w=1;
[dummy,overallInfected,dummy,dummy]=computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

c = [];     % Compute nonlinear inequalities at x.
ceq = infectedContraint-overallInfected;   % Compute nonlinear equalities at x.
return