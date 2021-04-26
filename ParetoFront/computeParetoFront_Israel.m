function computeParetoFront_Israel(vaccineRange)
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
addpath('../Codes/AllCountries/AuxilaryFunctions/');
%vaccineRange='above12';
susceptibilityFactor=1;
maxPrct=99.99;
betaVac=0.2;
nuVac=1;
effVac=0.75;
R0=3.5;

%% Parameters
infected_nv_prct=0;infected_v_prct=0;
fname=['./data/dataParetoFront_Israel_sus075'];

%% Select Country
country="ISRAEL";
vaccineRangeVec={'above20','above16','above12','above6','All'};
% vaccineRange='above16';
% vaccineRange='All';
%% Prepare data
%[uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,0.99*VcPrct,betaVac,effVac,maxPrct,true);
detRates=[0.4*ones(1,9),0.55,0.70,0.70,0.70,0.70,0.85,1]';

agDist=[0.1010 0.0952 0.0854 0.0790 0.0718 0.0681 0.0661 0.0642 0.0618 0.0562 0.0473 0.0438 0.0414 0.0392 0.0316 0.0480]';
N=agDist*9.3e6;%Nadult=N(4:16);
Ni=N;
IFR=[0,0,0,0,0,0,0,0,0,0.001,0.003,0.006,0.012,0.024,0.049,0.164]'./detRates;

%%
useCovidMatrix=true;
load IsraelContactMatrices.mat
if useCovidMatrix
    contactMatrix=CovidCij;
else
    contactMatrix=preCovidCij;
end

susVec=[0.5 0.5 0.5 ones(1,13)];susVec=[0.75 0.75 0.75 ones(1,13)];
infVec=[0.75 0.75 0.75 ones(1,13)];infVec=ones(size(infVec));
contactMatrix=diag(susVec)*contactMatrix*diag(infVec);
Mij=diag(Ni)*contactMatrix*diag(1./Ni);
[V,d]=eig(Mij);contactMatrix=contactMatrix/max(d(:));
%%
Cij=contactMatrix;
rKnown=[49967 68429 76767 97743 85793 67837 59189 52699 48466 45101 36691 29868 23161 17106 11127 15491]'./Ni;
r=rKnown./detRates;
runKnown=r-rKnown;
v_wasted=runKnown./(1-rKnown);
v=[320    156    350 326475 493972 451939 446919 445365 433317 437108 384319 341760 333801 322872 295377 448257]'./Ni;
v_actual=min(1,(1-v_wasted).*v);
v=v_actual;
infected0_nv=0*[0.045, 0.062, 0.079 ,0.127, 0.112 ,0.086 ,0.074, 0.067, 0.065, 0.063 ,0.053, 0.045 ,0.037, 0.029 ,0 ,0.035]'./Ni;
infected0_v=0*[0, 0, 0 ,0.127, 0.112 ,0.086 ,0.074, 0.067, 0.065, 0.063 ,0.053, 0.045 ,0.037, 0.029 ,0 ,0.035]'./Ni;

susceptibilityProfile=ones(16,1);susceptibilityProfile(1:3)=susceptibilityFactor;
Cij=diag(susceptibilityProfile.*Ni)*Cij*diag(1./Ni);[V,d]=eig(Cij);Cij=Cij/max(d(:));

upperboundFilter=ones(size(N));upperboundFilter(end)=0;
switch vaccineRange
    case 'above16'
        Nadult=N(4:14);
        upperboundFilter(4)=0.8;
        upperboundFilter(1:3)=0;
        sIdx=4:14;
        adultAges=4:14;
    case 'above12'
        Nadult=N(3:14);
        upperboundFilter(3)=0.6;
        upperboundFilter(1:2)=0;
        sIdx=3:14;
        adultAges=sIdx;
    case 'above6'
        Nadult=N(2:14);
        upperboundFilter(2)=0.8;
        upperboundFilter(1)=0;
        sIdx=2:14;
    case 'All'
        Nadult=N(1:14);
        sIdx=1:14;
end
adultAges=sIdx;
upperBound=upperboundFilter(adultAges);
vaccinesLeftToDistibute=5e5;

%% Define optimization problem
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',5e4,'ConstraintTolerance',1e-6,'StepTolerance',1e-10,'Display','none');%,'Algorithm','sqp','Display','none');%,'Display','iter');
problem.solver = 'fmincon';
problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
problem.A=[];problem.B=[];
problem.lb=0*upperBound;
problem.ub=upperBound;

%% Uniform
s=max(0,1-r-v-infected0_v-infected0_nv);
uniformAllocation =[s(sIdx).*Nadult/sum(s(sIdx).*Nadult)*problem.Beq]./Nadult;
spreadersAllocation=uniformAllocation;
idxVec(1:numel(uniformAllocation))=true;%idxVec(12:13)=false;
idx=find(uniformAllocation>problem.ub);

while numel(idx)>0
    vacExtra=(uniformAllocation(idx)-problem.ub(idx))'*Nadult(idx);
    uniformAllocation(idx)=max(problem.ub(idx)-1e-15,0);
    idxVec(idx)=false;
    n_idx=find(idxVec);
    uniformAllocation(n_idx)=uniformAllocation(n_idx)+(s(n_idx).*Nadult(n_idx)/sum(s(n_idx).*Nadult(n_idx))*vacExtra)./Nadult(n_idx);
    idx=find(uniformAllocation>problem.ub);
end
%%

defineColors;
%% Run loop
uniformAllocation
[dummy,overallInfected_uniform,overallMortality_uniform,data_uniform]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

Msample=5000;M=400;
%% First solve for infection minimizing allocation
w=0;
[result,overallInfected_spreaders,overallMortality_spreaders,data_spreaders]=computeFinalSize_generalized(spreadersAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
[x,fval,exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
[dummy,overallInfected(1),overallMortality(1),data{1}]=computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
distribution{1}=x;
parfor ix=1:Msample
    y=randomizePoint(distribution{1},Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
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
    y=randomizePoint(distribution{1},Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
    [aux,overallMortality_sample(kx),exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac),problem.options);
    if exitflag<=0
        overallMortality_sample(kx)=Inf;
    end
    %[dummy,overallInfected_sample(kx),overallMortality_sample(kx),dummy]=computeFinalSize_generalized(aux,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
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
[result,overallInfected_riskFocused,overallMortality_riskFocused,data_riskFocused]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
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
        y=randomizePoint(x,Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.01);
        [aux,overallMortality_sample(kx),exitflag,output]=fmincon(@(x)computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,@(x)mycon(x,infectedConstraint,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac),problem.options);
        if exitflag<=0
            overallMortality_sample(kx)=Inf;
        end
        sampleAllocation1{kx}=aux;
    end
    
    [minResult,minIdx]=min(overallMortality_sample(1:10));
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
    y=randomizePoint(distribution{1+mod(ix,M)},Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.01);
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