function computeParetoFrontUSexample(R0,VcPrct,betaVac,nuVac,effVac,equivFromScratch)
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
addpath('../Codes/AllCountries/AuxilaryFunctions/');

%% Parameters
infected_nv_prct=0;infected_v_prct=0;
fname=['./data/dataParetoFront_USexample'];
maxPrct=95;

%% Select Country
country="USA";
% vaccineRangeVec={'above20','above16','above12','above6','All'};
%% Prepare data
[uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters_US_Example(VcPrct,betaVac,effVac,maxPrct);

allKidsAllocation=0*uniformAllocation;allKidsAllocation(1)=sum(uniformAllocation.*Ni)/Ni(1)
susceptibilityFactor=1;
susceptibilityProfile=ones(9,1);
if susceptibilityFactor==1.5
    susceptibilityProfile(1:2)=[1 2];
elseif susceptibilityFactor==2
    susceptibilityProfile(1:2)=[2 2];
end
Cij=diag(susceptibilityProfile)*Cij;M2=diag(Ni)*Cij*diag(1./Ni);[V,d]=eig(M2);Cij=Cij/max(d(:));

%% Define optimization problem
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',5e4,'ConstraintTolerance',1e-6,'StepTolerance',1e-10,'Display','none');%,'Algorithm','sqp','Display','none');%,'Display','iter');     
problem.solver = 'fmincon';
problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
problem.A=[];problem.B=[];
problem.lb=0*upperBound;
problem.ub=upperBound;

defineColors;
%% Run loop
[dummy,overallInfected_uniform,overallMortality_uniform]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
[dummy2,overallInfected_uniform2,overallMortality_uniform2]=computeFinalSize(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
[dummy,overallInfected_allKids,overallMortality_allKids]=computeFinalSize_generalized(allKidsAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
[dummy,overallInfected_noKids,overallMortality_noKids]=computeFinalSize_generalized(0*allKidsAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

Msample=5000;M=800;
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
[result,overallInfected_riskFocused,overallMortality_riskFocused,data_riskFocused]=computeFinalSize_generalized(riskFocusedAllocation,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
% spreadersAllocation=[                0
%                    0
%                    0
%    0.538280868634493
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000
%    1.000000000000000];
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
    y=randomizePoint(distribution{1+mod(ix,M)},Ni,Nadult,upperBound,vaccinesLeftToDistibute,0.05+0.5*mod(ix,3).^2);
    sampleAllocation1{ix}=y;
    [dummy,sample_overallInfected1(ix),sample_overallMortality1(ix),dummy]=computeFinalSize_generalized(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
    
    y=randomizePoint(uniformAllocation,Ni,Nadult,upperBound,vaccinesLeftToDistibute,1);
    sampleAllocation2{ix}=y;
    [dummy,sample_overallInfected2(ix),sample_overallMortality2(ix),dummy]=computeFinalSize_generalized(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
end
sample_overallInfected=[sample_overallInfected1 sample_overallInfected2];
sample_overallMortality=[sample_overallMortality1 sample_overallMortality2];

save([fname,'_R=',num2str(R0)]);
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