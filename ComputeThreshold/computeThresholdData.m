function computeThresholdData(betaVac,susceptibilityFactor,recoveredprct)
%% Parameters
fname=['Vthreshold_betaVac',num2str(100*betaVac),'_susceptibilityFactor',num2str(susceptibilityFactor(1)*10),num2str(susceptibilityFactor(2)*10),'recoveredPrct_',num2str(10*recoveredprct)];

% Set efficacy in preventing disease so that overall protection in 95%
effVac=1-0.05/betaVac;nuVac=1;

% Set no active infected
infected_nv_prct=0;infected_v_prct=0;

% Maximal vaccination per age group
maxPrct=100;

%% Select Country
country="USA";
countryData=load(join(['../countryData/',country,'_data.mat'],''));

vaccineRangeVec={'above20','above10','All','above16','above12','above6'};
R0=1;

%% Define optimization problem
problem.options = optimoptions('fmincon','MaxFunctionEvaluations',1e4,'Display','none','Algorithm','sqp');%,'Display','iter');
problem.solver = 'fmincon';

problem.A=[];problem.B=[];

[uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,'All',recoveredprct,infected_nv_prct,infected_v_prct,1,betaVac,effVac,maxPrct,true);

VcPrctVec=linspace(0,100*sum(Ni)./sum(Ni(3:9)),150);
tic
sampleSize=100;


for kx=1:numel(VcPrctVec)
    for jx=1:numel(vaccineRangeVec)
        VcPrct=VcPrctVec(kx);
        
        
        vaccineRange=vaccineRangeVec{jx};
        
        
        %% Prepare data
        [uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,VcPrct,betaVac,effVac,maxPrct,true);
        
        % Take into account susceptibility 
        susceptibilityProfile=ones(1,9);susceptibilityProfile(1:2)=susceptibilityFactor;
        Cij=diag(susceptibilityProfile)*Cij;M2=diag(Ni)*Cij*diag(1./Ni);[V,d]=eig(M2);Cij =Cij/max(d(:));

        Aeq=Nadult';Beq=vaccinesLeftToDistibute;
        lb=0*upperBound;
        ub=upperBound;
        
        % Start scan from R=1 to R=4
        w=0;
        y=fmincon(@(x) computeMaxEig(x,adultAges,IFR,1,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,nuVac,effVac,1,1),uniformAllocation,problem.A,problem.B,Aeq,Beq,lb,ub,[],problem.options);
        distribution{jx,kx}=y;
        result=computeMaxEig(y,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,nuVac,effVac,1,1);
        Threshold(jx,kx)=1./result;

        %% Verify convergance to global minimum
        parfor ix=1:sampleSize
                z=randomizePoint(y,Ni,Ni(adultAges),upperBound,vaccinesLeftToDistibute,(ix>1).*(0.1+0.5*mod(ix,2)));
                z=fmincon(@(x) computeMaxEig(z,adultAges,IFR,1,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,nuVac,effVac,1,1),uniformAllocation,problem.A,problem.B,Aeq,Beq,lb,ub,[],problem.options);
                result=computeMaxEig(z,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,nuVac,effVac,1,1);
                if abs(result-1)<1e-15 & max(distribution{jx,kx}-z)>0
                    keyboard % Since we had not enocountered any problems with convergence to a non-global minimum we stop the run to examine any such case
                end
                if result>1+1e-15
                    keyboard % Since we had not enocountered any problems with convergence to a non-global minimum we stop the run to examine any such case
                end
        end
    end
    [uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,'All',recoveredprct,infected_nv_prct,infected_v_prct,VcPrctVec(kx),betaVac,effVac,maxPrct,true);
    VcGnrl(kx)=100*vaccinesLeftToDistibute/sum(Ni);
    
end


save(['data',fname]);

return


