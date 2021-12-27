[status,msg]=mkdir('./data');

% Compute pareto front Figure 5-7 for R0=4,6,8
susceptibilityFactor=1;maxPrct=100;VcPrct=73.2;betaVac=0.1;nuVac=1;effVac=1-0.05/betaVac;recoveredprct=0;
vaccineRangeVec={'All','above20','above10'};
for R0=[4 6 8]
    for ix=1:3
        %% Default scenario
        display('Default')
        computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);
    end
end

% Compute pareto fronts for SI (R0=6)
R0=6;
ix=1;
for ix=1:3
    %% All-or-none
    display('All-or-none')
    betaVacAllorNone=0;nuVacAllorNone=0.9;
    %computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacAllorNone,nuVacAllorNone,effVac,vaccineRangeVec{ix},R0,recoveredprct);

    %% Susceptibility
    display('Susceptibility')
    susceptibilityFactor15=1.5;
    susceptibilityFactor20=2;
    computeParetoFront(susceptibilityFactor15,maxPrct,VcPrct,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);
    computeParetoFront(susceptibilityFactor20,maxPrct,VcPrct,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);

end
ix=1;

%% Hesitancy
display('Hesitancy')
maxPrct90=90;
maxPrct80=80;

computeParetoFront(susceptibilityFactor,maxPrct90,VcPrct,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);
computeParetoFront(susceptibilityFactor,maxPrct80,VcPrct,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);

%% Efficacy
display('Efficacy')
betaVacVE85=0.15;effVacVE85=1-0.05/betaVacVE85;  % 85% efficacy
betaVacVE95=0.05;effVacVE95=1-0.05/betaVacVE95;  % 95% efficacy
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacVE85,nuVac,effVacVE85,vaccineRangeVec{ix},R0,recoveredprct);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacVE95,nuVac,effVacVE85,vaccineRangeVec{ix},R0,recoveredprct);

% Coverage
display('Coverage')

computeParetoFront(susceptibilityFactor,maxPrct,VcPrct+5,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct-5,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);

%% Recovered
display('Recovered')
recoveredprct10=10
recoveredprct20=20
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct20);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct-recoveredprct10,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct10);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct-recoveredprct20,betaVac,nuVac,effVac,vaccineRangeVec{ix},R0,recoveredprct);

produceAllGraphs(false)