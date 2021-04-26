susceptibilityFactor=1;maxPrct=100;VcPrct=73.2;betaVac=0.2;nuVac=1;effVac=1-0.05/betaVac;
% Default scenario
tic
computeParetoFront_sqp(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac);
toc

% All-or-none
tic
betaVacParm=0;nuVacParm=0.8; computeParetoFront_sqp(susceptibilityFactor,maxPrct,VcPrct,betaVacParm,nuVacParm,effVac);
toc

% Hesistancy
for maxPrctParm=[80 90]
    computeParetoFront_sqp(susceptibilityFactor,maxPrctParm,VcPrct,betaVac,nuVac,effVac);
end

% Efficacy
for betaVacParm=[0.15 0.25]%0.1 0.15 0.25 0.3]
    effVacParm=1-0.05/betaVacParm;
    computeParetoFront_sqp(susceptibilityFactor,maxPrct,VcPrct,betaVacParm,nuVac,effVac);
end

% Susceptibility
for susceptibilityFactorParm=[1.5 2]
    computeParetoFront_sqp(susceptibilityFactorParm,maxPrct,VcPrct,betaVac,nuVac,effVac);
end

% Vaccine coverage
for VcPrctParm=[66 80]
    computeParetoFront_sqp(susceptibilityFactor,maxPrct,VcPrctParm,betaVac,nuVac,effVac);
end


produceGraphs;
