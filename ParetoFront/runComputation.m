susceptibilityFactor=1;maxPrct=100;VcPrct=73.2;betaVac=0.2;nuVac=1;effVac=1-0.05/betaVac;

% Default scenario
%computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'All');
%computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above10');
%computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above20');
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above12');
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,'above16');
return
% All-or-none
betaVacParm=0;nuVacParm=0.8; computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacParm,nuVacParm,effVac);
betaVacParm=0;nuVacParm=0.8; computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacParm,nuVacParm,effVac,'above10');
betaVacParm=0;nuVacParm=0.8; computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacParm,nuVacParm,effVac,'above20');

% Hesistancy
for maxPrctParm=[80 90]
    computeParetoFront(susceptibilityFactor,maxPrctParm,VcPrct,betaVac,nuVac,effVac);
end

% Efficacy
for betaVacParm=[0.15 0.25]%0.1 0.15 0.25 0.3]
    effVacParm=1-0.05/betaVacParm;
    computeParetoFront(susceptibilityFactor,maxPrct,VcPrct,betaVacParm,nuVac,effVac);
end

% Susceptibility
for susceptibilityFactorParm=[1.5 2]
    computeParetoFront(susceptibilityFactorParm,maxPrct,VcPrct,betaVac,nuVac,effVac);
end

% Vaccine coverage
for VcPrctParm=[66 80]
    computeParetoFront(susceptibilityFactor,maxPrct,VcPrctParm,betaVac,nuVac,effVac);
end

% Preexisting immunity
VcPrct50=50;
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct50,betaVac,nuVac,effVac,'All',3,20);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct50,betaVac,nuVac,effVac,'above20',3,20);
computeParetoFront(susceptibilityFactor,maxPrct,VcPrct50,betaVac,nuVac,effVac,'above10',3,20);


