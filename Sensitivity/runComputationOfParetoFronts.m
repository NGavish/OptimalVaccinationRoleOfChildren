[status,msg]=mkdir('./data');

susceptibilityFactor=1;maxPrct=100;VcPrct=73.2;betaVac=0.1;nuVac=1;effVac=1-0.05/betaVac;recoveredprct=0;
vaccineRangeVec='All';R0=6;
for ix=11:40
        %% Default scenario
        display('Default')
        computeParetoFrontRnd(susceptibilityFactor,maxPrct,VcPrct,betaVac,nuVac,effVac,vaccineRangeVec,R0,recoveredprct,ix);
end