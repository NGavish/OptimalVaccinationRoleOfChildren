function produceAllGraphs(collectData)
collectData=false;

allOrNone=false;susceptibilityFactor=1;
outcomeAsFunctionofR(collectData,allOrNone,susceptibilityFactor) % default

allOrNone=true;
outcomeAsFunctionofR(collectData,allOrNone,susceptibilityFactor) % all-or-none

OutcomeAsFunctionOfR_susceptibility15(collectData);
OutcomeAsFunctionOfR_susceptibility20(collectData);
return