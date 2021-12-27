function produceAllGraphs(collectData)
collectData=false;

allOrNone=false;susceptibilityFactor=1;
outcomeAsFunctionofR(collectData,allOrNone,1.5) % susceptibilityFactor=1.5
outcomeAsFunctionofR(collectData,allOrNone,2) % susceptibilityFactor=2

allOrNone=true;
outcomeAsFunctionofR(collectData,allOrNone,susceptibilityFactor) % all-or-none

return