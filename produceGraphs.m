addpath('./AuxilaryFunctions/');
collectData=false;

%% Produce graphs for Figure 1 & related SM materials (Vaccination required for herd immunity)
display('Producing graphs for Figure 1 (Vaccination required for herd immunity)');
cd('./ComputeThreshold')
produceVThresholdGraphs(collectData);
close all;
cd('..')

%% Produce graphs for Figure 2 & related SM materials (Critical reproduction numbers)
display('Producing graphs for Figure 2 (Critical reproduction numbers for various countries)');
cd('./Reff4differentCountries')
produceAllGraphs;
close all;
cd('..')

%% Produce graphs for Figure 3 & related SM materials (Final size example)
display('Producing graphs for Figure 3 (Final size of epidemic)');
cd('./US_finalSize_Example')
allOrNone=true;produceGraph_US_finalSize(allOrNone);
allOrNone=false;produceGraph_US_finalSize(allOrNone);
close all;
cd('..')

%% Produce graphs for Figure 4 & related SM materials (Pareto front)
display('Producing graphs for Figure 4 (Pareto front)');
cd('./ParetoFront')
produceAllGraphs(collectData);
close all;
cd('..')


%% Produce graphs for Figure 5 & related SM materials (Pareto front)
display('Producing graphs for Figure 5 (Impact of change in reproduction number)');
cd('./OutcomesAsFunctionOfR_kids')
produceAllGraphs(collectData);
close all;
cd('..')