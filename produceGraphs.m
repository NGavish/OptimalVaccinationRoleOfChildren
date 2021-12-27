addpath('./AuxilaryFunctions/');
collectData=false;

%% Produce graphs for Figure 1 & related SM materials (Final size example)
display('Producing graphs for Figure 1 (Final size of epidemic)');
cd('./US_finalSize_Example')
allOrNone=true;produceGraph_US_finalSize(allOrNone);
allOrNone=false;produceGraph_US_finalSize(allOrNone);
close all;
cd('..')

%% Produce graphs for Figure 2 & related SM materials (Vaccination required for herd immunity)
display('Producing graphs for Figure 2 (Vaccination required for herd immunity)');
cd('./ComputeThreshold');
produceVThresholdGraphs(collectData);

close all;
cd('..')

%% Produce graphs for Figure 3 & related SM materials (Critical reproduction numbers)
display('Producing graphs for Figure 2 (Critical reproduction numbers for various countries)');
cd('./Reff4differentCountries')
produceAllGraphs;
close all;
cd('..')

%% Produce graphs for Figure 4 & related SM materials (Pareto front)
display('Producing graphs for Figure 4 (Impact of change in reproduction number)');
cd('./OutcomesAsFunctionOfR_kids')
produceAllGraphs(collectData);
close all;
cd('..')

%% Produce graphs for Figure 5-7 & related SM materials (Pareto front)
display('Producing graphs for Figures 5-7 (Pareto front)');
cd('./ParetoFront')
produceAllGraphs(collectData);
close all;
cd('..')

%% Produce graphs for Figure 11
display('Pareto front - vaccination of an additional 300,000 people in Israel');
cd('./IsraelTestCase')
produceAllGraphs(collectData);
close all;
cd('..')

%% Produce graphs for Figure 12
display('Producing graphs for SI Figure 12 (Effect of perturbations in contact matrix)');
cd('./Sensitivity')
presentData(collectData);
close all;
cd('..')

