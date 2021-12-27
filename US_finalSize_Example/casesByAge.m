load("cases_by_age.mat")

bar([data,infectedDistributionUSA,infectedDistributionUS2]);
agDistLabels = {"0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"};
set(gca,'xticklabel',agDistLabels);xlabel('Age groups');
ylabel('Distribution of cases');

legend({'Actual case distribution','Prem et al.','Vespignani et al.'});

print -depsc ageDist
shg