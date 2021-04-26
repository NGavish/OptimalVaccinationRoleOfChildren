function produceVThresholdGraphs(collectData)
if collectData
    display('Collecting data');
    % Default graph
    susceptibilityFactor=1;recoveredPrct=0;betaVac=0.2;
    computeThresholdData(betaVac,susceptibilityFactor,recoveredPrct);
    
    % Efficacy (SM)
    susceptibilityFactor=1;recoveredPrct=0;betaVac=0.2;
    for betaVac=[0.15 0.25]
        computeThresholdData(betaVac,susceptibilityFactor,recoveredPrct)
    end
    
    % Susceptibility (SM)
    susceptibilityFactor=1;recoveredPrct=0;betaVac=0.2;
    for susceptibilityFactor=[1.5 2]
        computeThresholdData(betaVac,susceptibilityFactor,recoveredPrct)
    end
    
    % Recoverd (SM)
    susceptibilityFactor=1;recoveredPrct=0;betaVac=0.2;
    
    for recoveredPrct=[10 20]
        computeThresholdData(betaVac,susceptibilityFactor,recoveredPrct)
    end
end

% Default graph
susceptibilityFactor=1;recoveredPrct=0;betaVac=0.2;
presentDataDefault(betaVac,susceptibilityFactor,recoveredPrct)

% Default graph
presentDataDefault1612(betaVac,susceptibilityFactor,recoveredPrct);shg
return
% Recovered
presentData_recovered;

% Efficacy
presentData_efficacy;

% Susceptibility
presentData_susceptibility;

return

function presentDataDefault(betaVac,susceptibilityFactor,recoveredprct)
defineColors;

fname=['Vthreshold_betaVac',num2str(100*betaVac),'_susceptibilityFactor',num2str(susceptibilityFactor*10),'recoveredPrct_',num2str(10*recoveredprct)];
load(['data',fname]);

tiledlayout(3,2, 'TileSpacing', 'compact','Padding','none'); axVec{1}=nexttile(1);

prct20plus=100*sum(Ni(3:9))/sum(Ni);a=min(find(VcGnrl>prct20plus));
plot(Threshold(1,1:a),VcGnrl(1:a),'-.','Color',blue,'linewidth',1.5);hold on;
Rcritical20=Threshold(1,a);
prct10plus=100*sum(Ni(2:9))/sum(Ni);a=min(find(VcGnrl>prct10plus));
plot(Threshold(2,1:a),VcGnrl(1:a),'-','Color',orange,'linewidth',1.5);
Rcritical10=Threshold(2,a);
plot(Threshold(3,:),VcGnrl(:),'--','Color',green,'linewidth',1.5);

s=max(1-(1-betaVac)*VcGnrl/100-recoveredprct/100,0);
plot(1./s,VcGnrl,'k:','linewidth',2)
ylim([0 100]);

Rthreshold=Threshold(1,end);

ylabel('V_{threshold} (% of all population)')
xlabel('R_0');
ytickformat( 'percentage');
legend('Eligibilty - ages 20 and older','Eligibilty - ages 10 and older','Eligibilty - all ages','Homogeneous Threshold','location','best','autoupdate','off');
yyaxis right;plot(1./s,VcPrctVec,'k:');ylim(100*[0 sum(Ni)./sum(Ni(3:9))]);
ylabel('V_{threshold} (% of adult population)')
ax = gca;
ax.YAxis(2).Color = 'k';

ax1 = gca; % position of first axes

grid on;
set(gcf,'Position',[520 246 659 551]);
ax.YAxis(2).TickValues=[0:20:100];
xlim([1 5]);
ytickformat('percentage');

axVec{3}=nexttile(3);
R=Threshold(3,:);
Mat=cell2mat(distribution(3,:));
hndl=area(Threshold(3,:),100*(Mat.*Ni)'/sum(Ni)) ;   %title([vaccineRange,' - along w=0']);
xlabel('R_0');ylabel('Vaccine supply (percent of population)');ytickformat( 'percentage');
dataMat=100*(Mat.*Ni)'/sum(Ni);
a=130;
textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
currY=0;
for ix=1:6
    currY=currY+dataMat(a,ix);
    if dataMat(a,ix)>sum(dataMat(a,:))/50
        text(R(a),currY,textAges{ix},'verticalAlign','top','color','w','fontweight','bold')
    end
end

text(R(a),currY+11,textAges{7},'verticalAlign','top','color','w','fontweight','bold')

for ix=1:9
    hndl(ix).FaceColor=lineColors(ix,:);
end

grid on
ylim([0 100]);
yyaxis right;plot(1./s,VcPrctVec,'w');ylim(100*[0 sum(Ni)./sum(Ni(3:9))]);
ax = gca;
ax.YAxis(2).Color = 'k';
ylabel('V_{threshold} (% of adult population)');
ytickformat( 'percentage');
ax.YAxis(2).TickValues=[0:20:100];
linkaxes([axVec{1} axVec{3}],'xy');
text(axVec{1},0.05,0.9,'A','units','normalized','fontsize',13);
text(axVec{3},0.05,0.9,'B','units','normalized','fontsize',13);

% Create textbox
annotation('textbox',...
    [0.379578696678468 0.590634414029361 0.0289351851851849 0.0374753451676533],...
    'String',{'80+'},...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textarrow
annotation('textarrow',[0.462165598831 0.496887821053222],...
    [0.61717316193974 0.600421682649799]);

% Create textarrow
annotation('textarrow',[0.412089726296857 0.446811948519079],...
    [0.606283869743732 0.589532390453791]);

% Create textbox
annotation('textbox',...
    [0.414506407013994 0.596634414029361 0.041666666666666 0.0374753451676533],...
    'String','70-79',...
    'FitBoxToText','off',...
    'EdgeColor','none');

xlim([1 5]);

% Take care of xlabel 
set(axVec{1},'xtick',sort([1:0.5:1.5 Rcritical20 2.5 3 Rcritical10 4 4.5 5]));
xtickList=axVec{1}.XTick;
xtickLabelList=axVec{1}.XTickLabel;
xtickLabelList{3}='          R^{20+}_{critical}';
xtickLabelList{6}='     R^{10+}_{critical}';
set(axVec{1},'xtick',xtickList,'xticklabel',xtickLabelList);axVec{1}.XTickLabelRotation=0;
set(axVec{3},'xtick',xtickList,'xticklabel',xtickLabelList);axVec{3}.XTickLabelRotation=0;

printGraph(['../graphs/',fname]);
return

function presentDataDefault1612(betaVac,susceptibilityFactor,recoveredprct)
defineColors;

fname=['Vthreshold_betaVac',num2str(100*betaVac),'_susceptibilityFactor',num2str(susceptibilityFactor*10),'recoveredPrct_',num2str(10*recoveredprct)];
load(['data',fname]);

tiledlayout(2,2, 'TileSpacing', 'compact','Padding','none'); axVec{1}=nexttile([1 2]);
gray=0.7*[1 1 1];
% above 20
prct20plus=100*sum(Ni(3:9))/sum(Ni);a=min(find(VcGnrl>prct20plus));
plot(Threshold(1,1:a),VcGnrl(1:a),'-.','Color',gray,'linewidth',1.5);hold on;
Rcritical20=Threshold(1,a);
% above 16
a=120
plot(Threshold(4,1:a),VcGnrl(1:a),'-.','Color',blue,'linewidth',1.5);
Rcritical16=Threshold(4,a);
% above 12
a=129
plot(Threshold(5,1:a),VcGnrl(1:a),'-','Color',orange,'linewidth',1.5);
Rcritical12=Threshold(5,a);

% above 10
prct10plus=100*sum(Ni(2:9))/sum(Ni);a=min(find(VcGnrl>prct10plus));
plot(Threshold(2,1:a),VcGnrl(1:a),'-','Color',gray,'linewidth',1.5);
Rcritical10=Threshold(2,a);
% 'All'
a=150
plot(Threshold(3,1:a),VcGnrl(1:a),'--','Color',gray,'linewidth',1.5);

ylim([0 100]);

Rthreshold=Threshold(1,end);

ylabel('V_{threshold} (% of all population)')
xlabel('R_0');
ytickformat( 'percentage');
legend('Eligibilty - ages 20 and older','Eligibilty - ages 16 and older','Eligibilty - ages 12 and older','Eligibilty - ages 10 and older','Eligibilty - all ages','location','best','autoupdate','off');
% plot(1./s,VcPrctVec,'k:');
yyaxis right;ylim(100*[0 sum(Ni)./sum(Ni(3:9))]);
ylabel('V_{threshold} (% of adult population)')
ax = gca;
ax.YAxis(2).Color = 'k';

ax1 = gca; % position of first axes

grid on;
set(gcf,'Position',[520 246 659 551]);
ax.YAxis(2).TickValues=[0:20:100];
xlim([1 5]);
ytickformat('percentage');

% Take care of xlabel 
set(axVec{1},'xtick',sort([1:0.5:1.5 Rcritical20 Rcritical16 Rcritical12  Rcritical10 4 4.5 5]));
xtickList=axVec{1}.XTick;
xtickLabelList=axVec{1}.XTickLabel;
xtickLabelList{3}='     R^{20+}_{critical}';
xtickLabelList{4}='     R^{16+}_{critical}';
xtickLabelList{5}='     R^{12+}_{critical}';
xtickLabelList{6}='       R^{10+}_{critical}';
set(axVec{1},'xtick',xtickList,'xticklabel',xtickLabelList);axVec{1}.XTickLabelRotation=0;
% set(axVec{3},'xtick',xtickList,'xticklabel',xtickLabelList);axVec{3}.XTickLabelRotation=0;


axVec{3}=nexttile(3);
R=Threshold(4,:);
Mat=cell2mat(distribution(4,:));
areadata=100*(Mat.*Ni(2:9))'/sum(Ni(2:9));

title('Allocation - vaccine eligibililty age 16 and older');
hndl=area(Threshold(4,:),areadata) ;   %title([vaccineRange,' - along w=0']);
hndl(1).FaceColor = 'k';
xlabel('R_0');ylabel('Vaccine supply (percent of population)');ytickformat( 'percentage');
xlim([1 Rcritical16]);
set(axVec{3},'xtick',sort([1:0.5:1.5 Rcritical20 Rcritical16 Rcritical12]));
xtickLabelList{3}='     R^{20+}_{critical}';
xtickLabelList{4}='     R^{16+}_{critical}';
xtickLabelList{5}='     R^{12+}_{critical}';
set(axVec{3},'xtick',xtickList,'xticklabel',xtickLabelList);axVec{3}.XTickLabelRotation=0;

   for ix=1:8
            hndl(ix).FaceColor=lineColors(ix+1,:);
        end
%         hndl(2).FaceColor=blue;
        
        title('Vaccine eligibililty - ages 16 and older');
        
axVec{4}=nexttile(4);
R=Threshold(5,:);
Mat=cell2mat(distribution(5,:));
areadata=100*(Mat.*Ni(2:9))'/sum(Ni(2:9));
areadata=[areadata(:,1)/2 areadata(:,1)/2 areadata(:,2:end)];



hndl=area(Threshold(5,:),areadata) ;   %title([vaccineRange,' - along w=0']);
xlabel('R_0');ylabel('Vaccine supply (percent of population)');ytickformat( 'percentage');
xlim([1 Rcritical12]);
set(axVec{4},'xtick',sort([1:0.5:1.5 Rcritical20 Rcritical16 Rcritical12]));
xtickLabelList{3}='     R^{20+}_{critical}';
xtickLabelList{4}='     R^{16+}_{critical}';
xtickLabelList{5}='     R^{12+}_{critical}';
set(axVec{4},'xtick',xtickList,'xticklabel',xtickLabelList);axVec{4}.XTickLabelRotation=0;
title('Vaccine eligibililty - ages 12 and older');
        for ix=3:9
            hndl(ix).FaceColor=lineColors(ix,:);
        end
        hndl(1).FaceColor=gray;
        hndl(2).FaceColor=blue;
textAges={'12-15','16-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};

lg = legend(axVec{4},textAges,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
 lg.Layout.Tile='south'
% dataMat=100*(Mat.*Ni)'/sum(Ni);
printGraph(['../graphs/Vthreshold_ages1216']);

return