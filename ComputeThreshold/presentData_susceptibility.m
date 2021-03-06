function presentData_susceptibility
% This function prepares SI Figure -
% "Effect of change in childrens' susceptibility on vaccination coverage required for herd immunity."
close all;
defineColors;
tiledlayout(2,2, 'TileSpacing', 'compact','Padding','none');
susceptiblityTxt={'1020','2020'};


for jx_idx=1:2
    % Load computed values
    load(['./data/dataVthreshold_betaVac10_susceptibilityFactor',susceptiblityTxt{jx_idx},'recoveredPrct_0']);
    jx=jx_idx;
    axVec{jx}=nexttile(jx);

    % Present graph of vaccination coverage required for herd immunity
    prct20plus=100*sum(Ni(3:9))/sum(Ni);a=min(find(VcGnrl>prct20plus));
    plot(Threshold(1,1:a),VcGnrl(1:a),'-.','Color',blue,'linewidth',1.5);hold on;
    Rcritical20=Threshold(1,a);
    prct10plus=100*sum(Ni(2:9))/sum(Ni);a=min(find(VcGnrl>prct10plus));
    plot(Threshold(2,1:a),VcGnrl(1:a),'-','Color',orange,'linewidth',1.5);
    Rcritical10=Threshold(2,a);
    plot(Threshold(3,:),VcGnrl,'--','Color',green,'linewidth',1.5);

    % Add graph of vaccination coverage required for herd immunity in the homogenous case
    s=max(1-(1-betaVac)*VcGnrl/100-recoveredprct/100,0);
    plot(1./s,VcGnrl,':','linewidth',2)
    xlim([1 1/betaVac]);ylim([0 100]);

    Rthreshold=Threshold(1,end);

    % Handle legend, left and right axis
    ylabel('V_{threshold} (% of all population)')
    xlabel('R_0');
    ytickformat( 'percentage');
    legend('Allocation to ages 20 and older','Allocation to ages 10 and older','Allocation to all ages','Homogeneous Threshold','location','best','autoupdate','off');
    yyaxis right;plot(1./s,VcPrctVec,'k:');
    ylim(100*[0 sum(Ni)./sum(Ni(3:9))]);
    ylabel('V_{threshold} (% of adult population)')
    ax = gca;
    ax.YAxis(2).Color = 'k';
    ax.YAxis(2).TickValues=[0:20:120];grid on;
    ytickformat('percentage');

    axVec{jx+2}=nexttile(jx+2);
    % Prepare graph of present the optimal allocations that achieve herd
    % immunity when all population is eligibile for vaccination
    R=Threshold(3,:);
    Mat=cell2mat(distribution(3,:));
    hndl=area(Threshold(3,:),100*(Mat.*Ni)'/sum(Ni)) ;   %title([vaccineRange,' - along w=0']);
    xlabel('R_0');ylabel('Vaccine supply (percent of population)');ytickformat( 'percentage');
    dataMat=100*(Mat.*Ni)'/sum(Ni);

    % Add age group text labels to graph
    a=134;
    textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
    currY=-1.5;
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

    % Handle left and right axis 
    grid on
    yyaxis left;xlim([0 1./betaVac]);
    ylim([ 0 100 ]);

    yyaxis right;
    xlim([0 1./betaVac]);
    ylim(100*[0 sum(Ni)./sum(Ni(3:9))]);
    ax = gca;
    ax.YAxis(2).Color = 'k';
    ylabel('V_{threshold} (% of adult population)');
    ytickformat( 'percentage');
    ax.YAxis(2).TickValues=[0:20:120];


    xlim([1 1/betaVac]);
end

% Superimpose baseline graph data
gray=0.8*[1 1 1];
for jx_idx=1:2

    load('./data/dataVthreshold_betaVac10_susceptibilityFactor1010recoveredPrct_0');hold on;
    jx=jx_idx;
    axVec{jx}=nexttile(jx);

    prct20plus=100*sum(Ni(3:9))/sum(Ni);a=min(find(VcGnrl>prct20plus));
    yyaxis left;plot(Threshold(1,1:a),VcGnrl(1:a),'-.','Color',gray,'linewidth',1.5);hold on;
    Rcritical20=Threshold(1,a);ylim([0 100]);
    prct10plus=100*sum(Ni(2:9))/sum(Ni);a=min(find(VcGnrl>prct10plus));
    yyaxis left;plot(Threshold(2,1:a),VcGnrl(1:a),'-','Color',gray,'linewidth',1.5);
    %Rcritical10=Threshold(2,a);
    plot(Threshold(3,:),VcGnrl,'--','Color',gray,'linewidth',1.5);

    s=max(1-(1-betaVac)*VcGnrl/100-recoveredprct/100,0);

    Rthreshold=Threshold(1,end);
end

% Link axis and add text label to graphs
linkaxes([ axVec{3} axVec{4} ],'xy');
Ages={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
lg = legend(axVec{4},Ages,'fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
lg.Layout.Tile='south'

% Add legend & subplot reference text
set(gcf,'Position',[520 426 942 486]);
title(axVec{1},'Adolescents (ages 10-19) are equally susceptible as adults');
title(axVec{2},'Children (ages 0-19) are equally susceptible as adults');
text(axVec{1},0.05,0.9,'A','units','normalized','fontweight','bold','fontsize',12);
text(axVec{2},0.05,0.9,'C','units','normalized','fontweight','bold','fontsize',12);
text(axVec{3},0.05,0.9,'B','units','normalized','fontweight','bold','fontsize',12);
text(axVec{4},0.05,0.9,'D','units','normalized','fontweight','bold','fontsize',12);

% Save results
printGraph('../graphs/Vthreshold_with_susceptibility');