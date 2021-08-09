function outcomeAsFunctionofR(collectData,allOrNone,susceptibilityFactor)
collectData=false;
allOrNone=false;
susceptibilityFactor=1;
addpath('../Codes/AllCountries/AuxilaryFunctions/');
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
close all
%% Parameters

if allOrNone
    nuVac=0.8;betaVac=0;
    fname='all_or_none'
    a2=270
else
    nuVac=1;betaVac=0.2; 
    fname='default'
    a2=270
end

effVac=0.75;recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;
susceptibilityProfile=ones(1,9);susceptibilityProfile(1:2)=susceptibilityFactor;
switch susceptibilityFactor
    case 1.5
        a2=350
    case 2
        a2=450;
end
susFactorText=['susceptiblityFactor=',num2str(susceptibilityFactor*10)];
% The percent of 20+ population that have vaccines avaliable
maxPrct=100; % Maximal vaccination per age group

%% Select Country
country="USA";
countryData=load(join(['../countryData/',country,'_data.mat'],''));

VcPrctVec=73.2;
vaccineRangeVec={'All','above10','above20'};%'above16','above12','above6'
%% Run loop
if collectData
    sampleSize=5000;
    Rvec=unique([linspace(1,4,200) linspace(1.74,1.86,100) linspace(2.5,2.7,200) linspace(2.7,3.2,200) linspace(2.8,2.87,50) linspace(3,3.05,50) linspace(3,3.4,100)]);
    count=0
    M=numel(Rvec);
    for kx=1:numel(VcPrctVec)
        for jx=1:numel(vaccineRangeVec)
            tic
            VcPrct=VcPrctVec(kx);
            vaccineRange=vaccineRangeVec{jx};
            
            %% Prepare data
            [uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,VcPrct,betaVac,effVac,maxPrct,true);
            Cij=diag(susceptibilityProfile)*Cij;M2=diag(Ni)*Cij*diag(1./Ni);[V,d]=eig(M2);Cij =Cij/max(d(:));

            xw0=uniformAllocation;%spreadersAllocation;
            xw1=riskFocusedAllocation;
            
            %% Define optimization problem
            problem.options = optimoptions('fmincon','MaxFunctionEvaluations',1e4,'Display','none');%,'Algorithm','sqp');%,'Display','iter');
            problem.solver = 'fmincon';
            problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
            problem.A=[];problem.B=[];
            problem.lb=0*upperBound;
            problem.ub=upperBound;
            
            % Start low-res scan from R=1 to R=4
            w=0
            xw0=uniformAllocation;
            xw1=riskFocusedAllocation;
            for ix=1:M
                R0=Rvec(ix);
%                             [result,overallInfected_nv_uniform(ix),overallFatality_uniform(ix),data_uniform{ix}]=computeFinalSize(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                [result,overallInfected_nv_uniform(ix),overallFatality_uniform(ix),data_uniform{ix}]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
                [resultA,overallInfected_nv_w0(ix),overallFatality_w0(ix),data_w0{ix}]=computeFinalSize_generalized(xw0,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
                w=0;y=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),xw0,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                [resultB,overallInfected_nv_w0(ix),overallFatality_w0(ix),data_w0{ix}]=computeFinalSize_generalized(y,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
                if resultB<resultA
                    xw0=y;
                    count=count+1;
                    [ix count]
                end
                parfor jx=1:sampleSize
                    y=randomizePoint(xw0,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.1+0.3*mod(jx,2)));
                    w=0;[sample_result,sample_overallInfected_nv(jx),sample_overallFatality(jx),sample_data{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                    sample_distribution{jx}=y;
                    
                    y=randomizePoint(xw1,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.1+0.3*mod(jx,2)));
                    w=1;[sample_result(jx),sample_overallInfected_nv2(jx),sample_overallFatality2(jx),sample_data2{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                    sample_distribution2{jx}=y;
                end
                sample_overallInfected_nv(sampleSize+1:2*sampleSize)=sample_overallInfected_nv2;
                sample_overallFatality(sampleSize+1:2*sampleSize)=sample_overallFatality2;
                sample_distributionCombined=[sample_distribution sample_distribution2];
                sample_dataCombined=[sample_data sample_data2];
                
                minIdx=find(sample_overallInfected_nv>0);
                [minResultLTR,minIdx]=min(sample_overallInfected_nv(minIdx));
                overallInfected_nv_w0_LTR(ix)=sample_overallInfected_nv(minIdx);
                overallFatality_w0_LTR(ix)=sample_overallFatality(minIdx);
                data_w0_LTR{ix}=sample_dataCombined{minIdx};
                distribution_w0_LTR{ix}=sample_distributionCombined{minIdx};
                
                minIdx=find(sample_overallFatality>0);
                [minResultLTR,minIdx]=min(sample_overallFatality(minIdx));
                overallInfected_nv_w1_LTR(ix)=sample_overallInfected_nv(minIdx);
                overallFatality_w1_LTR(ix)=sample_overallFatality(minIdx);
                data_w1_LTR{ix}=sample_dataCombined{minIdx};
                distribution_w1_LTR{ix}=sample_distributionCombined{minIdx};
                
                xw0=distribution_w0_LTR{ix};
                xw1=distribution_w1_LTR{ix};
            end
            
            % Low-res scan from R=4 to R=1
            for ix=fliplr(1:M)
                R0=Rvec(ix);
                parfor jx=1:sampleSize
                    y=randomizePoint(xw0,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.1+0.3*mod(jx,2)));
                    w=0;[sample_result(jx),sample_overallInfected_nv(jx),sample_overallFatality(jx),sample_data{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                    sample_distribution{jx}=y;
                    
                    y=randomizePoint(xw1,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.1+0.3*mod(jx,2)));
                    w=1;[sample_result2(jx),sample_overallInfected_nv2(jx),sample_overallFatality2(jx),sample_data2{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                    sample_distribution2{jx}=y;
                end
                sample_overallInfected_nv(sampleSize+1:2*sampleSize)=sample_overallInfected_nv2;
                sample_overallFatality(sampleSize+1:2*sampleSize)=sample_overallFatality2;
                sample_distributionCombined=[sample_distribution sample_distribution2];
                sample_dataCombined=[sample_data sample_data2];
                
                minIdx=sample_overallInfected_nv>0;
                [minResultRTL(ix),minIdx]=min(sample_overallInfected_nv(minIdx));
                if overallInfected_nv_w0_LTR(ix)>minResultRTL(ix)
                    overallInfected_nv_w0(ix)=sample_overallInfected_nv(minIdx);
                    overallFatality_w0(ix)=sample_overallFatality(minIdx);
                    data_w0{ix}=sample_dataCombined{minIdx};
                    distribution_w0{ix}=sample_distributionCombined{minIdx};
                else
                    overallInfected_nv_w0(ix)=overallInfected_nv_w0_LTR(ix);
                    overallFatality_w0(ix)=overallFatality_w0_LTR(ix);
                    data_w0{ix}=data_w0_LTR{ix};
                    distribution_w0{ix}=distribution_w0_LTR{ix};
                end
                
                [minResultRTL(ix),minIdx]=min(sample_overallFatality);
                if overallFatality_w1_LTR(ix)>minResultRTL(ix)
                    overallInfected_nv_w1(ix)=sample_overallInfected_nv(minIdx);
                    overallFatality_w1(ix)=sample_overallFatality(minIdx);
                    data_w1{ix}=sample_dataCombined{minIdx};
                    distribution_w1{ix}=sample_distributionCombined{minIdx};
                else
                    overallInfected_nv_w1(ix)=overallInfected_nv_w1_LTR(ix);
                    overallFatality_w1(ix)=overallFatality_w1_LTR(ix);
                    data_w1{ix}=data_w1_LTR{ix};
                    distribution_w1{ix}=distribution_w1_LTR{ix};
                end
                xw0=distribution_w0{ix};
                xw1=distribution_w1{ix};
            end
            
            save(['../data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
            toc
        end
    end
end

lineStyle={'-','--','-.'}

for kx=1:numel(VcPrctVec)
    close all;
t=tiledlayout(3,2, 'TileSpacing','compact','Padding','none')
    VcPrct=VcPrctVec(kx);

    for ix=1:numel(vaccineRangeVec)
        
        vaccineRange=vaccineRangeVec{ix};
        data=load(['../data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
        VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
        Rvec=data.Rvec;
        %% Present Data
        green=[0.4660, 0.6740, 0.1880];orange=[0.8500, 0.3250, 0.0980];blue=[0, 0.4470, 0.7410];gray=0.5*[1 1 1];darkgray=0.4*[1 1 1];
        overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
        overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
        overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;
        
        R_threshold(ix)=Rvec(min(find(overallInfected_nv_w0>0.5)));
        
        axVec{1}=nexttile(1);plot(Rvec,overallInfected_nv_w0,'linewidth',2,'linestyle',lineStyle{ix});hold on;
        grid on;box on;ax = gca;ax.YRuler.Exponent = 0;
        xlabel('R_0');ylabel('Overall infected');
        ytickformat('%gM');%ylim([0 85])
        axVec{2}=nexttile(2);plot(Rvec,overallInfected_nv_w1,'linewidth',2,'linestyle',lineStyle{ix});hold on;
        grid on;box on;ax = gca;ax.YRuler.Exponent = 0;
        xlabel('R_0');ylabel('Overall infected');
        ytickformat('%gM');%ylim([0 85])
        
        box on;
        
        axVec{3}=nexttile(3);plot(Rvec,data.overallFatality_w0/1e6,'linewidth',2,'linestyle',lineStyle{ix});hold on;grid on;
        xlabel('R_0');ylabel('Mortality');
        ytickformat('%0.1gM');ax.YRuler.Exponent = 0;box on;
        axVec{4}=nexttile(4);plot(Rvec,data.overallFatality_w1/1e6,'linewidth',2,'linestyle',lineStyle{ix});hold on;grid on;
        xlabel('R_0');ylabel('Mortality');
        ytickformat('%0.1gM');ax.YRuler.Exponent = 0;box on;
    end
            %% Load data for 'All'
        vaccineRange=vaccineRangeVec{1};
        data=load(['../data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
        VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
        Rvec=data.Rvec;

        overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
        overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
        overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;

        
    defineColors;
    linkaxes([axVec{1} axVec{2}],'y');ylabel(axVec{1},'Overall infected','fontsize',11);
    linkaxes([axVec{3} axVec{4}],'y');ylabel(axVec{2},'Overall infected','fontsize',11);
    linkaxes([axVec{1} axVec{3}],'x');ylabel(axVec{3},'Mortality','fontsize',11);
    linkaxes([axVec{2} axVec{4}],'x');ylabel(axVec{4},'Mortality','fontsize',11);


    Rvec=data.Rvec;
    a=100;
    axVec{5}=nexttile(5);
         Mat=cell2mat(data.distribution_w0)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;  % title([vaccineRange,' - along w=0']);
          title('Allocations minimizing infections - all ages eligible');
         xlabel('R_0');ylabel('Number of vaccines');
         dataMat=(Mat.*data.Nadult)';
%         
         textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
        for ix=1:9
            hndl(ix).FaceColor=lineColors(ix,:);
        end
                      currY=0;

         for ix=1:7
             currY=currY+dataMat(a2,ix);
             if dataMat(a2,ix)>sum(dataMat(a2,:))/50
                 text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
             end
         end
                  hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)
         end
          ytickformat('%4gM')   
          xlim([R_threshold(1) 4]);ylim([0 182]);

    axVec{6}=nexttile(6);
         Mat=cell2mat(data.distribution_w1)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;  
         title('Allocations minimizing mortality - all ages eligible');
         title(axVec{2},'Allocations minimizing mortality');
         xlabel('R_0');ylabel('Number of vaccines');
         dataMat=(Mat.*data.Nadult)';    hold on;     plot([data.Rvec(a) data.Rvec(a)],[0 100],'w:')

 ytickformat('%4gM')        
         textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
         currY=0;
         for ix=2:9
             currY=currY+dataMat(a2,ix);
             if dataMat(a2,ix)>sum(dataMat(a2,:))/50
                 t(ix)=text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold','fontsize',10)
             end
         end
                for ix=1:9
            hndl(ix).FaceColor=lineColors(ix,:);
                end
        
         hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)
          xlim([R_threshold(1)-1e-4 4]);ylim([0 182.2]);
    lg = legend(axVec{1},'All','Above 10','Above 20','fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    lg.Layout.Tile='north'
    set(gcf, 'Position', [72 130 831 600]);
    sgtitle(['Vaccination of ',num2str(100*VcGnrl,2),'% of total population (',num2str(VcPrct),'% of adult population)']);
    
    switch susceptibilityFactor
        case 1
    set(axVec{1},'xtick',[1 1.5 R_threshold(3) 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','','R_{threshold}^{20+}','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{3},'xtick',[1 1.5 R_threshold(3) 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','','R_{threshold}^{20+}','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{2},'xtick',[1 1.5 R_threshold(3) 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','','R_{threshold}^{20+}','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{4},'xtick',[1 1.5 R_threshold(3) 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','','R_{threshold}^{20+}','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);

    
   set(axVec{5},'xtick',[R_threshold(1) 3 3.5 4],'xticklabel',{'            R_{threshold}','3','3.5','4'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 4]);
   set(axVec{6},'xtick',[R_threshold(1) 3 3.5 4],'xticklabel',{'            R_{threshold}','3','3.5','4'},'XTickLabelRotation',0);xlim(axVec{6},[R_threshold(1) 4]);
   set(t(2),'Position',[2.5653 15 0]);
set(t(4),'Position',[2.5653 26 0]);
set(t(9),'Position',[2.5653 185 0]);
        case 1.5
                set(axVec{1},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{3},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{2},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{4},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
set(axVec{5},'xtick',[R_threshold(1) 3 3.5 4],'xticklabel',{'            R_{threshold}','3','3.5','4'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 4]);
   set(axVec{6},'xtick',[R_threshold(1) 3 3.5 4],'xticklabel',{'            R_{threshold}','3','3.5','4'},'XTickLabelRotation',0);xlim(axVec{6},[R_threshold(1) 4]);
   set(t(9),'Position',[2.5462 185 0]);

        case 2
            
                set(axVec{1},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{3},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{2},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);
    set(axVec{4},'xtick',[1 R_threshold(3) 1.5 2 R_threshold(2) 3:0.5:4],'xticklabel',{'1','R_{threshold}^{20+}','','2','  R_{threshold}^{10+}','3','3.5','4'},'XTickLabelRotation',0);    
   set(axVec{5},'xtick',[R_threshold(1) 3 3.5 4],'xticklabel',{'            R_{threshold}','3','3.5','4'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 4]);
   set(axVec{6},'xtick',[R_threshold(1) 3 3.5 4],'xticklabel',{'            R_{threshold}','3','3.5','4'},'XTickLabelRotation',0);xlim(axVec{6},[R_threshold(1) 4]);

    end
shg


         title(axVec{1},'Allocations minizing infections');

         
    printGraph(['../graphs/OutcomesAsFuncOfR_',fname,',VcPrct_',num2str(10*VcPrct)]);
    shg
    

    return
    
