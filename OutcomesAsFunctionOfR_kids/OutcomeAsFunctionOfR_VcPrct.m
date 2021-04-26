function main
addpath('../Codes/AllCountries/AuxilaryFunctions/');
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
close all
%% Parameters
collectData=;
fname='data'

betaVac=0.2; effVac=0.8;
recoveredprct=16;infected_nv_prct=0.0045;infected_v_prct=0.0045;fname='data_r16_95';
recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;fname='data_no_r_no_i_95'
% The percent of 20+ population that have vaccines avaliable
maxPrct=99.9;%9.99; % Maximal vaccination per age group

%% Select Country
country="USA";
countryData=load(join(['../countryData/',country,'_data.mat'],''));

VcPrctVec=sort([35 50 65 70 75 85 30 40 60 80 90]);
vaccineRangeVec={'above20','above10','All'};%'above16','above12','above6'

%% Run loop
if collectData
    sampleSize=2500;
    Rvec=linspace(1,4,200);M=numel(Rvec);
    for kx=1:numel(VcPrctVec)
        for jx=1:numel(vaccineRangeVec)
            tic
            VcPrct=VcPrctVec(kx);
            vaccineRange=vaccineRangeVec{jx};
            
            %% Prepare data
            [uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,VcPrct,betaVac,effVac,maxPrct,true);
            
            xw0=uniformAllocation;%spreadersAllocation;
            xw1=uniformAllocation;%riskFocusedAllocation;
            
            %% Define optimization problem
            problem.options = optimoptions('fmincon','MaxFunctionEvaluations',1e4,'Display','none');%,'Display','iter');
            problem.solver = 'fmincon';
            problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
            problem.A=[];problem.B=[];
            problem.lb=0*upperBound;
            problem.ub=upperBound;
            
            % Start low-res scan from R=1 to R=4
            w=0
            xw0=uniformAllocation;
            xw1=uniformAllocation;
            for ix=1:M
                R0=Rvec(ix);
                [result,overallInfected_nv_uniform(ix),overallFatality_uniform(ix),data_uniform{ix}]=computeFinalSize(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                
                parfor jx=1:sampleSize
                    y=randomizePoint(xw0,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.2+0.8*mod(jx,2)));
                    [sample_result,sample_overallInfected_nv(jx),sample_overallFatality(jx),sample_data{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                    sample_distribution{jx}=y;
                    
                    y=randomizePoint(xw1,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.2+0.8*mod(jx,2)));
                    [sample_result(jx),sample_overallInfected_nv2(jx),sample_overallFatality2(jx),sample_data2{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
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
                    y=randomizePoint(xw0,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.2+0.8*mod(jx,2)));
                    [sample_result(jx),sample_overallInfected_nv(jx),sample_overallFatality(jx),sample_data{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
                    sample_distribution{jx}=y;
                    
                    y=randomizePoint(xw1,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.2+0.8*mod(jx,2)));
                    [sample_result2(jx),sample_overallInfected_nv2(jx),sample_overallFatality2(jx),sample_data2{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
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
            
            
            %          for ix=1:M
            %             R0=Rvec(ix);
            %
            %             %% Compute Uniform
            %             [result,overallInfected_nv_uniform(ix),overallFatality_uniform(ix),data_uniform{ix}]=computeFinalSize(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
            %
            %             %% Compute xw0
            %             w=0;
            %             parfor jx=1:sampleSize
            %                 y=randomizePoint(xw0,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.2+0.8*mod(jx,2)));
            %                 [sample_result(jx),sample_overallInfected_nv(jx),sample_overallFatality(jx),sample_data{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
            %                 sample_distribution{jx}=y;
            %             end
            %             [minResult,minIdx]=min(sample_result);
            %             y=fmincon(@(x) computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1),sample_distribution{minIdx},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
            %             [result,overallInfected_nv_w0(ix),overallFatality_w0(ix),data_w0{ix}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
            %             distribution_w0{ix}=y;
            %
            %             if minResult<result
            %                 display('Optimization wasn''t optimal');
            %                 overallInfected_nv_w0(ix)=sample_overallInfected_nv(minIdx);
            %                 overallFatality_w0(ix)=sample_overallFatality(minIdx);
            %                 data_w0{ix}=sample_data{minIdx};
            %                 distribution_w0{ix}=sample_distribution{minIdx};
            %             end
            %
            %             xw0=distribution_w0{ix};
            %             %% Compute xw1
            %             w=1;
            %             parfor jx=1:sampleSize
            %                 y=randomizePoint(xw1,Ni,Nadult,upperBound,vaccinesLeftToDistibute,(jx>1).*(0.2+0.8*mod(jx,2)));
            %                 [sample_result(jx),sample_overallInfected_nv(jx),sample_overallFatality(jx),sample_data{jx}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
            %                 sample_distribution{jx}=y;
            %             end
            %             [minResult,minIdx]=min(sample_result);
            %             y=fmincon(@(x) computeFinalSize(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1),sample_distribution{minIdx},problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
            %             [result,overallInfected_nv_w1(ix),overallFatality_w1(ix),data_w1{ix}]=computeFinalSize(y,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,betaVac,effVac,1,1);
            %             distribution_w1{ix}=y;
            %
            %             if minResult<result
            %                 display('Optimization wasn''t optimal');
            %                 overallInfected_nv_w1(ix)=sample_overallInfected_nv(minIdx);
            %                 overallFatality_w1(ix)=sample_overallFatality(minIdx);
            %                 data_w1{ix}=sample_data{minIdx};
            %                 distribution_w1{ix}=sample_distribution{minIdx};
            %             end
            %
            %             xw1=distribution_w1{ix};
            %     end
            
            
            save([fname,'VcPrct_',num2str(VcPrct),'_',vaccineRange]);
            toc
        end
    end
end

lineStyle={'-','-','-'}

% vaccineRange=vaccineRangeVec{1};VcPrctVec(1)
% data=load([fname,'VcPrct_',num2str(VcPrctVec(1)),'_',vaccineRange]);
% overallInfected_nv_w1_base=data.overallInfected_nv_w1/1e6;

close all;
t=tiledlayout(3,6, 'TileSpacing','compact','Padding','none')
for kx=9 %;1:numel(VcPrctVec)
    VcPrct=VcPrctVec(kx);
%     overallFatality_w1_base=data.overallFatality_w1/1e6;
    for ix=1:numel(vaccineRangeVec)
        
        data=load([fname,'VcPrct_',num2str(VcPrct),'_',vaccineRangeVec{1}]);
        overallInfected_nv_w1_base=data.overallInfected_nv_w1/1e6;
        overallFatality_w1_base=data.overallFatality_w1/1e6;
        vaccineRange=vaccineRangeVec{ix};
        data=load([fname,'VcPrct_',num2str(VcPrct),'_',vaccineRange]);
        VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
        Rvec=data.Rvec;
        %% Present Data
        green=[0.4660, 0.6740, 0.1880];orange=[0.8500, 0.3250, 0.0980];blue=[0, 0.4470, 0.7410];gray=0.5*[1 1 1];darkgray=0.4*[1 1 1];
        axVec{ix}=nexttile(1+2*(ix-1),[1 2]);
        %     subplot(5,2,1+2*(ix-1));
        overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
        overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
        overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;
        
        R_threshold=Rvec(min(find(overallInfected_nv_w0>0.5)));
        
        h=patch([Rvec Rvec(end:-1:1)],[overallInfected_nv_w0 overallInfected_nv_w1(end:-1:1)],1,'FaceColor',gray,'EdgeColor','none');hold on;
        plot(Rvec,overallInfected_nv_w1,'color',orange','linewidth',2,'linestyle','-');hold on;
        plot(Rvec,overallInfected_nv_w0,'color',green,'linewidth',2,'linestyle','--');hold on;
        plot(Rvec,overallInfected_nv_uniform,'color',blue,'linewidth',2,'linestyle','-.');hold on;
        plot(Rvec,overallInfected_nv_w1_base,'color',darkgray','linewidth',2,'linestyle',':');
        alpha(0.2); grid on;box on;ax = gca;ax.YRuler.Exponent = 0;
        xlabel('R');%ylabel('Non-vaccinated infectives');
        ytickformat('%gM');%ylim([0 85])
        %title('Additional non-vaccinated infectives');
        
        % legend('Allocations along Pareto front','Uniform allocations','location','best');
        box on;
        %     subplot(5,2,2+2*(ix-1));
        %subplot(3,4,ix+4);
        axVec{ix+3}=nexttile(7+2*(ix-1),[1 2]);
        %axVec{ix+3}=nexttile(ix+3);
        h=patch([Rvec Rvec(end:-1:1)],[data.overallFatality_w0 data.overallFatality_w1(end:-1:1)]/1e6,1,'FaceColor',gray,'EdgeColor','none');hold on;
        
        plot(Rvec,data.overallFatality_w1/1e6,'color',orange','linewidth',2,'linestyle','-');hold on;
        plot(Rvec,data.overallFatality_w0/1e6,'color',green,'linewidth',2,'linestyle','--');hold on;
        plot(Rvec,data.overallFatality_uniform/1e6,'color',blue,'linewidth',2,'linestyle','-.');hold on;
        
        
        plot(Rvec,overallFatality_w1_base,'color',darkgray,'linewidth',2,'linestyle',':');
        alpha(0.2); grid on;box on;ax = gca;ax.YRuler.Exponent = 0;
        xlabel('R');%ylabel('Mortality');
        ytickformat('%0.1gM');ax.YRuler.Exponent = 0;box on;
        %ytickformat('%,4.0f');ax.YRuler.Exponent = 0;
        %ylim([0 1.09]);
        % legend('Allocations along Pareto front','w=0','Uniform allocation','w=1','location','best','fontsize',10);
%         figure(2);
%         Rvec=data.Rvec;
%         a=round((R_threshold+Rvec(end))*numel(Rvec)/Rvec(end)/2)
%         subplot(1,2,1);
%         
%         Mat=cell2mat(data.distribution_w0);area(data.Rvec,(Mat.*data.Nadult)') ;   title([vaccineRange,' - along w=0']);
%         xlabel('R');ylabel('Number of vaccines');
%         dataMat=(Mat.*data.Nadult)';
%         
%         textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
%         currY=0;
%         for ix=1:7
%             currY=currY+dataMat(a,ix);
%             if dataMat(a,ix)>sum(dataMat(a,:))/50
%                 text(data.Rvec(a),currY,[textAges{data.adultAges(ix)},' (',num2str(100*Mat(ix,a),3),'%)'],'verticalAlign','top')
%             end
%         end
%         %      legend('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+');
%          xlim([R_threshold 4]);
%         subplot(1,2,2);
%         
%         Mat=cell2mat(data.distribution_w1);area(data.Rvec,(Mat.*data.Nadult)') ;   title([vaccineRange,' - along w=1']);
%         xlabel('R');ylabel('Number of vaccines');
%         dataMat=(Mat.*data.Nadult)';
%         
%         textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
%         
%         currY=0;
%         for ix=1:numel(data.adultAges)
%             currY=currY+dataMat(a,ix);
%             if dataMat(a,ix)>sum(dataMat(a,:))/50
%                 text(data.Rvec(a),currY,[textAges{data.adultAges(ix)},' (',num2str(100*Mat(ix,a),3),'%)'],'verticalAlign','top')
%             end
%         end
%  xlim([R_threshold 4]);
        
        
        %     subplot(3,4,1);title('Ages 20 and older');ylabel('Infected');
        %     subplot(3,4,2);title('Ages 10 and older');
        %     subplot(3,4,3);title('All ages');
        %     subplot(3,4,5);ylabel('Mortality');
        %     legend('Allocations along Pareto front','w=1 (minimum moratility)','w=0 (minimum infections)','Uniform allocation','Baseline curve','fontsize',11);
    %end 
    %sgtitle(['Vaccination of ',num2str(VcPrct),'% of adult population (',num2str(100*VcGnrl,3),'% of total population)']);
    %printGraph([fname,'allocationsAsFunctionOfR_,VcPrct_',num2str(VcPrct)]);
    end
    linkaxes([axVec{1} axVec{2} axVec{3}],'y');ylabel(axVec{1},'Non-vaccinated Infected','fontsize',11);
    linkaxes([axVec{4} axVec{5} axVec{6}],'y');ylabel(axVec{4},'Mortality','fontsize',11);
    title(axVec{1},'Ages 20 and older');
    title(axVec{2},'Ages 10 and older');
    title(axVec{3},'All Ages');
    
    for ix=1:3
        text(axVec{ix},0.1,0.9,char(2*(ix-1)+'A'),'units','normalized','fontsize',13)
        text(axVec{ix+3},0.1,0.9,char(2*(ix-1)+'B'),'units','normalized','fontsize',13)
    end
    
    Rvec=data.Rvec;
    a=180;
    axVec{7}=nexttile(13,[1,3]);
%         a=round((R_threshold+Rvec(end))*numel(Rvec)/Rvec(end)/2)
%         subplot(1,2,1);
%         
         Mat=cell2mat(data.distribution_w0)/1e6;area(data.Rvec,(Mat.*data.Nadult)') ;   title([vaccineRange,' - along w=0']);
         xlabel('R_0');ylabel('Number of vaccines');
         dataMat=(Mat.*data.Nadult)';
%         
         textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
         currY=0;
         for ix=1:7
             currY=currY+dataMat(a,ix);
             if dataMat(a,ix)>sum(dataMat(a,:))/50
                 text(data.Rvec(a),currY,[textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
             end
         end
          ytickformat('%4gM')   
%         %      legend('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+');
          xlim([R_threshold 4]);

    axVec{8}=nexttile(16,[1,3]);
%         a=round((R_threshold+Rvec(end))*numel(Rvec)/Rvec(end)/2)
%         subplot(1,2,1);
%         
         Mat=cell2mat(data.distribution_w1)/1e6;area(data.Rvec,(Mat.*data.Nadult)') ;   title([vaccineRange,' - along w=0']);
         xlabel('R_0');ylabel('Number of vaccines');
         dataMat=(Mat.*data.Nadult)';
 ytickformat('%4gM')        
         textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
         currY=0;
         for ix=1:9
             currY=currY+dataMat(a,ix);
             if dataMat(a,ix)>sum(dataMat(a,:))/50
                 text(data.Rvec(a),currY,[textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
             end
         end
%         %      legend('0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+');
          xlim([R_threshold 4]);

    lg = legend(nexttile(3),'Allocations along Pareto front','w=1 (minimum moratility)','w=0 (minimum infections)','Uniform allocation','Baseline curve','fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
    lg.Location = 'northoutside';
    set(gcf, 'Position', [72 130 766 563]);
    sgtitle(['Vaccination of ',num2str(100*VcGnrl,2),'% of total population (',num2str(VcPrct),'% of adult population)']);
    printGraph([fname,'outComesAsFunctionOfR,VcPrct_',num2str(VcPrct)]);
    % fixepsbbox(outcomeAsFunctionOfR);
    shg
    end

    return
    
