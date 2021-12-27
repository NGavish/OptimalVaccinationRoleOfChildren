function outcomeAsFunctionofR(collectData,allOrNone,susceptibilityFactor)
% This function produce Figure 4 - "Impact of change in reproduction number." - and related figures in main text and SI
addpath('../Codes/AllCountries/AuxilaryFunctions/');
defineColors
close all

%% Parameters
sampleSize=2000;

if allOrNone
    nuVac=0.9;betaVac=0;
    fname='all_or_none'
    a2=400; % a2 determines the location of age groups label in the graphs
    Rvec=unique([8 linspace(1,10,sampleSize/2)]);
else
    nuVac=1;betaVac=0.1;
    fname='default'
    a2=900; % a2 determines the location of age groups label in the graphs
    Rvec=unique([linspace(1,10,sampleSize) 7.3+linspace(-0.05,0.05,101) 8.8+linspace(-0.1,0.1,401) 9.5+linspace(-0.05,0.05,51) 5.45+linspace(-0.05,0.05,51) 8.6+linspace(-0.05,0.05,51)]);
end

effVac=0.5;recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;
susFactorText=['susceptiblityFactor=',num2str(susceptibilityFactor*10)];

maxPrct=100; % Maximal vaccination per age group

%% Select Country
country="USA";
countryData=load(join(['../countryData/',country,'_data.mat'],''));

VcPrctVec=73.2;
vaccineRangeVec={'All','above10','above20'};

%% Run loop
if collectData
    [a,b]=mkdir('./data');
    count=0;
    M=numel(Rvec);
    for kx=1:numel(VcPrctVec)
        for jx=1:numel(vaccineRangeVec)
            tic
            VcPrct=VcPrctVec(kx);
            vaccineRange=vaccineRangeVec{jx};

            %% Prepare situation specific data for computation
            [uniformAllocation,riskFocusedAllocation,spreadersAllocation,upperBound,vaccinesLeftToDistibute,adultAges,ageAbove60,age40to60,age20to40,IFR,Cij,Ni,Nadult,r,v,infected0_v,infected0_nv]=prepareFinalSizeParameters(country,vaccineRange,recoveredprct,infected_nv_prct,infected_v_prct,VcPrct,betaVac,effVac,maxPrct,true);
            susceptibilityProfile=ones(size(Ni));
            if susceptibilityFactor==1.5
                susceptibilityProfile(2)=2;
            end
            if susceptibilityFactor==2
                susceptibilityProfile(1:2)=2;
                a2=100;
            end

            % Normalize contact matrix
            Cij=diag(susceptibilityProfile)*Cij;M2=diag(Ni)*Cij*diag(1./Ni);[V,d]=eig(M2);Cij =Cij/max(d(:));

            %% Define optimization problem
            problem.options = optimoptions('fmincon','MaxFunctionEvaluations',5e4,'ConstraintTolerance',1e-6,'StepTolerance',1e-10,'Display','none');%,'Algorithm','sqp','Display','none');%,'Display','iter');
            problem.solver = 'fmincon';
            problem.Aeq=Nadult';problem.Beq=vaccinesLeftToDistibute;
            problem.A=[];problem.B=[];
            problem.lb=0*upperBound;
            problem.ub=upperBound;

            %% Start scan from R=1 to R=10
            w=0;
            xw0=spreadersAllocation;
            xw1=riskFocusedAllocation;
            for ix=1:M
                R0=Rvec(ix)
                % First compute result for uniform allocation
                [dummy,overallInfected_nv_uniform(ix),overallFatality_uniform(ix),data_uniform{ix}]=computeFinalSize_generalized(uniformAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
                [dummy,overallInfected_nv_spreaders(ix),overallFatality_spreaders(ix),data_spreaders{ix}]=computeFinalSize_generalized(spreadersAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
                [dummy,overallInfected_nv_riskFocused(ix),overallFatality_riskFocused(ix),data_riskFocused{ix}]=computeFinalSize_generalized(riskFocusedAllocation,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

                % Find optimal infection minimizing allocation where
                % spreadersAllocation or riskFocusedAllocation serve as an initial guess
                w=0; y_infected_IC1=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),spreadersAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                %                 y_infected_IC1=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y_infected_IC1,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                [dummy,overallInfected_nv_w0_IC1,overallFatality_w0_IC1,data_w0_IC1]=computeFinalSize_generalized(y_infected_IC1,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

                w=1;y_mortality_IC1=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),riskFocusedAllocation,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                y_mortality_IC1=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y_mortality_IC1,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                [dummy,overallInfected_nv_w1_IC1,overallFatality_w1_IC1,data_w1_IC1]=computeFinalSize_generalized(y_mortality_IC1,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

                % Find optimal infection minimizing allocation where previous point which serves as an initial guess
                w=0; y_infected_IC2=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),xw0,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                [dummy,overallInfected_nv_w0_IC2,overallFatality_w0_IC2,data_w0_IC2]=computeFinalSize_generalized(y_infected_IC2,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

                w=1;y_mortality_IC2=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),xw1,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                [dummy,overallInfected_nv_w1_IC2,overallFatality_w1_IC2,data_w1_IC2]=computeFinalSize_generalized(y_mortality_IC2,adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);

                % Take the minimum between the two - infections
                if overallInfected_nv_w0_IC1<overallInfected_nv_w0_IC2
                    overallInfected_nv_w0(ix)=overallInfected_nv_w0_IC1;
                    overallFatality_w0(ix)=overallFatality_w0_IC1;
                    data_w0{ix}=data_w0_IC1;
                    xw0=y_infected_IC1;
                else
                    overallInfected_nv_w0(ix)=overallInfected_nv_w0_IC2;
                    overallFatality_w0(ix)=overallFatality_w0_IC2;
                    data_w0{ix}=data_w0_IC2;
                    xw0=y_infected_IC2;
                end
                distribution_w0{ix}=xw0;
                
                % In case of a sudden jump - verify that minimum solution
                % is indeed reached
                if ix>1 & (abs(overallFatality_w0(ix-1)-overallFatality_w0(ix))>2e5 | max(abs(distribution_w0{ix-1}-distribution_w0{ix}))>0.1)
                    parfor lx=1:100
                        y=randomizePoint(distribution_w0{ix},Ni,Ni(adultAges),upperBound,vaccinesLeftToDistibute,0.05+mod(lx,3)/10);
                        [x_aux{lx},overallInfections_sample(lx),exitflag,output]=fmincon(@(x) computeFinalSize_generalized(x,adultAges,IFR,0,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1),y,problem.A,problem.B,problem.Aeq,problem.Beq,problem.lb,problem.ub,[],problem.options);
                        if exitflag<=0
                            overallInfections_sample(lx)=Inf;
                        end
                    end

                    [C,idx]=min(overallInfections_sample);
                    if overallInfections_sample(idx)<overallInfected_nv_w0(ix)
                        1-overallInfections_sample(idx)/overallInfected_nv_w0(ix)
                        optimalAllocation{ix}=x_aux{idx};
                        [dummy,overallInfected_nv_w0(ix),overallFatality_w0(ix),data_w0{ix}]=computeFinalSize_generalized(optimalAllocation{ix},adultAges,IFR,w,R0,Cij,Ni,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,1,1);
                        idx
                    end
                end

                % Take the minimum between the two - mortality
                if overallFatality_w1_IC1<overallFatality_w1_IC2
                    overallInfected_nv_w1(ix)=overallInfected_nv_w1_IC1;
                    overallFatality_w1(ix)=overallFatality_w1_IC1;
                    data_w1{ix}=data_w1_IC1;
                    xw1=y_mortality_IC1;
                else
                    overallInfected_nv_w1(ix)=overallInfected_nv_w1_IC2;
                    overallFatality_w1(ix)=overallFatality_w1_IC2;
                    data_w1{ix}=data_w1_IC2;
                    xw1=y_mortality_IC2;
                end
                distribution_w1{ix}=xw1;
            end

            save(['./data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
            toc
        end
    end
end

% Produce graphs
lineStyle={'-','--','-.'}
a=100
for kx=1:numel(VcPrctVec)
    close all;
    t=tiledlayout(5,2, 'TileSpacing','compact','Padding','none')
    VcPrct=VcPrctVec(kx);

    for ix=1:numel(vaccineRangeVec)

        vaccineRange=vaccineRangeVec{ix};
        data=load(['./data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
        VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
        Rvec=data.Rvec;

        %% Present Data
        green=[0.4660, 0.6740, 0.1880];orange=[0.8500, 0.3250, 0.0980];blue=[0, 0.4470, 0.7410];gray=0.5*[1 1 1];darkgray=0.4*[1 1 1];
        overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
        overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
        overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;
        overallInfected_nv_spreaders=data.overallInfected_nv_spreaders/1e6;

        R_threshold(ix)=Rvec(min(find(overallInfected_nv_w0>0.5)));

        % Infection minimizing
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
    data=load(['./data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
    VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
    Rvec=data.Rvec;

    overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
    overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
    overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;
    overallInfected_nv_spreaders=data.overallInfected_nv_spreaders/1e6;


    defineColors;
    linkaxes([axVec{1} axVec{2}],'y');ylabel(axVec{1},'Overall infected','fontsize',11);
    linkaxes([axVec{3} axVec{4}],'y');ylabel(axVec{2},'Overall infected','fontsize',11);
    linkaxes([axVec{1} axVec{3}],'x');ylabel(axVec{3},'Mortality','fontsize',11);
    linkaxes([axVec{2} axVec{4}],'x');ylabel(axVec{4},'Mortality','fontsize',11);


    Rvec=data.Rvec;
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
    currY=4;

    for ix=2:7
        currY=currY+dataMat(a2,ix);
        if dataMat(a2,ix)>sum(dataMat(a2,:))/50
            text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
        end
    end
    hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)
end
ytickformat('%4gM')

axVec{6}=nexttile(6);
Mat=cell2mat(data.distribution_w1)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;
title('Allocations minimizing mortality - all ages eligible');
title(axVec{2},'Allocations minimizing mortality');
xlabel('R_0');ylabel('Number of vaccines');
dataMat=(Mat.*data.Nadult)';    hold on;     plot([data.Rvec(a) data.Rvec(a)],[0 100],'w:')

ytickformat('%4gM')
textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
currY=5;
for ix=1:8
    currY=currY+dataMat(a2,ix);
    if dataMat(a2,ix)>sum(dataMat(a2,:))/40
        t(ix)=text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold','fontsize',10)
    end
end
for ix=1:9
    hndl(ix).FaceColor=lineColors(ix,:);
end

hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)
xlim([R_threshold(1)-1e-4 10]);ylim([0 182.2]);


%% Load data for 'above20'
vaccineRange=vaccineRangeVec{3};
data=load(['./data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
Rvec=data.Rvec;

overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;


Rvec=data.Rvec;
a=100;
axVec{9}=nexttile(9);
Mat=cell2mat(data.distribution_w0)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;  % title([vaccineRange,' - along w=0']);
title('Allocations minimizing infections - ages 20 and older eligible');
xlabel('R_0');ylabel('Number of vaccines');
dataMat=(Mat.*data.Nadult)';
%
textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
for ix=1:7
    hndl(ix).FaceColor=lineColors(2+ix,:);
end
currY=-8;

for ix=1:4
    currY=currY+dataMat(a2,ix);
    if dataMat(a2,ix)>sum(dataMat(a2,:))/30
        text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
    end
end
hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)

ytickformat('%4gM')
xlim([R_threshold(3) 10]);ylim([0 182]);

axVec{10}=nexttile(10);
Mat=cell2mat(data.distribution_w1)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;
title('Allocations minimizing mortality - ages 20 and older eligible');
title(axVec{2},'Allocations minimizing mortality');
xlabel('R_0');ylabel('Number of vaccines');
dataMat=(Mat.*data.Nadult)';    hold on;     plot([data.Rvec(a) data.Rvec(a)],[0 100],'w:')

ytickformat('%4gM')
textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
currY=5;
for ix=1:6
    currY=currY+dataMat(a2,ix);
    if dataMat(a2,ix)>sum(dataMat(a2,:))/50
        t2(ix)=text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold','fontsize',10)
    end
end
for ix=1:7
    hndl(ix).FaceColor=lineColors(ix+2,:);
end

hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)
xlim([R_threshold(3)-1e-4 10]);ylim([0 182.2]);

%% Load data for 'above10'
vaccineRange=vaccineRangeVec{2};
data=load(['./data/OutcomesAsFuncOfR_',fname,'VcPrct_',num2str(10*VcPrct),'_',vaccineRange,'_',susFactorText]);
VcGnrl=data.vaccinesLeftToDistibute/sum(data.Ni);
Rvec=data.Rvec;

overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;


Rvec=data.Rvec;

axVec{7}=nexttile(7);
Mat=cell2mat(data.distribution_w0)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;  % title([vaccineRange,' - along w=0']);
title('Allocations minimizing infections - ages 10 and older eligible');
xlabel('R_0');ylabel('Number of vaccines');
dataMat=(Mat.*data.Nadult)';
%
textAges={'0-9','10-19','20-29','30-39','40-49','50-59','60-69','70-79','80+'};
for ix=2:8
    hndl(ix).FaceColor=lineColors(1+ix,:);
end
currY=0;

for ix=1:4
    currY=currY+dataMat(a2,ix);
    if dataMat(a2,ix)>sum(dataMat(a2,:))/30
        text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold')
    end
end
hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)

ytickformat('%4gM')
xlim([R_threshold(1) 10]);ylim([0 182]);

axVec{8}=nexttile(8);
Mat=cell2mat(data.distribution_w1)/1e6;hndl=area(data.Rvec,(Mat.*data.Nadult)') ;
title('Allocations minimizing mortality - ages 10 and older eligible');
title(axVec{2},'Allocations minimizing mortality');
xlabel('R_0');ylabel('Number of vaccines');
dataMat=(Mat.*data.Nadult)';    hold on;     plot([data.Rvec(a) data.Rvec(a)],[0 100],'w:')

ytickformat('%4gM')

currY=5;
for ix=1:7
    currY=currY+dataMat(a2,ix);
    if dataMat(a2,ix)>sum(dataMat(a2,:))/50
        t2(ix)=text(data.Rvec(a2),currY,['  ',textAges{data.adultAges(ix)},' (',num2str(100*1e6*Mat(ix,a2),3),'%)'],'verticalAlign','top','Color','w','fontweight','bold','fontsize',10)
    end
end
for ix=2:8
    hndl(ix).FaceColor=lineColors(1+ix,:);
end

hold on;plot([data.Rvec(a2) data.Rvec(a2)],[0 182],'w:','linewidth',2)
xlim([R_threshold(1)-1e-4 10]);ylim([0 182.2]);

%%%%%%%%%%%%%%%%
lg = legend(axVec{1},'All','Above 10','Above 20','fontsize',10,'Orientation','horizontal');%legend(nexttile(2), [line1,line2,line3,line4]);
lg.Layout.Tile='north'

lg = legend(axVec{6},textAges,'Orientation','horizontal')
lg.Layout.Tile='south'
set(gcf, 'Position', [66 72 831 725]);
sgtitle(['Vaccination of ',num2str(100*VcGnrl,2),'% of total population (',num2str(VcPrct),'% of adult population)']);

switch susceptibilityFactor
    case 1
        for ix=1:4
                        xlim(axVec{ix},[R_threshold(2)  10])

            set(axVec{ix},'xtick',[1 R_threshold(3) R_threshold(2) 4:10],'xticklabel',{'1','R_{threshold}^{20+}','R_{threshold}^{10+}','4','5','6','7','8','9','10'},'XTickLabelRotation',0);xlim(axVec{ix},[1,10])
        end
        for ix=5:8
            xlim(axVec{ix},[R_threshold(2)  10])
            set(axVec{ix},'xtick',[R_threshold(2) 4:10],'xticklabel',{'         R_{threshold}^{10+}','4','5','6','7','8','9','10'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 10]);
        end

        for ix=9:10
            set(axVec{ix},'xtick',[R_threshold(3) 3:10],'xticklabel',{'         R_{threshold}^{20+}','3','4','5','6','7','8','9','10'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 10]);
        end
        % set(t(2),'Position',[2.5653 15 0]);
        % set(t(4),'Position',[2.5653 26 0]);
        % set(t(9),'Position',[2.5653 185 0]);delete(t(2));
    case 1.5

            xlim(axVec{5},[R_threshold(2)  10])

sgtitle(['Adolescents (age group 10-19) are equally susceptible as adults']);
for ix=1:10
         set(axVec{ix},'xtick',[1 R_threshold(3) 2 R_threshold(2)+1e-2 4:10],'xticklabel',{'','R_{threshold}^{20+}','2','R_{threshold}^{10+}','4','5','6','7','8','9','10'},'XTickLabelRotation',0);
end

delete(axVec{5}.Children(2));
    case 2
       for ix=1:4
           xlim(axVec{ix},[1  10])
            set(axVec{ix},'xtick',[1 R_threshold(3) 2 R_threshold(2) R_threshold(1) 4:10],'xticklabel',{'','R_{threshold}^{20+}','','','R_{threshold}','4','5','6','7','8','9','10'},'XTickLabelRotation',0);
        end
        for ix=5:6
            xlim(axVec{ix},[R_threshold(1)  10])
            set(axVec{ix},'xtick',[R_threshold(1) 4:10],'xticklabel',{'R_{threshold}','4','5','6','7','8','9','10'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 10]);ylim(axVec{ix},[0 data.vaccinesLeftToDistibute/1e6]);
        end
        for ix=7:8
            xlim(axVec{ix},[R_threshold(2)  10])
            set(axVec{ix},'xtick',[R_threshold(2) R_threshold(1) 4:10],'xticklabel',{'R_{threshold}^{10+}','','4','5','6','7','8','9','10'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 10]);ylim(axVec{ix},[0 data.vaccinesLeftToDistibute/1e6]);
        end


        for ix=9:10
            set(axVec{ix},'xtick',[R_threshold(3) 3:10],'xticklabel',{'         R_{threshold}^{20+}','3','4','5','6','7','8','9','10'},'XTickLabelRotation',0);xlim(axVec{5},[R_threshold(1) 10]);;ylim(axVec{ix},[0 data.vaccinesLeftToDistibute/1e6]);
        end

        axVec{5}.Children(2).Position=axVec{5}.Children(2).Position+[0 15 0]
        axVec{5}.Children(3).Position=axVec{5}.Children(3).Position+[0 10 0]
        axVec{5}.Children(4).Position=axVec{5}.Children(4).Position+[0 20 0]
        axVec{5}.Children(5).Position=axVec{5}.Children(5).Position+[0 15 0]


        sgtitle(['Children (age group 0-19) are equally susceptible as adults']);
end
shg


title(axVec{1},'Allocations minimizing infections');
ylim(axVec{5},[0 182.2])

printGraph(['../graphs/OutcomesAsFuncOfR_',fname,',VcPrct_',num2str(10*VcPrct),'_',susFactorText]);
shg


return

