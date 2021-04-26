function main
addpath('../Codes/AllCountries/AuxilaryFunctions/');
addpath('/Users/nirgavish/Dropbox (Technion Dropbox)/My Shared Functions')
close all
%% Parameters

betaVac=0.2; effVac=0.8;
recoveredprct=16;infected_nv_prct=0.0045;infected_v_prct=0.0045;fname='data_r16_95';

recoveredprct=8;infected_nv_prct=0.0045;infected_v_prct=0.0045;fname='data_r8_95';
% recoveredprct=0;infected_nv_prct=0;infected_v_prct=0;fname='data_no_r_no_i_95'
% The percent of 20+ population that have vaccines avaliable
maxPrct=95;%9.99; % Maximal vaccination per age group

%% Select Country
country="USA";
countryData=load(join(['../countryData/',country,'_data.mat'],''));

VcPrctVec=sort([35 50 65 70 75 85 30 40 60 80 90]);
vaccineRangeVec={'above20','above10','All'};%'above16','above12','above6'

%% Run loop
%     close all;

for kx=1:numel(VcPrctVec)
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
        %     subplot(5,2,1+2*(ix-1));
        overallInfected_nv_w1=data.overallInfected_nv_w1/1e6;
        overallInfected_nv_w0=data.overallInfected_nv_w0/1e6;
        overallInfected_nv_uniform=data.overallInfected_nv_uniform/1e6;
        
        R_threshold(kx,ix)=NaN;
        idx=min(find(overallInfected_nv_w0>0.5));
        if ~isempty(idx)
            R_threshold(kx,ix)=Rvec(min(find(overallInfected_nv_w0>0.5)));
        end
    end
end
save(['summary',fname]);
plot(R_threshold,'linewidth',1);legend('Above 20','Above 10','All');

    return
    
