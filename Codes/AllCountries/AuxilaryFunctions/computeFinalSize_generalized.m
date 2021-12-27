function [result,overallInfected,overallFatality,data]=computeFinalSize_generalized(vacDistParm,adultAges,IFR,w,R0,C,N,r,v,infected0_v,infected0_nv,nuVac,betaVac,effVac,overallInfected_uniform,overallFatality_uniform)
vacDist=0*r;vacDist(adultAges)=vacDistParm;

infected=infected0_nv+infected0_v;

v_protected=nuVac*(v+vacDist);
v_vunerable=(1-nuVac)*(v+vacDist);
A=R0*C;
s=max(1e-16,1-r-v_protected-infected);
s_nv=max(1e-16,1-r-v_protected-v_vunerable-infected);

%r_inf=real(fsolve(@(z) g(z,A,r,v_protected,infected,betaVac),(1+r)/2,optimset('Display','none','TolX',1e-8)));
r_inf=real(fsolve(@(z) g(z,A,r,v_protected,infected,betaVac),1+0*r,optimset('Display','none','TolX',1e-8)));
s_nv_inf=s_nv.*exp(-A*(r_inf-r));
v_inf=v_protected.*exp(-betaVac*A*(r_inf-r))+v_vunerable.*exp(-A*(r_inf-r));
infected_nv=s_nv-s_nv_inf+infected0_nv;
infected_v=vacDist+v-v_inf+infected0_v;

sus_v=v_inf;
sus_nv=s_nv_inf;
overallInfected_nv=(infected_nv+0*infected_v)'*N;
overallInfected=(infected_nv+infected_v)'*N;
overallFatality=((infected_nv+nuVac*(1-effVac)*infected_v+(1-nuVac)*infected_v).*N)'*IFR; % Overall number of deaths

result=(1-w)*overallInfected/overallInfected_uniform+w*overallFatality/overallFatality_uniform;

data=[sus_v.*N infected_v.*N r.*N infected_nv.*N sus_nv.*N];data=[data N-sum(data,2)];
%NG
return

function y = g(rinf,A,r,v,infected,betaVac)
s=1-r-v-infected;
y=rinf+s.*exp(-A*(rinf-r))+v.*exp(-betaVac*A*(rinf-r))-1;
% s_inf_normalized=exp(-A*(z-r./s));
% v_inf_normalized=v./s.*s_inf_normalized.^betaVac;
% y=1./s-z-s_inf_normalized-v_inf_normalized;
return

