function [result,overallInfected,overallFatality,data]=computeFinalSize(vacDistParm,adultAges,IFR,w,R0,C,N,r,v,infected0_v,infected0_nv,betaVac,effVac,overallInfected_uniform,overallFatality_uniform)
vacDist=0*r;vacDist(adultAges)=vacDistParm;

infected=infected0_nv+infected0_v;
s=1-r-infected-v-vacDist;
A=R0*C;
r_inf=real(fsolve(@(z) g(z,A,r,v+vacDist,infected,betaVac),r+0.8*s,optimset('Display','none','TolX',1e-8)));
s_inf=s.*exp(-A*(r_inf-r));
v_inf=(v+vacDist).*exp(-betaVac*A*(r_inf-r));
infected_nv=s-s_inf+infected0_nv;infected_v=vacDist+v-v_inf+infected0_v;
sus_v=v_inf;
sus_nv=s_inf;
overallInfected_nv=(infected_nv+0*infected_v)'*N;
overallInfected=(infected_nv+infected_v)'*N;
overallFatality=((infected_nv+(1-effVac)*infected_v).*N)'*IFR; % Overall number of deaths

result=(1-w)*overallInfected/overallInfected_uniform+w*overallFatality/overallFatality_uniform;

data=[sus_v.*N infected_v.*N r.*N infected_nv.*N sus_nv.*N];data=[data N-sum(data,2)];

return

function y = g(rinf,A,r,v,infected,betaVac)
s=1-r-v-infected;
y=rinf+s.*exp(-A*(rinf-r))+v.*exp(-betaVac*A*(rinf-r))-1;
% s_inf_normalized=exp(-A*(z-r./s));
% v_inf_normalized=v./s.*s_inf_normalized.^betaVac;
% y=1./s-z-s_inf_normalized-v_inf_normalized;
return

