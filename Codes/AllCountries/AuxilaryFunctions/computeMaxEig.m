function result=computeMaxEig(vacDistParm,adultAges,IFR,w,R0,C,N,r,v,infected0_v,infected0_nv,betaVac,nuVac,effVac,overallInfected_nv_uniform,overallFatality_uniform)
vacDist=0*r;vacDist(adultAges)=vacDistParm;

infected=infected0_nv+infected0_v;
s=1-r-infected-nuVac*(1-betaVac)*(v+vacDist);
A=diag(s)*C;[V,d]=eig(A);

result=max(d(:));
return
