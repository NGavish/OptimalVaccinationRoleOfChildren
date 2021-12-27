function uniformAllocation=computeUniformAllocation(s,adultAges,Nadult,vaccinesLeftToDistibute,upperBound)

uniformAllocation=[s(adultAges).*Nadult/sum(s(adultAges).*Nadult)*vaccinesLeftToDistibute]./Nadult;

idxVec=ones(1,numel(adultAges));
idx=find(uniformAllocation>upperBound);

while numel(idx)>0
    vacExtra=(uniformAllocation(idx)-upperBound(idx))'*Nadult(idx);
    uniformAllocation(idx)=upperBound(idx)-1e-15;
    idxVec(idx)=false;
    n_idx=find(idxVec);
    uniformAllocation(n_idx)=uniformAllocation(n_idx)+(s(n_idx).*Nadult(n_idx)/sum(s(n_idx).*Nadult(n_idx))*vacExtra)./Nadult(n_idx);
    idx=find(uniformAllocation>upperBound);
end

uniformAllocation=max(uniformAllocation,0);

return

