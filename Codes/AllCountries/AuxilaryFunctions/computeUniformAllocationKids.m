function uniformAllocation=computeUniformAllocationKids(s,N,vaccinesLeftToDistibute,upperBound,ageFilter)

ageFilterIdx=find(ageFilter>0);
sFiltered=s(ageFilterIdx).*ageFilter(ageFilterIdx);
Nfiltered=N(ageFilterIdx);
upperBoundFiltered=upperBound(ageFilterIdx);
uniformAllocation=[sFiltered.*Nfiltered/sum(sFiltered.*Nfiltered)*vaccinesLeftToDistibute]./Nfiltered;

idxVec=ones(1,numel(ageFilterIdx));
idx=find(uniformAllocation>upperBoundFiltered);

while numel(idx)>0
    vacExtra=(uniformAllocation(idx)-upperBoundFiltered(idx))'*Nfiltered(idx);
    uniformAllocation(idx)=max(upperBound(idx)-1e-15,0);
    idxVec(idx)=false;
    n_idx=find(idxVec);
    uniformAllocation(n_idx)=uniformAllocation(n_idx)+(sFiltered(n_idx).*Nfiltered(n_idx)/sum(sFiltered(n_idx).*Nfiltered(n_idx))*vacExtra)./Nfiltered(n_idx);
    idx=find(uniformAllocation>upperBoundFiltered);
end

uniformAllocation(ageFilterIdx)=uniformAllocation;
uniformAllocation(find(ageFilter==0))=0;
return

