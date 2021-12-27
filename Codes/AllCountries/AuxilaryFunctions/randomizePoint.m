function x=randomizePoint(xSeed,N,Nadult,upperBound,vaccinesLeftToDistibute,delta)

% while sum(x)==0 
x=min(max(xSeed+delta*(rand(size(xSeed))-0.5),0),upperBound);
% end
x=vaccinesLeftToDistibute/sum(x.*Nadult)*x;

idxVec=upperBound>0;
idx=find(x>upperBound & upperBound>0);

while numel(idx)>0
    vacExtra=(x(idx)-upperBound(idx))'*Nadult(idx);
    x(idx)=max(upperBound(idx)-1e-15,0);
    idxVec(idx)=false;
    n_idx=find(idxVec);
    x(n_idx)=x(n_idx)+(upperBound(n_idx).*Nadult(n_idx)/sum(upperBound(n_idx).*Nadult(n_idx))*vacExtra)./Nadult(n_idx);
    idx=find(x>upperBound & upperBound>0);
end

if isnan(x(1))
    x=xSeed;
end
return