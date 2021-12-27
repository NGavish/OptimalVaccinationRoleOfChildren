function riskFocusedAllocation=computeRiskFocusedAllocation(s,adultAges,Nadult,vaccinesLeftToDistibute,upperBound)

riskFocusedAllocation=Nadult*0;
ix=numel(adultAges);
while vaccinesLeftToDistibute>0 & ix>0
    jx=adultAges(ix);
    currVacDist=min([s(jx).*Nadult(ix),upperBound(ix).*Nadult(ix),vaccinesLeftToDistibute]);
    riskFocusedAllocation(ix)=currVacDist./Nadult(ix);
    vaccinesLeftToDistibute=vaccinesLeftToDistibute-currVacDist;
    ix=ix-1;
end
% warning([num2str(vaccinesLeftToDistibute),' were not allocated']);
return

