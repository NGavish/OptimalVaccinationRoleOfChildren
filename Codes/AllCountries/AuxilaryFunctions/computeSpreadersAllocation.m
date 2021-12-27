function spreadersAllocation=computeSpreadersAllocation(s,adultAges,Nadult,vaccinesLeftToDistibute,upperBound)

spreadersAllocation=Nadult*0;
ix=1;
while vaccinesLeftToDistibute>0 & ix<numel(adultAges)
    jx=adultAges(ix);
    currVacDist=min([s(jx).*Nadult(ix),upperBound(ix).*Nadult(ix),vaccinesLeftToDistibute]);
    spreadersAllocation(ix)=currVacDist./Nadult(ix);
    vaccinesLeftToDistibute=vaccinesLeftToDistibute-currVacDist;
    ix=ix+1;
end
return

