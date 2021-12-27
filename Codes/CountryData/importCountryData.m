countryList={"BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"}
for ix=1:numel(countryList)
    contactMatrix = readmatrix(join([countryList{ix},'contactMatrix.csv'],''));
    contactMatrix(:,1)=[];
    
    agDist = readmatrix(join([countryList{ix},'agDist.csv'],''));
    N=agDist(10,2)
    agDist=agDist(1:9,2);
    
   [V,d]=eig(contactMatrix );contactMatrix =contactMatrix /max(d(:));
   %susVec=ones(1,9);susVec(1)=0.7;
   contactMatrix=contactMatrix*diag(susVec);
   IFR=[9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01,4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01]/100;
   YLL_vec=[76.296132 66.912821 57.170992 47.969662 38.675122 29.722525 21.803262 14.672852  7.771498];

   save(join([countryList{ix},'_data.mat'],''),'contactMatrix','N','agDist','IFR','YLL_vec')
end