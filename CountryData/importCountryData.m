countryList={"USA","US2","ISR","BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"}
for ix=1:numel(countryList)
    contactMatrix = readmatrix(join([countryList{ix},'contactMatrix.csv'],''));
    contactMatrix(:,1)=[];if isnan(contactMatrix(1,1)) contactMatrix(1,:)=[];end
    
    agDist = readmatrix(join([countryList{ix},'agDist.csv'],''));
    N=agDist(11,2)
    agDist=agDist(2:10,2);
    
    susVec=[0.4, 0.38, 0.79,0.86,0.8, 0.82,0.88,0.74,0.74];
    contactMatrix=diag(susVec)*contactMatrix;
    Mij=diag(agDist)*contactMatrix*diag(1./agDist);
    [V,d]=eig(Mij);contactMatrix =contactMatrix/max(d(:));
    
    IFR=[9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01,4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01]/100;
    YLL_vec=[76.296132 66.912821 57.170992 47.969662 38.675122 29.722525 21.803262 14.672852  7.771498];
    
    save(join([countryList{ix},'_data.mat'],''),'contactMatrix','N','agDist','IFR','YLL_vec')
end