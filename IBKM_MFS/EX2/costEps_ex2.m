function ceps=costEps_ex2(ep,coll,source0,lambda,u,fuds,psr,noise,center)
source=source0.*ep.*(1+psr*noise)+center;
xs=source(:,1);ys=source(:,2);zs=source(:,3); 
DM=DistanceMatrix(coll,[xs ys zs]); 
BKM=fuds(DM,lambda);
BKM(isnan(BKM))=0;
BC=u(coll(:,1),coll(:,2),coll(:,3));
%% Effective condition Method
 S=svd(BKM); s_min=S(end);
%s_min=1;
coef=lsqminnorm(BKM,BC); 
eff=1/s_min*norm(BC)/norm(coef);
ceps=1/eff;
end