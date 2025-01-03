function ceps =costEps_ex3(ep,coll,source1,lambda,u,fuds,psr,noise)
% source=getSolidSphere_3D(ns,ep,coll,center);
xs=source1(:,1)*ep.*(1+psr*noise);  
ys=source1(:,2)*ep.*(1+psr*noise);
zs=source1(:,3)*ep.*(1+psr*noise);
xb=coll(:,1);yb=coll(:,2);zb=coll(:,3);
DM=DistanceMatrix([xb yb zb],[xs ys zs]);

BKM=fuds(DM,lambda);
BKM(isnan(BKM))=0;
%%RHS
BC=u(xb,yb,zb);
%% Effective condition Method
% S=svd(BKM); s_min=S(end);
  s_min=1;
coef=lsqminnorm(BKM,BC); %
eff=1/s_min*norm(BC)/norm(coef);
ceps=1/eff;
end