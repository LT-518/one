function ceps =ccostEps_ex4(ep,center,coll,bdpt_s,nb,lambda,exact,exact_nu,exact_L,fuds_1,fuds_2,fuds_3,...
    fuds_1laplace,fuds_2laplace,fuds_3laplace,dfuds_1,dfuds_2,dfuds_3,Nor)
source1=ep*bdpt_s+center;
source=[coll;source1];
DM_cacu=pdist2(coll,source);
%%
D_BKM_1=fuds_1(DM_cacu,lambda);
D_BKM_2=fuds_2(DM_cacu,lambda);
D_BKM_3=fuds_3(DM_cacu,lambda);
D_BKM=[D_BKM_1,D_BKM_2,D_BKM_3];
D_BKM(isnan(D_BKM))=0;

L_BKM_1=fuds_1laplace(DM_cacu,lambda);
L_BKM_2=fuds_2laplace(DM_cacu,lambda);
L_BKM_3=fuds_3laplace(DM_cacu,lambda);
L_BKM=[L_BKM_1,L_BKM_2,L_BKM_3];
L_BKM(isnan(L_BKM))=0;
%%
DMx2=coll(:,1)-repmat(source(:,1)',nb,1); %DMx2=DifferenceMatrix(xb,xs); %  samething
DMy2=coll(:,2)-repmat(source(:,2)',nb,1); %DifferenceMatrix(yb,ys);
DMz2=coll(:,3)-repmat(source(:,3)',nb,1); %DifferenceMatrix(zb,zs);

Hx_1=dfuds_1(DM_cacu,lambda,DMx2); Hy_1=dfuds_1(DM_cacu,lambda,DMy2); Hz_1=dfuds_1(DM_cacu,lambda,DMz2);
Hx_2=dfuds_2(DM_cacu,lambda,DMx2); Hy_2=dfuds_2(DM_cacu,lambda,DMy2); Hz_2=dfuds_2(DM_cacu,lambda,DMz2);
Hx_3=dfuds_3(DM_cacu,lambda,DMx2); Hy_3=dfuds_3(DM_cacu,lambda,DMy2); Hz_3=dfuds_3(DM_cacu,lambda,DMz2);

N_BKM_1=Hx_1.*Nor(:,1) + Hy_1.*Nor(:,2)+Hz_1.*Nor(:,3);
N_BKM_2=Hx_2.*Nor(:,1) + Hy_2.*Nor(:,2)+Hz_2.*Nor(:,3);
N_BKM_3=Hx_3.*Nor(:,1) + Hy_3.*Nor(:,2)+Hz_3.*Nor(:,3);
N_BKM=[N_BKM_1,N_BKM_2,N_BKM_3];
N_BKM(isnan(N_BKM))=0;

A=[D_BKM;N_BKM;L_BKM];
A(isnan(A))=0;

%%
D_BC=exact(coll(:,1),coll(:,2),coll(:,3));%Dirichlet boundary conditions
N_BC=exact_nu(coll(:,1),coll(:,2),coll(:,3),Nor(:,1),Nor(:,2),Nor(:,3));%Neumann boundary conditions
L_BC=exact_L(coll(:,1),coll(:,2),coll(:,3));%%\Delta u(x,y)
rhs=[D_BC;N_BC;L_BC];%boundary conditions
%% loocv
% invA=pinv(A);
% errorvector= (invA*rhs)./diag(invA);
% ceps=norm(errorvector);

%% Effective condition Method
 % S=svd(A); s_min=S(end);
   s_min=1;
coef=lsqminnorm(A,rhs); %
eff=1/s_min*norm(rhs)/norm(coef);
ceps=1/eff;

end