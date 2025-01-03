function ceps =costEps_ex1(ep,coll_sph,NBH,nort1,nort2,up_coll,lo_coll,l_coll,r_coll,other_coll,...
    y_oth,x_oth,xb_up,xb_lo,lambda,fuds,fuds_r)
s_sph1=coll_sph(1:NBH,:)    +   ep.*nort1;
s_sph2=coll_sph(NBH+1:end,:)+   ep.*nort2;

s_up=up_coll+ep*[0,1];
s_lo=lo_coll(2:end-1,:)+ep*[0,-1];
s_l=l_coll(1:end-1,:)+ep*[-1,0];
s_r=r_coll(1:end-1,:)+ep*[1,0];
a_llo=[0.1,0.1]+ep*[-sqrt(2)/2,-sqrt(2)/2];
a_lup=[0.1,1.6]+ep*[-sqrt(2)/2,sqrt(2)/2];
a_rup=[1.1,1.6]+ep*[sqrt(2)/2,sqrt(2)/2];
a_rlo=[1.1,0.1]+ep*[sqrt(2)/2,-sqrt(2)/2];
source=[s_sph1;s_sph2;s_up;s_l;s_lo;s_r;a_llo;a_lup;a_rup;a_rlo];

%%coll norm
ns=length(source);
NBE=length(y_oth);
M=source(:,1);N=source(:,2);
m1=repmat(y_oth,1,ns);m2=repmat(N',NBE,1);
n1=repmat(x_oth,1,ns);n2=repmat(M',NBE,1);
m3=m1-m2; n3=n1-n2;
dx=gradient(x_oth);dy=gradient(y_oth);
nx=dy./sqrt(dx.^2+dy.^2);ny=-dx./sqrt(dx.^2+dy.^2);
nx1=repmat(nx,1,ns);ny1=repmat(ny,1,ns);
nor=n3.*nx1+m3.*ny1;
%% Calculate the approximate solution
% RHS
BC_1=ones(length(xb_up), 1);
BC_2=0*ones(length(xb_lo), 1);
BC_3=0*ones(length(x_oth), 1);
BC=[BC_1;BC_2;BC_3];
%LHS
DM_cacu_1=DistanceMatrix(up_coll,source);
D_BKM_1=fuds(DM_cacu_1,lambda);
DM_cacu_2=DistanceMatrix(lo_coll,source);
D_BKM_2=fuds(DM_cacu_2,lambda);
D_BKM=[D_BKM_1;D_BKM_2];
%coll
DM_cacu_3=DistanceMatrix(other_coll,source);
N_BKM=fuds_r(DM_cacu_3,lambda);

N_BKM(isnan(N_BKM))=0;
N_BKM=N_BKM.*nor;
BKM=[D_BKM;N_BKM];
%% Effective condition Method
% S=svd(BKM); s_min=S(end);
  s_min=1;
coef=lsqminnorm(BKM,BC);
eff=1/s_min*norm(BC)/norm(coef);
ceps=1/eff;

end