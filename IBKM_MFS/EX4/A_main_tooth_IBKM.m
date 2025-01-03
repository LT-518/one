warning off all; 
clear all;
close all;
tic
%%
source_type=deal(2);  %  0_coll, 1_ghost_loop 2_cloud
[r1,r2]=deal(3,30);  
j1=1;

for r=18   %r1:j1:r2
    radius_s=r;
%%
xi1=2;
BigN=3;
lambda=1;
%% 
exact=@(x,y,z) exp(x)-2*exp(y)+exp(z); 
exact_nu=@(x,y,z,nx,ny,nz) exp(x).* nx - 2*exp(y).* ny + exp(z).* nz;   
exact_L=@(x,y,z) exp(x)-2*exp(y)+exp(z);
%%
fuds_1=@(r,lambda) besseli(1/2,lambda*r)./(lambda*r).^(1/2);
fuds_2=@(r,lambda) besseli(3/2,lambda*r).*(r).^(1/2)./(lambda).^(3/2);
fuds_3=@(r,lambda) besseli(5/2,lambda*r).*(r).^(3/2)./(lambda).^(5/2);

fuds_1r=@(r,lambda) sqrt(0.2e1).*(cosh(lambda.*r).*lambda.*r - sinh(lambda.*r))./ lambda./ r.^2.*pi.^(-0.1e1./0.2e1);
fuds_2r=@(r,lambda) sqrt(0.2e1).*((lambda.*r).^(-0.1e1./0.2e1)).*(lambda.^(-0.5e1./0.2e1)).*(r.^(-0.3e1./0.2e1)).*...
    ((lambda.^ 2.*r.^ 2 + 1).*sinh((lambda.*r))-cosh((lambda.*r)).*lambda.*r).*pi.^(-0.1e1./0.2e1);
fuds_3r=@(r,lambda) sqrt(0.2e1).*(cosh(lambda.*r).*lambda.^ 3.*r.^ 3 - 0.2e1.*lambda.^ 2.*sinh(lambda.*r).*r.^ 2 +...
    0.3e1.*cosh(lambda.*r).*lambda.*r - 0.3e1.* sinh(lambda.*r)).*lambda.^ (-0.9e1./ 0.2e1).* r.^ (-0.3e1./ 0.2e1).*...
    (lambda.* r).^(-0.1e1./0.2e1).* pi.^ (-0.1e1./ 0.2e1);

fuds_1laplace=@(r,lambda) 0.1e1./ r.*sqrt(0.2e1).* pi.^(-0.1e1./0.2e1).*lambda.*sinh(lambda.*r);
fuds_2laplace=@(r,lambda) r.^(-0.1e1./0.2e1).*sqrt(0.2e1).*lambda.^(-0.1e1./0.2e1).*(cosh(lambda.*r).*lambda.* r + ...
    sinh(lambda.*r)).*pi.^(-0.1e1./0.2e1).*(lambda.*r).^(-0.1e1./0.2e1);
fuds_3laplace=@(r,lambda) sqrt(0.2e1).* ((lambda.* r).^ (-0.1e1./ 0.2e1)).* (lambda.^ (-0.5e1./ 0.2e1)).* (r.^ (-0.1e1./ 0.2e1)).*...
    pi.^ (-0.1e1./ 0.2e1).* ((lambda.^ 2.* r.^ 2 - 1).* sinh((lambda.* r)) + cosh((lambda.* r)).* lambda.* r);

dR_x=@(r,lambda,dd) dd./r;
dfuds_1=@(r,lambda,dd) fuds_1r(r,lambda).*dR_x(r,lambda,dd);
dfuds_2=@(r,lambda,dd) fuds_2r(r,lambda).*dR_x(r,lambda,dd);
dfuds_3=@(r,lambda,dd) fuds_3r(r,lambda).*dR_x(r,lambda,dd);
%% Generate collocation points
NB1 = 800;   NBH = 60;   
NT1 =1000;  NTH =300;  

[collc,colls,colls1,colls2,colls3,coll,tspt]=IBKM_tooth_hollow(NB1,NBH,NTH);
nt=size(tspt,1);
nb=size(coll,1);
%%
T = delaunay(coll(:,1),coll(:,2)); 
TR = triangulation(T,coll(:,1),coll(:,2),coll(:,3)); 
Nor= vertexNormal(TR); 
%% source or ghost
noise=0.1;
psr=-1+2*rand(nb,1);
bmin=min(coll,[],1);bmax=max(coll,[],1);  
center=(bmin+bmax)/2;
 switch source_type
   case 0
      source=coll;
      radius_s=1.0;
    case 1   
      ns=nb; 
      [bdpt_s,~]=generateB_LEI(ns,1,'sphere');
      %%%ecn
%       LL=2; RR=20;
%       option=optimset('TolX',1e-4);
%       option=optimset('TolX',1e-30,'MaxIter',10);  %%,'PlotFcns',@optimplotfval,'Display','iter'
%       radius_s=fminbnd(@(ep)costEps_ex4(ep,center,coll,bdpt_s,psr,noise,nb,lambda,exact,exact_nu,exact_L,fuds_1,...
%           fuds_2,fuds_3,fuds_1laplace,fuds_2laplace,fuds_3laplace,dfuds_1,dfuds_2,dfuds_3,Nor),LL,RR,option);      
      source1=radius_s*bdpt_s'+center;   
      source2=xi1*radius_s*bdpt_s'+center;
      source=[coll;source1.*(1+psr*noise);source2.*(1+psr*noise)];
    case 2
         ns=2*nb;
      m=nb;
      [xs,ys,zs]=pss3D(m,ns);
      bdpt_s=[xs ys zs];
      %ecn
     % LL=2; RR=20;
     % option=optimset('TolX',1e-4,'MaxIter',9);  %%,'PlotFcns',@optimplotfval,'Display','iter'
     % radius_s=fminbnd(@(ep)ccostEps_ex4(ep,center,coll,bdpt_s,nb,lambda,exact,exact_nu,exact_L,fuds_1,...
     %      fuds_2,fuds_3,fuds_1laplace,fuds_2laplace,fuds_3laplace,dfuds_1,dfuds_2,dfuds_3,Nor),LL,RR,option);
      source1=radius_s*bdpt_s+center;   
      source=[coll;source1];
 end
%% ghost
%{
figure(2)  
scatter3(colls(:,1),colls(:,2),colls(:,3),5,'fill','bo');
hold on;
scatter3(source1(:,1),source1(:,2),source1(:,3),5,'filled','MarkerFaceColor',[0,1,0],'MarkerFaceAlpha',0.9);
% scatter3(source2(:,1),source2(:,2),source2(:,3),5,'filled','MarkerFaceColor',[1,0,0],'MarkerFaceAlpha',0.8);
scatter3=scatter3(collc(:,1),collc(:,2),collc(:,3),4,'filled','bo');
scatter3.MarkerFaceAlpha = 0.2;
hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid off
hold on;
%}
%% Calculate the approximate solution
D_BC=exact(coll(:,1),coll(:,2),coll(:,3));
N_BC=exact_nu(coll(:,1),coll(:,2),coll(:,3),Nor(:,1),Nor(:,2),Nor(:,3));%Neumann boundary conditions
L_BC=exact_L(coll(:,1),coll(:,2),coll(:,3));%%\Delta u(x,y)
BC=[D_BC;N_BC;L_BC];
%LHS
DM_cacu=pdist2(coll,source);
DM_tst=pdist2(tspt,source);

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
DMx2=coll(:,1)-repmat(source(:,1)',nb,1); %DMx2=DifferenceMatrix(xb,xs); %samething
DMy2=coll(:,2)-repmat(source(:,2)',nb,1);
DMz2=coll(:,3)-repmat(source(:,3)',nb,1); 

Hx_1=dfuds_1(DM_cacu,lambda,DMx2); Hy_1=dfuds_1(DM_cacu,lambda,DMy2); Hz_1=dfuds_1(DM_cacu,lambda,DMz2);
Hx_2=dfuds_2(DM_cacu,lambda,DMx2); Hy_2=dfuds_2(DM_cacu,lambda,DMy2); Hz_2=dfuds_2(DM_cacu,lambda,DMz2);
Hx_3=dfuds_3(DM_cacu,lambda,DMx2); Hy_3=dfuds_3(DM_cacu,lambda,DMy2); Hz_3=dfuds_3(DM_cacu,lambda,DMz2);

N_BKM_1=Hx_1.*Nor(:,1) + Hy_1.*Nor(:,2)+Hz_1.*Nor(:,3);
N_BKM_2=Hx_2.*Nor(:,1) + Hy_2.*Nor(:,2)+Hz_2.*Nor(:,3);
N_BKM_3=Hx_3.*Nor(:,1) + Hy_3.*Nor(:,2)+Hz_3.*Nor(:,3);
N_BKM=[N_BKM_1,N_BKM_2,N_BKM_3];
N_BKM(isnan(N_BKM))=0;

BKM=[D_BKM;N_BKM;L_BKM];
BKM(isnan(BKM))=0;

coef=lsqminnorm(BKM,BC);
%% error
exact_tst=exact(tspt(:,1),tspt(:,2),tspt(:,3));
D_BKMt=[fuds_1(DM_tst,lambda),fuds_2(DM_tst,lambda),fuds_3(DM_tst,lambda)];
D_BKMt(isnan(D_BKMt))=0;
approx_tst=D_BKMt*coef;
RMSE=norm(abs(exact_tst-approx_tst),2)/sqrt(nt);

% fprintf('LL=%6.2f,RR=%6.2f\n',LL,RR)
fprintf('NB1= %3d,NBH= %3d,nb= %3d,nt= %3d，',NB1,NBH,nb,nt);
fprintf('RMSE = %8.4e, radius_s =%6.2f\n',RMSE,radius_s)
end
toc
fprintf('====================================\n')