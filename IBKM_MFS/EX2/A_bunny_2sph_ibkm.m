clear; warning off all; tic 
%%
source_type=deal(2);  %0_coll, 1_ghost_loop, 2_ecn
[r1,r2]=deal(80,200);
j1=1;
load('bunny_2sphere.mat')

for r=80   % r1:j1:r2
    radius_s=r/10;
%%
lambda=1;
fuds=@(r,lambda)  1/(4*pi).*sinh(lambda*r)./r;%mh
u= @(x,y,z) sin(y+z).*cosh(sqrt(3).*x)-cos(x+z).*sinh(sqrt(3).*y)+exp(sqrt(3).*z).*sin(x+y)+1/exp(3*pi);

NB1= 400;   
NB2= 80;   %80,40 
NB=NB1+2*NB2;

[collb]=sel_point(bunny,NB1);   
[colls_1]=sel_point(sphere1,NB2);  %r=0.2
[colls_2]=sel_point(sphere2,NB2);  
colls=[colls_1;colls_2];
coll=[collb;colls]; 
xb=coll(:,1);yb=coll(:,2);zb=coll(:,3);
%%
ns=length(coll);
bmin=min(coll,[],1);bmax=max(coll,[],1);
center=(bmin+bmax)/2;
noise=0.1;
psr=-1+2*rand(NB,1);
switch source_type
    case 0
    source=coll;
    xs=source(:,1);ys=source(:,2);zs=source(:,3);
   case 1 
    source=getSolidSphere_3D(ns,radius_s,coll,center);
    xs=source(:,1).*(1+psr*noise);ys=source(:,2).*(1+psr*noise);zs=source(:,3).*(1+psr*noise);  %.*(1+psr*noise)
    case 2
    LL=3; RR=20;
    radius_s=1;
    source0=getSolidSphere_3D(ns,radius_s,coll,[0,0,0]);
    radius_s=fminbnd(@(ep)costEps_ex2(ep,coll,source0,lambda,u,fuds,psr,noise,center),LL,RR,optimset('TolX',1e-4,'MaxIter',15));  %fix
    source=source0.*radius_s.*(1+psr*noise)+center;
    xs=source(:,1);ys=source(:,2);zs=source(:,3);   
end

%{
 C=linspecer(8);
figure(1)
scatter3(collb(:,1),collb(:,2),collb(:,3),3,'filled','MarkerFaceColor',C(7,:),'MarkerFaceAlpha',0.5) ; 
hold on
C=linspecer(5);
scatter3(colls(:,1),colls(:,2),colls(:,3),4,'filled', 'MarkerFaceColor',C(2,:),'MarkerFaceAlpha',0.7) 
hold on
scatter3(xs(1:3:end,:),ys(1:3:end,:),zs(1:3:end,:),2,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.9) ; 
xlabel('X');
ylabel('Y');
zlabel('Z');
grid of
%}
%%
DM=DistanceMatrix(coll,[xs ys zs]); 
BKM=fuds(DM,lambda);
BKM(isnan(BKM))=0;
BC=u(xb,yb,zb);
% coef=BKM\BC;
coef=lsqminnorm(BKM,BC);
%%
Nt=1500; 
load('u_FEM.mat')
t_sq=randperm(size(data, 1), Nt);
sel_t = data(t_sq, :);
tspt=sel_t(:,1:3);
xt=tspt(:,1);yt=tspt(:,2);zt=tspt(:,3);
ut = sel_t(:,4);
NT=length(ut);

DM_tst=DistanceMatrix(tspt,[xs ys zs]);
BKMt=fuds(DM_tst,lambda);
BKMt(isnan(BKMt))=0; 

approx_tst = BKMt*coef;
rmse=norm(ut-approx_tst,2)/sqrt(NT);

% fprintf('LL=%6.2f,RR=%6.2f\n',LL,RR)
fprintf('NB1= %3d,NB2= %3d,NT= %3d, rmse = %8.4e,radius_s = %6.2f\n',NB1,NB2,NT,rmse,radius_s);
end
toc