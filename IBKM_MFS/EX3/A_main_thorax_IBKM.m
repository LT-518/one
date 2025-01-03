clear; warning off all; tic 
%% 
source_type=deal(1);  % 0_coll, 1_ghost_loop, 2_loocv
[r1,r2]=deal(35,500);
j1=5;

for r=60  %r1:j1:r2
    radius_s=r/10;
%%   
lambda=sqrt(3);
u= @(x,y,z) cos(x+y).*sin(z); 
fuds=@(r,lambda)  1/(4*pi).*sin(lambda*r)./r;  %h
%%
NB1 = 600;  
NBH = 120; 

NT1 =1000;  
NTH =300;     
%%
[collt,collv,colllu,collru,coll,tspt]=Thorax_twolung_ventricle(NB1,NBH);
nb=size(coll,1);
nt=size(tspt,1);
coll=0.01*coll;
tspt=0.01*tspt;

xb=coll(:,1);yb=coll(:,2);zb=coll(:,3);
xt=tspt(:,1);yt=tspt(:,2);zt=tspt(:,3);
colll=[colllu;collru]; 
%% plot domain
%{
figure (1)
scatter1=scatter3(0.01*collt(:,1),0.01*collt(:,2),0.01*collt(:,3),8,'filled','MarkerFaceColor',[0.78431,0.63529,0.78431]) ;
scatter1.MarkerFaceAlpha = 0.3;  %% 0.3
scatter1.MarkerEdgeAlpha = 0.3;
hold on;                                                                
scatter3(0.01*collv(:,1),0.01*collv(:,2),0.01*collv(:,3),5,'filled', 'MarkerFaceColor',[0.8500,0.3250,0.0980]) 
scatter3(0.01*colll(:,1),0.01*colll(:,2),0.01*colll(:,3),5,'filled','MarkerFaceColor',[0.4660,0.6740,0.1880])  
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on;
%}
%% Ghost point
ns=nb;
bmin=min(coll,[],1);bmax=max(coll,[],1);
center=(bmin+bmax)/2;
noise=0.1;
psr=-1+2*rand(nb,1); 
switch source_type
    case 0
    source=coll;
    xs=source(:,1);ys=source(:,2);zs=source(:,3);
   case 1 
    source=getSolidSphere_3D(ns,radius_s,coll,center);
    xs=source(:,1).*(1+psr*noise);ys=source(:,2).*(1+psr*noise);zs=source(:,3).*(1+psr*noise);
    case 2      
    LL=3; RR= 16;  
    radius_s=1;
    source1=getSolidSphere_3D(ns,radius_s,coll,center);
    option=optimset('TolX',1e-5,'MaxIter',25);
    radius_s=fminbnd(@(ep)costEps_ex3(ep,coll,source1,lambda,u,fuds,psr,noise),LL,RR,option);

    xs=source1(:,1)*radius_s.*(1+psr*noise);
    ys=source1(:,2)*radius_s.*(1+psr*noise);
    zs=source1(:,3)*radius_s.*(1+psr*noise);    
end
%% plot ghost
%{
scatter3(xs,ys,zs,2,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.9)
legend({'thorax','ventricle','lung','ghost point'},'FontSize',12,'Fontname', 'Times New Roman')
hold on
grid of
xlabel('X');
ylabel('Y');
zlabel('Z');
%}
%% 
DM=DistanceMatrix([xb yb zb],[xs ys zs]); 

BKM=fuds(DM,lambda);
BKM(isnan(BKM))=0;
%% RHS
BC=u(xb,yb,zb);
coef=lsqminnorm(BKM,BC);
%%
ut= u(xt,yt,zt);  
DM_tst=DistanceMatrix([xt yt zt],[xs ys zs]);
BKMt=fuds(DM_tst,lambda);
BKMt(isnan(BKMt))=0; 

approx_tst = BKMt*coef;   % approximate solution
rmse=norm(ut-approx_tst,2)/sqrt(nt);
max_err=max(abs(ut-approx_tst));

% fprintf('LL=%6.2f,RR=%6.2f\n',LL,RR)
fprintf('NB1=%3d,NBH=%3d,nb=%3d,nt=%3d, rmse= %8.4e, radius_s= %4.4f\n ',NB1,NBH,nb,nt,rmse,radius_s)
end
toc
fprintf('====================================\n')