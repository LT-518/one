clear all; warning off all; tic 
%%  
load('bunny_2sphere.mat')

%%%%%r(max)=sqrt(3)=1.732
%{
     radius_us=linspace(1.8,10,83);   %0.1   
for  radius_u0=1:size(radius_us,2)
     R1=radius_us(radius_u0);

      radius_is=linspace(0.001,0.2,200); %0.001
for   radius_i0=1:size(radius_is,2)
      R2=radius_is(radius_i0);
 %}
R1=4.5;R2=0.008;
lambda=1;  
fuds=@(r,lambda)  1/(4*pi).*exp(-lambda*r)./r;%mh
u= @(x,y,z) sin(y+z).*cosh(sqrt(3).*x)-cos(x+z).*sinh(sqrt(3).*y)+exp(sqrt(3).*z).*sin(x+y)+1/exp(3*pi); 

NB1=300; 
NB2=40;    
NB=NB1+2*NB2;

[collb]=sel_point(bunny,NB1);   
[colls_1]=sel_point(sphere1,NB2);  %r=0.2
[colls_2]=sel_point(sphere2,NB2);  
colls=[colls_1;colls_2];
coll=[collb;colls]; 
xb=coll(:,1);yb=coll(:,2);zb=coll(:,3);
%%
bmin=min(coll,[],1);bmax=max(coll,[],1);
center=(bmin+bmax)/2;
center2=[0.8,0.8,1.1]; center3=[1.3,0.8,1.1];

s_ucm=RandSampleSphere(NB1,'uniform');  
s_uc=s_ucm*R1+center;

[s_s,~]=generateB_LEI(NB2,1,'sphere');
s_1=R2*s_s'+center2;
s_2=R2*s_s'+center3;
s_ic=[s_1;s_2];
s=[s_uc;s_ic];

%{
 C=linspecer(8);
figure(1)
scatter3(collb(:,1),collb(:,2),collb(:,3),3,'filled','MarkerFaceColor',C(7,:),'MarkerFaceAlpha',0.5) ; 
hold on
C=linspecer(5);
scatter3(colls(:,1),colls(:,2),colls(:,3),4,'filled', 'MarkerFaceColor',C(2,:),'MarkerFaceAlpha',0.6) 
hold on
scatter3(s_uc(:,1),s_uc(:,2),s_uc(:,3),3,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.6) ;
hold on;
scatter3(s_ic(:,1),s_ic(:,2),s_ic(:,3),3,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.9) ;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid of
%}
%%
DM=DistanceMatrix(coll,s); 
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

DM_tst=DistanceMatrix(tspt,s);
BKMt=fuds(DM_tst,lambda);
BKMt(isnan(BKMt))=0; 
approx_tst = BKMt*coef;
rmse=norm(ut-approx_tst,2)/sqrt(NT);

% err(radius_u0+1,1)=R1;
% err(1,radius_i0+1)=R2;
% err(radius_u0+1,radius_i0+1)=rmse;
fprintf('NB1= %3d,NB2= %3d,NT= %3d, rmse = %8.4e, R1= %4.3f,R2= %4.3f\n',NB1,NB2,NT,rmse,R1,R2);
% end
% end
toc
fprintf('====================================\n')