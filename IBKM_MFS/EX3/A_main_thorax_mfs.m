% Example 6:  Trihamonic 3D airplane
clear; warning off all; tic 
%%

%{
     radius_ts=linspace(3,10,71);   %%0.1
for  radius_t0=1:size(radius_ts,2)
     R1=radius_ts(radius_t0);

      radius_ls=linspace(0.01,0.6,60); %%0.01
for   radius_l0=1:size(radius_ls,2)
      R2=radius_ls(radius_l0);
           
     radius_vs=linspace(0.01,1,100);   %0.01
for  radius_v0=1:size(radius_vs,2)  
     R3=radius_vs(radius_v0);
     
%}
R1=6.5; 
R2=0.08; R3=0.02; 
%% 
lambda=sqrt(3);
u= @(x,y,z) cos(x+y).*sin(z);
%%
fuds=@(r,lambda)  1/(4*pi).*cos(lambda*r)./r; 
%%
NB1 = 800; NBH = 90;   %90/120
NT1 =1000;  NTH =300;     
%%
[collt,collv,colllu,collru,coll,tspt,nort,norv,norlu,norru,nor]=Thorax_twolung_ventricle(NB1,NBH);
nb=size(coll,1);
nt=size(tspt,1);
coll=0.01*coll;  collt=0.01*collt;  collv=0.01*collv; colllu=0.01*colllu; collru=0.01*collru;
tspt=0.01*tspt;  

xb=coll(:,1);yb=coll(:,2);zb=coll(:,3);
xt=tspt(:,1);yt=tspt(:,2);zt=tspt(:,3);
%% Ghost point ns=NB
s_t =collt - nort.*R1;  %%-out
s_lu =colllu + norlu.*R2;  %%+in
s_ru =collru + norru.*R2;
s_v =collv + norv.*R3;
s=[s_t;s_lu;s_ru;s_v];
%% polt ghost
s_u=[s_lu;s_ru];
collu=[colllu;collru];
collv=collv+[0,0,-0.45];
s_v=s_v+[0,0,-0.45];

%%plot domain and source
%{
scatter1=scatter3(collt(:,1),collt(:,2),collt(:,3),8,'filled','MarkerFaceColor',[0.78431,0.63529,0.78431]) ;
scatter1.MarkerFaceAlpha = 0.3;  
scatter1.MarkerEdgeAlpha = 0.3;
hold on;                                                               
scatter3(collv(:,1),collv(:,2),collv(:,3),5,'filled', 'MarkerFaceColor',[0.8500,0.3250,0.0980]) 
scatter3(collu(:,1),collu(:,2),collu(:,3),5,'filled','MarkerFaceColor',[0.4660,0.6740,0.1880],'MarkerFaceAlpha',0.8)  
%source
scatter3(s_t(1:3:end,1),s_t(1:3:end,2),s_t(1:3:end,3),5,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.4)
scatter3(s_v(:,1),s_v(:,2),s_v(:,3),5,'filled', 'MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.5)
scatter3(s_u(:,1),s_u(:,2),s_u(:,3),5,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.5)  
legend({'thorax','ventricle','lung','ghost point'},'FontSize',12,'Fontname', 'Times New Roman')
hold on
xlabel('X');    
ylabel('Y');
zlabel('Z');
% grid of
%}
%% 
DM=DistanceMatrix(coll,s); 
BKM=fuds(DM,lambda);
BKM(isnan(BKM))=0;
%% RHS
BC=u(xb,yb,zb);
coef=lsqminnorm(BKM,BC);
%%
ut= u(xt,yt,zt);
DM_tst=DistanceMatrix(tspt,s);
BKMt=fuds(DM_tst,lambda);
BKMt(isnan(BKMt))=0; 

approx_tst = BKMt*coef;
rmse=norm(ut-approx_tst,2)/sqrt(nt);

fprintf('NB1=%3d,NBH=%3d,nb=%3d,nt=%3d, rmse= %8.4e, R1= %4.3f,R2= %4.3f,R3= %4.3f\n ',NB1,NBH,nb,nt,rmse,R1,R2,R3)
%{
end
end
end
%}
toc
fprintf('====================================\n')