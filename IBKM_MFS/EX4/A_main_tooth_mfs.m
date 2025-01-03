warning off all; 
clear all;
close all;
tic
%%
%{
R1i=[];R2i=[];R3i=[];R4i=[];RMSEi=[];
     radius_ts=linspace(3.5,5.5,3);   %%0.5  r1>2 2-
for  radius_t0=1:size(radius_ts,2)
     R1=radius_ts(radius_t0);

      radius_as=linspace(0.01,0.20,20); %%0.01  r2=0.3
for   radius_a0=1:size(radius_as,2)
      R2=radius_as(radius_a0);

     radius_bs=linspace(0.01,0.13,13);   %0.01; r3=0.2
for  radius_b0=1:size(radius_bs,2)  
     R3=radius_bs(radius_b0);

     radius_cs=linspace(0.01,0.07,7);   %0.01
for  radius_c0=1:size(radius_cs,2)  
     R4=radius_cs(radius_c0);  
%}
R1=3.5;
R2=0.18;R3=0.12; R4=0.05; 
%% 
lambda=1;
exact=@(x,y,z) exp(x)-2*exp(y)+exp(z);  
exact_nu=@(x,y,z,nx,ny,nz) exp(x).* nx - 2*exp(y).* ny + exp(z).* nz;   
exact_L=@(x,y,z) exp(x)-2*exp(y)+exp(z);
%%
fuds_1=@(r,lambda) 1/(4*pi).*exp(-lambda*r)./r;
fuds_2=@(r,lambda) 1/(8*pi).*exp(-lambda*r)./lambda;
fuds_3=@(r,lambda) (1 + r.* lambda).* exp(-(r * lambda))./ pi.* (lambda.^ 3)./0.32e2;

fuds_1r=@(r,lambda) -lambda.* exp(-lambda.* r)./ r./ pi./ 0.4e1 - exp(-lambda.* r)./ r.^ 2./ pi./ 0.4e1;
fuds_2r=@(r,lambda) -exp(-lambda.* r)./ pi./ 0.8e1;
fuds_3r=@(r,lambda) -exp(-lambda.* r)./ lambda.* r./ pi./ 0.32e2;

fuds_1laplace=@(r,lambda) lambda.^ 2.* exp(-lambda.* r)./ pi./ r./ 0.4e1;
fuds_2laplace=@(r,lambda) 0.1e1./ r.* exp(-lambda.* r).* (lambda.* r - 0.2e1)./ pi./ 0.8e1;
fuds_3laplace=@(r,lambda) exp(-lambda.* r)./ lambda.* (lambda.* r - 0.3e1)./ pi./ 0.32e2;

dR_x=@(r,lambda,dd) dd./r;
dfuds_1=@(r,lambda,dd) fuds_1r(r,lambda).*dR_x(r,lambda,dd);
dfuds_2=@(r,lambda,dd) fuds_2r(r,lambda).*dR_x(r,lambda,dd);
dfuds_3=@(r,lambda,dd) fuds_3r(r,lambda).*dR_x(r,lambda,dd);
%% Generate collocation points
NB1 =800;   NBH =60;    
NT1 =1000;  NTH =300; 
%%%%%%
[collc,colls,colls1,colls2,colls3,coll,tspt]=IBKM_tooth_hollow(NB1,NBH,NTH);
nt=size(tspt,1);
nb=size(coll,1);  
%%
T = delaunay(coll(:,1),coll(:,2)); 
TR = triangulation(T,coll(:,1),coll(:,2),coll(:,3)); 
Nor= vertexNormal(TR); 
%% outsider ghost
s_ucm =RandSampleSphere(NB1,'uniform');  
bmin=min(coll,[],1);bmax=max(coll,[],1);
center=(bmin+bmax)/2;
s_t=s_ucm*R1+center;
%%inn ghost
s_icm =RandSampleSphere(NBH,'uniform');
s_a =s_icm*R2+[0,0,0.4]; 
s_b =s_icm*R3+[0.3,0.4,0.1];
s_c =s_icm*R4+[-0.4,-0.4,0];
s=[s_t;s_a;s_b;s_c];
%% polt source
%{
figure (1)
scatter1=scatter3(collc(:,1),collc(:,2),collc(:,3),5,'filled','MarkerFaceColor',[0.78431,0.63529,0.78431]) ;
scatter1.MarkerFaceAlpha = 0.6; 
scatter1.MarkerEdgeAlpha = 0.6;
hold on;  
scatter3(colls1(:,1),colls1(:,2)-0.13,colls1(:,3),5,'filled', 'MarkerFaceColor',[0.8500,0.3250,0.0980],'MarkerFaceAlpha',0.6) 
scatter3(colls2(:,1),colls2(:,2),colls2(:,3),5,'filled','MarkerFaceColor',[0.4660,0.6740,0.1880],'MarkerFaceAlpha',0.6)
scatter3(colls3(:,1),colls3(:,2),colls3(:,3),3,'filled','MarkerFaceColor',[0.9290,0.6940,0.1250],'MarkerFaceAlpha',0.6)

scatter3(s_a(:,1),s_a(:,2)-0.13,s_a(:,3),5,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.6)
scatter3(s_t(:,1),s_t(:,2),s_t(:,3),5,'filled', 'MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.2)
scatter3(s_b(:,1),s_b(:,2),s_b(:,3),5,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.5)
scatter3(s_c(:,1),s_c(:,2),s_c(:,3),5,'filled','MarkerFaceColor',[0,0,1],'MarkerFaceAlpha',0.5)
legend({'tooth','cavity I','cavity II','cavity III','ghost point'},'FontSize',12,'Fontname', 'Times New Roman')
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on;
%}
%% Calculate the approximate solution
D_BC=exact(coll(:,1),coll(:,2),coll(:,3));
N_BC=exact_nu(coll(:,1),coll(:,2),coll(:,3),Nor(:,1),Nor(:,2),Nor(:,3));%Neumann boundary conditions
L_BC=exact_L(coll(:,1),coll(:,2),coll(:,3));
BC=[D_BC;N_BC;L_BC];
%LHS
DM_cacu=pdist2(coll,s);
DM_tst=pdist2(tspt,s);

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
DMx2=coll(:,1)-repmat(s(:,1)',nb,1); 
DMy2=coll(:,2)-repmat(s(:,2)',nb,1); 
DMz2=coll(:,3)-repmat(s(:,3)',nb,1); 

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
%%
% R1i=[R1i;R1];R2i=[R2i;R2];R3i=[R3i;R3];R4i=[R4i;R4];RMSEi=[RMSEi;RMSE];
% save('mfs_err.mat','R1i','R2i','R3i','R4i','RMSEi'); 
fprintf('NB1=%3d,NBH=%3d,nb=%3d,nt=%3d,RMSE= %8.4e,R1= %4.2f,R2= %4.2f,R3= %4.2f,R4= %4.2f\n ',NB1,NBH,nb,nt,RMSE,R1,R2,R3,R4)

%{
end
end
end
end
%}
toc