warning off all;
clear;
tic
%% 
source_type=deal(1);  %0_coll 1_GHOST
%{
      radius_us=linspace(0.01,1,100);    
for   radius_u0=1:size(radius_us,2)
      R=radius_us(radius_u0);
%}
R=0.06;
%% exact solution and RBF
lambda=1;
exact_1=@(x,y)  1;   
exact_2=@(x,y)  0;
exact_3=@(x,y)  0;

fuds=@(r,lambda)  besselh(0, lambda * r);   
fuds_r=@(r,lambda) -lambda*besselh(1, lambda*r)./r;          
%% Generate collocation points on boundaryn % square + two sphere
%%
NB1=200;   %%40,80,100,200
NBH=60;  NB2=2*NBH;   %%60,80
NB=NB1+NB2; 

[up_coll,lo_coll,lr_coll,l_coll,r_coll,coll_sph,coll] = generate_points1(NB1,NBH);
xb_up=up_coll(:,1);yb_up=up_coll(:,2);
xb_lo=lo_coll(:,1);yb_lo=lo_coll(:,2);
other_coll=[lr_coll;coll_sph];
x_oth=other_coll(:,1);y_oth=other_coll(:,2);
%% source or ghost
switch source_type
        case 0
            [up_s,lo_s,lr_s,l_s,r_s,s_sph,source] = generate_points1(60,61);
            DMM=DistanceMatrix(coll,source);
            [rows, cols] = find(DMM < 0.01);
            keep = true(size(source, 1), 1);
            keep(cols) = false;
            filtered_s = source(keep, :);
            source=filtered_s;
        case 1        
            vectors_1 = [0.50, 1.2]- coll_sph(1:NBH,:);   
            lengths_1 = sqrt(sum(vectors_1.^2, 2));  
            nort1 = vectors_1 ./ lengths_1(:, ones(1, 2));  
            
            vectors_2 = [0.70,0.50]- coll_sph(NBH+1:end,:);
            lengths_2 = sqrt(sum(vectors_2.^2, 2));
            nort2 = vectors_2 ./ lengths_2(:, ones(1, 2));

            % LL=0.01;RR=0.2; %%ecn and eecn
            % option=optimset('TolX',1e-4,'MaxIter',13);
            % R=fminbnd(@(ep)costEps_ex1(ep,coll_sph,NBH,nort1,nort2,up_coll,lo_coll,l_coll,r_coll,other_coll,...
            %     y_oth,x_oth,xb_up,xb_lo,lambda,fuds,fuds_r),LL,RR,option);

            s_sph1=coll_sph(1:NBH,:)    +   R.*nort1;
            s_sph2=coll_sph(NBH+1:end,:)+   R.*nort2;  
            
            s_up=up_coll+R*[0,1];
            s_lo=lo_coll(2:end-1,:)+ R*[0,-1];
            s_l=l_coll(1:end-1,:)+ R*[-1,0];
            s_r=r_coll(1:end-1,:)+ R*[1,0];
            a_llo=[0.1,0.1]+R*[-sqrt(2)/2,-sqrt(2)/2];
            a_lup=[0.1,1.6]+R*[-sqrt(2)/2,sqrt(2)/2];
            a_rup=[1.1,1.6]+R*[sqrt(2)/2,sqrt(2)/2];
            a_rlo=[1.1,0.1]+R*[sqrt(2)/2,-sqrt(2)/2];
            source=[s_sph1;s_sph2;s_up;s_l;s_lo;s_r;a_llo;a_lup;a_rup;a_rlo];

end

% plot(coll(:,1),coll(:,2),'r.',source(:,1),source(:,2),'b*');   
% hold on;
% xlabel('X');
% ylabel('Y');
% axis equal
% axis ([-0.5 2 -0.5 2])
%% Radial basis function normal
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
BC_1=ones(length(xb_up), 1);
BC_2=0*ones(length(xb_lo), 1);
BC_3=0*ones(length(x_oth), 1);
BC=[BC_1;BC_2;BC_3];
%LHS
DM_cacu_1=DistanceMatrix([xb_up,yb_up],[source(:,1),source(:,2)]);
D_BKM_1=fuds(DM_cacu_1,lambda);
DM_cacu_2=DistanceMatrix([xb_lo,yb_lo],[source(:,1),source(:,2)]);
D_BKM_2=fuds(DM_cacu_2,lambda);
D_BKM=[D_BKM_1;D_BKM_2];
%coll
DM_cacu_3=DistanceMatrix([x_oth,y_oth],[source(:,1),source(:,2)]);
N_BKM=fuds_r(DM_cacu_3,lambda);

N_BKM(isnan(N_BKM))=0;
N_BKM=N_BKM.*nor;
BKM=[D_BKM;N_BKM];
coef=lsqminnorm(BKM,BC);
%% error
%%
%{
Nt=2000;    
load('mesh.mat');
t_tst = randperm(size(data, 1), Nt);
tst = data(t_tst, :);
tst_point=tst(:,1:2);

exact_tst=tst(:,3);
DM_t=DistanceMatrix(tst_point,source);
approx_tst=fuds(DM_t,lambda)*coef;        
error=norm(exact_tst-approx_tst,2)/sqrt(Nt);  
%}

% {
load('mesh_uy.mat'); %%only at A and B
UExact = scatteredInterpolant(data(:,1),data(:,2),data(:,3),'linear');
A=[0.3,1.2];B=[0.5,0.5];
test=[A;B];
exact_tst = UExact(test(:,1),test(:,2));
DM_cacu_t=DistanceMatrix(test,[source(:,1),source(:,2)]);
t1=repmat(test(:,2),1,ns);t2=repmat(N',2,1);
t3=t1-t2;
a=fuds_r(DM_cacu_t,lambda).*t3;
approx_tst=a*coef;
Nt=length(test);
error = abs(exact_tst-approx_tst);
fprintf('error: ');
for i = 1:length(error)
    fprintf('error %d: %.4e\n', i, error(i));
end
%}
%% err save a matrix n*n*n
% err(radius_u0+1,1)=R;
% err(radius_u0+1,2)=error;
%% fprintf 
fprintf('NB1= %3d,NBH= %3d,NB= %3d,Nt=%3d,',NB1,NBH,NB,Nt);
fprintf('R= %4.3f,error = %8.4e\n',R,error); 
% end
toc