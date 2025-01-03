function [collc,colls,colls1,colls2,colls3,coll,tspt]=IBKM_tooth_hollow(NB1,NBH,NTH)
%%inn sphere boundary
[colls,~]=generateB_LEI(NBH,1,'sphere');
rh1 = 0.3; rh2 = 0.2; rh3 = 0.15;% radius of hollow sphere
colls1= rh1 *colls'+[0,0,0.4]; 
colls2= rh2 *colls'+[0.3,0.4,0.1];%%50ghost
% [colls3,~]=generateB_LEI(70,1,'sphere');
colls3= rh3 *colls'+[-0.4,-0.4,0];%%40ghost
colls=[colls1;colls2;colls3];
%%%test
[tspts,~]=generateB_LEI(NTH,1,'sphere');
tspts1= rh1 *tspts'+[0,0,0.4];
tspts2= rh2 *tspts'+[0.3,0.4,0.1];
tspts3= rh3 *tspts'+[-0.4,-0.4,0]; 
tspts=[tspts1;tspts2;tspts3];
% plot3(xbh,ybh,zbh,'r.'); 
% plot3(xth,yth,zth,'y.')
% hold on
%% outside tooth boundary
[collc,tsptc]=tooth_surface(NB1); %%%NT1=1000;通过子代码设置
coll=[collc;colls];   tspt=[tsptc;tspts];
%%
% plot3(colls(:,1),colls(:,2),colls(:,3),'r*'); %%,tspt(:,1),tspt(:,2),tspt(:,3),'b*'
% hold on;
% plot3(collc(:,1),collc(:,2),collc(:,3),'b.');  %%tsptc(:,1),tsptc(:,2),tsptc(:,3),'b*'); 
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% hold on;
%  figure ('color',[1 1 1]);
end
%%
function [coll,tspt,nor]=tooth_surface(nb)
load('tooth_surface.mat','surface')  %%thorax测试点设为NT1=1000 %%%coll、tspt、nor、center
%% generate normal vectors by delaunay
T = delaunay(surface(:,1),surface(:,2));   %%%对于边缘光滑、图形简单的问题域可以用这个三角剖分
TR = triangulation(T,surface(:,1),surface(:,2),surface(:,3)); 
normals= vertexNormal(TR);

dsites=surface;
idx_test=floor(1:3.2:size(dsites,1)-17); %%%决定测试点的数量和所在行数（由不同文件数据数量导致不同）
tspt=dsites(idx_test,:);
dsites(idx_test,:)=[];      %%%%%去除数据中的测试点，使配置点与测试点不同
normals(idx_test,:)=[];
% % make sure the dsites and nb have same size
tmp_nb1=floor(size(dsites,1)./nb);
tmp_nb2=size(dsites,1)-tmp_nb1*nb;
if tmp_nb2==0
    coll=dsites(1:tmp_nb1:end,:);
    nor=normals(1:tmp_nb1:end,:);
else
idx_delate=floor(size(dsites,1)./tmp_nb2);
dsites(1:idx_delate:idx_delate*tmp_nb2,:)=[];
normals(1:idx_delate:idx_delate*tmp_nb2,:)=[];

coll=dsites(1:tmp_nb1:end,:);
nor=normals(1:tmp_nb1:end,:);
end
bmin=min(dsites,[],1);bmax=max(dsites,[],1);
center=(bmin+bmax)/2;
end
