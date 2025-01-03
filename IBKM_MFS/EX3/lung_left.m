function [coll,tspt,nor]=lung_left(nb)
load('tvll.mat','lung_left')  %%lung_left测试点设为NT1=300
%%
m=lung_left(1:2:end,:);
p=lung_left(2:2:end,:);
normals=m;
dsites=p;
idx_test=floor(1:14:size(dsites,1)-58); %%%决定测试点的数量和所在行数（由不同文件数据数量导致不同）
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