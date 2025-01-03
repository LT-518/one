function [coll,tspt,nor]=thorax(nb)
load('tvll.mat','thorax') 
%%
m=thorax(1:2:end,:);  %%%nor
p=thorax(2:2:end,:);  %%coll
normals=m;
dsites=p;
idx_test=floor(1:54.1:size(dsites,1)-76); 
tspt=dsites(idx_test,:);
dsites(idx_test,:)=[];   
normals(idx_test,:)=[];

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