function [coll,tspt,nor]=ventricle(nb)
load('tvll.mat','ventricle') 
%%
m=ventricle(1:2:end,:);
p=ventricle(2:2:end,:);
normals=m;
dsites=p;
idx_test=floor(1:2.67:size(dsites,1)-2); 
tspt=dsites(idx_test,:);
dsites(idx_test,:)=[];      
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