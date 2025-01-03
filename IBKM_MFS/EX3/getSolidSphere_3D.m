function epoints2=getSolidSphere_3D(ns,radius_s,coll,center)
% GETSOLIDSPHERE_3D GETSOLIDSPHERE3 m1=2;
[bdpt_1,~]=generateB_LEI(ns,1,'sphere');
bdpt=radius_s*bdpt_1';
neval=50;
bmin1=min(bdpt,[],1); bmax1=max(bdpt,[],1); %bounding box
xgrid=linspace(bmin1(1),bmax1(1),neval);
ygrid=linspace(bmin1(2),bmax1(2),neval);
zgrid=linspace(bmin1(3),bmax1(3),neval);
[xe1,ye1,ze1]=meshgrid(xgrid,ygrid,zgrid);       % Meshgrid points
epoints0=[xe1(:),ye1(:),ze1(:)]; 
save('epoints0.mat','epoints0')
load('epoints0.mat')
% scatter3(epoints(:,1),epoints(:,2),epoints(:,3))
% {
% % %inner nodes of this sphere
radi=epoints0(:,1).^2+epoints0(:,2).^2+epoints0(:,3).^2;
[idx_outsphere,~]=find(radi>radius_s.^ 2);
epoints0(idx_outsphere,:)=[];
% % make sure the dsites and the ghost points have same size
%  coll=coll(1:ns,:);  %%%%这一行的作用是令coll与ghost点数量不一样，ghost只取ns个
tmp1=floor(size(epoints0,1)./size(coll,1));
tmp2=size(epoints0,1)-tmp1*size(coll,1);
if tmp2==0
    epoints1=epoints0(1:tmp1:end,:);
else
idx_coll=floor(size(epoints0,1)./tmp2);
epoints0(1:idx_coll:idx_coll*tmp2,:)=[];
epoints1=epoints0(1:tmp1:end,:);
end
epoints2=epoints1+center;
end