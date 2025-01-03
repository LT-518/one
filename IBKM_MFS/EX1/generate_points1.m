function [up_coll,lo_coll,lr_coll,l_coll,r_coll,coll_sph,coll] = generate_points1(NB1,NBH)
%% outside boundary coll
n1=1*NB1/5+1;   %%squ
n2=3*NB1/10+1;   
x=linspace(0.1,1.1,n1);
y=linspace(0.1,1.6,n2);
xb_up=x(2:n1-1)';  yb_up=1.6*ones(1,length(x)-2)';  up_coll=[xb_up,yb_up];
xb_lo=x';  yb_lo= 0.1*ones(1,length(x))';           lo_coll=[xb_lo,yb_lo];
xb_l=0.1*ones(1,length(y)-1)';  xb_r=1.1*ones(1,length(y)-1)';
yb_l=y(2:n2)'; yb_r= y(2:n2)';  
l_coll=[xb_l,yb_l]; r_coll=[xb_r,yb_r];
lr_coll=[l_coll;r_coll];
%%
[xb2_0,yb2_0]=pol2cart((1:NBH)/NBH*2*pi,0.2);  % circle 

xb2_l = xb2_0+ 0.5; yb2_l = yb2_0+ 1.2;
xb2_r = xb2_0+ 0.7; yb2_r = yb2_0+ 0.5;
xb2_l=xb2_l';
xb2_r=xb2_r';
xb2=[xb2_l; xb2_r];

yb2_l=yb2_l';
yb2_r=yb2_r';
yb2=[yb2_l; yb2_r];
coll_sph=[xb2,yb2];
coll=[up_coll;r_coll;lo_coll;l_coll;coll_sph];

end