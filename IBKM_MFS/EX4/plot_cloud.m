m=2000; ns=1/30*m;
t=(2*pi/m)*(0:m-1)';
% rt=0.275*sqrt(1+3*(cos(t).^2));
rt=(1+(cos(3/2*t).^2))/2;
[x,y]=pol2cart(t,rt);
dx=gradient(x);             dy=gradient(y);
xn=dy./sqrt(dx.^2+dy.^2);    yn=-dx./sqrt(dx.^2+dy.^2); 
%%%%%%%%%%
xr=cos(t);  yr=sin(t);
p=haltonset(2,"Leap",1000); q=net(p,m*25);
qq=q*2-1; %%creat a square domian points
in1=inpolygon(qq(:,1),qq(:,2),xr,yr);
qq=qq(in1,:);  %%choice circle-domian points
in2=inpolygon(qq(:,1),qq(:,2),xr*0.84,yr*0.84);
qq=qq(~in2,:);
xs=qq(1:ns,1); ys=qq(1:ns,2);
hh= plot(xs,ys,'bo','markersize',4,'markerfacecolor',[0 0 1]); 
 set(hh, 'Color', [59/255, 66/255, 147/255]);
axis equal; axis tight;
hold on;
h1= plot(xr*0.81,yr*0.81,'m-','LineWidth',1.5); 
set(h1, 'Color', [253/255, 141/255, 60/255]);
hold on;
h2= plot(xr*1.05,yr*1.05,'m-','LineWidth',1.5); 
set(h2, 'Color', [253/255, 141/255, 60/255]);
xlim([-1.5, 1.5]);
ylim([-1.5, 1.5]);
hold on;
% 定义长方形的中心点、宽度和高度
center_x = 0; % x坐标中心
center_y = 0; % y坐标中心
width = 0.8;  % 长方形的宽度
height = 0.45; % 长方形的高度

% 使用rectangle函数绘制长方形
% 'Position'参数定义了长方形的位置和大小
% 'EdgeColor'设置为'none'以隐藏边框
% 'LineStyle'设置为'--'以绘制虚线
rectangle('Position', [center_x-width/2, center_y-height/2, width, height], ...
          'EdgeColor',[23/255, 134/255, 66/255], 'LineStyle', '--','LineWidth',1.5);  %[23/255, 134/255, 66/255]
text(0.5, 0.5, '\Omega', 'FontSize', 24, 'Color', 'k', 'FontName', 'Symbol');
ax = gca;
set(ax, 'LineWidth', 1.5); 
a=in2;
b=0.5*in2;