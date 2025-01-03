function [xs, ys,zs] = pss3D(m, ns)
    % 生成立方体内的均匀分布点
    p = haltonset(3, "Leap", 1000);
    q = net(p, m*25);
    qq = q*6-3; % 将点映射到[-1, 1]立方体内
    % 计算每个点到原点的距离
    dist = sqrt(sum(qq.^2, 2));       
    % 筛选出距离大于R且小于或等于2R的点，即位于半径为2的球与单位球之间的圆环中的点
%   R=1;
    inRing = (dist >1) & (dist <= 2);
    ringPoints = qq(inRing, :);  
    % 提取圆环中的点的坐标
    xs = ringPoints(1:ns,1);
    ys = ringPoints(1:ns,2);
    zs = ringPoints(1:ns,3); 
    
% %     innRing = (dist <1);
% %     rirclePoints = qq(innRing, :); 
% %     xss = rirclePoints(1:250,1);
% %     yss = rirclePoints(1:250,2);
% %     zss = rirclePoints(1:250,3);   
% %     
% %     iRing = (dist >2) & (dist <= 3);
% %     Points = qq(iRing, :); 
% %     xsss = Points(1:250,1);
% %     ysss = Points(1:250,2);
% %     zsss = Points(1:250,3);   
% scatter3(xs,ys,zs,10,'filled', 'MarkerFaceColor',[0.4660,0.6740,0.1880])
% hold on;
% % scatter3(xss,yss,zss,10,'filled', 'MarkerFaceColor',[0.9290,0.6940,0.1250])
% % hold on;
% % scatter3(xsss,ysss,zsss,10,'filled', 'MarkerFaceColor','b')
% axis equal; axis tight;  
% xlim([-3, 3]);
% ylim([-3, 3]);
% zlim([-3, 3]); 
end     