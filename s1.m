%clear;clc;close all;
Box = [-500,-500;-500,500;500,-500;500,500]; %仿真范围
theata = 0:2*pi/600:2*pi;
R = 8.55; %%细胞半径
d_m = 0.005; %%膜厚度
%内膜
cm_n(:,1) = R * cos(theata);%%x
cm_n(:,2) = R * sin(theata);%%y
%外膜
cm_w(:,1) = (R + d_m) .* cos(theata);
cm_w(:,2) = (R + d_m) .* sin(theata);
%固定节点
pfix = [Box;cm_n;cm_w];
%%三角剖分
figure(1)
hold on;
fd = @(p) drectangle(p,-500,500,-500,500);
fh = @(p) (0.3 + 0.3*abs( dcircle(p,0,0,8.55) ));
[p,t] =  distmesh2d(fd,fh,0.2,[-500,-500;500,500],pfix);

%%配色
patch('vertices', p, 'faces', t, 'facecolor', [.9, .9, .9] );



