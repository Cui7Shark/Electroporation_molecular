clc;clear;tic;
global A B C Cap ;
load('p.mat');
load('Boun.mat');
load('cm_suoyin.mat');
load('Cap.mat');
load('Res.mat');
load('Rus.mat');
load('c1.mat');
p=roundn(p,-15);
Inner_n_idex = find(roundn(abs(p(:,2)),-4)~=500);  %内节点索引，排除正负极板
A = size(Inner_n_idex,1); %存储所有节点电位
B = size(cm_suoyin,1);    %存储膜节点的孔密度、孔半径 
C = size(Inner_n_idex,1); %存储所有节点分子浓度
%%
%%计算膜节点角度并排序  
Mem_n_idex=zeros(length(cm_suoyin),4);  %膜节点索引  
for i=1:size(cm_suoyin,1)
    for j = 1:size(cm_suoyin,2)
        Mem_n_idex(i,j) = find(Inner_n_idex == cm_suoyin(i,j));
    end
    Mem_n_idex(i,3) = Boun(cm_suoyin(i,1),cm_suoyin(i,2));
    if abs(atand(p(Inner_n_idex(Mem_n_idex(i,1)),2)/p(Inner_n_idex(Mem_n_idex(i,1)),1))) == 90
        if p(Inner_n_idex(Mem_n_idex(i,1)),2)>0
            Mem_n_idex(i,4) = 90;
        elseif p(Inner_n_idex(Mem_n_idex(i,1)),2)<0
            Mem_n_idex(i,4) = 270;
        end
    elseif p(Inner_n_idex(Mem_n_idex(i,1)),1) > 0 && p(Inner_n_idex(Mem_n_idex(i,1)),2) > 0
        Mem_n_idex(i,4) = atand(p(Inner_n_idex(Mem_n_idex(i,1)),2) / p(Inner_n_idex(Mem_n_idex(i,1)),1));
    elseif p(Inner_n_idex(Mem_n_idex(i,1)),1) < 0 && p(Inner_n_idex(Mem_n_idex(i,1)),2) > 0
        Mem_n_idex(i,4) = atand(p(Inner_n_idex(Mem_n_idex(i,1)),2) / p(Inner_n_idex(Mem_n_idex(i,1)),1)) + 180;
    elseif p(Inner_n_idex(Mem_n_idex(i,1)),1) < 0 && p(Inner_n_idex(Mem_n_idex(i,1)),2) < 0
        Mem_n_idex(i,4) = atand(p(Inner_n_idex(Mem_n_idex(i,1)),2) / p(Inner_n_idex(Mem_n_idex(i,1)),1)) + 180;
    elseif p(Inner_n_idex(Mem_n_idex(i,1)),1) > 0 && p(Inner_n_idex(Mem_n_idex(i,1)),2) < 0
        Mem_n_idex(i,4) = atand(p(Inner_n_idex(Mem_n_idex(i,1)),2) / p(Inner_n_idex(Mem_n_idex(i,1)),1)) + 360;
    end
end
Mem_n_idex=sortrows(Mem_n_idex,4);

%%
%matrix coefficient %%（算一算Cap的秩或者行列式 确定乘的系数的大小）--如果M的行列式为0或inf 后面无法计算微分方程
Matrix_coe = -8.37e13;
M = eye(A+B+C);                %% A-电传输 B-孔传输 C-分子传输
M(1:A,1:A) = Matrix_coe * Cap;       
opts = odeset('Mass',M);
%%%%初始条件
u0(1:A) = 0;
u0(Mem_n_idex(:,1)) = -50e-3;       %%细胞膜内膜电位-50mv
u0(A+1:A+B) = 1.5e9;        %%初始孔隙密度N0
u0(A+B+1:A+B+C) = c1;       %%初始分子浓度c1  %可能错了 维度问题
%%
%%解微分方程
tspan = [0 1e-3]; %%0-1ms
[t,u] = ode15s(@transmem,tspan,u0,opts);
%
disp('完成S5');
toc
 



