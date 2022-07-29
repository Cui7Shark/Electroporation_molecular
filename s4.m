clc;clear;close all;tic;
load('p.mat');
%%Adj.mat存储邻接矩阵，Boun.mat存储点间距离
load('Adj.mat');
load('Boun.mat');
load('cm_n.mat');
load('cm_w.mat');
p=roundn(p,-100);
cm_n=roundn(cm_n,-100);
cm_w=roundn(cm_w,-100);

%===============================分配传输特性========================================%
VC_th = 8.55;  %VC深度(微米)
R = 8.55;      %细胞半径(微米)
d_M = 0.005;   %膜厚度(微米)
%%  电导率(s/m)
sigma_i = 0.3;
sigma_m = 9.5e-9;
sigma_e = 1.58;
sigma_p = 2*sigma_i * sigma_e / (sigma_i + sigma_e);
%%  介电常数(F/m)
epsilon_i = 6.38e-10;
epsilon_m = 4.43e-11;
epsilon_e = 6.38e-10;
%%  溶质扩散系数(m^2/s)
Ds_e = 47.7e-11;
Ds_i = 11.9e-11;
Ds_p = 2*Ds_i * Ds_e / (Ds_i + Ds_e);
%%  初始分子浓度
c1_i0 = 0;       % 胞内初始分子浓度 [mol/m^3]
c1_e0 = 1.0;     % 胞外初始分子浓度 [mol/m^3]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cm 在膜上的P点的索引列向量
cm_nei_suoyin = zeros(length(cm_n),1);
for i = 1:length(p)
   for j=1:length(cm_n)
    if p(i,:) == cm_n(j,:)
        cm_nei_suoyin(i,1) = i;
    end
   end
end
cm_nei_suoyin(cm_nei_suoyin == 0) =[];

cm_wai_suoyin = zeros(length(cm_w),1);
for i = 1:length(p)
   for j=1:length(cm_w)
    if p(i,:) == cm_w(j,:)
        cm_wai_suoyin(i,1) = i;
    end
   end
end
cm_wai_suoyin(cm_wai_suoyin == 0) =[];
%%计算膜上点的索引 --- 便于下面给区域分配介电参数和电导率
%%%CM
cm_suoyin =zeros(length(cm_nei_suoyin),2);
for i=1:size(cm_nei_suoyin,1)
    for j=1:size(cm_wai_suoyin,1)
        D=sqrt((p(cm_nei_suoyin(i,1),1)-p(cm_wai_suoyin(j,1),1))^2+(p(cm_nei_suoyin(i,1),2)-p(cm_wai_suoyin(j,1),2))^2);%% 求D膜内外两点间的距离即等于膜的厚度5nm
        if roundn(D,-4) == d_M
            cm_suoyin(i,1)=cm_nei_suoyin(i);
            cm_suoyin(i,2)=cm_wai_suoyin(j);%Mem n*2的数组 第一列为内膜点索引，第二列为外膜点索引 //(i,j)是一对在膜上的节点
            break
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%分配电导率、介电常数、溶质扩散系数%%%%%%%%%%%%%%%%%%%%%
sigma=zeros(size(p,1),1);
epsilon=zeros(size(p,1),1);
Ds = zeros(size(p,1),1);
c1 = zeros(size(p,1),1);  %节点的分子初始浓度
for i=1:size(p,1)
    if (size(find(cm_suoyin == i),1) == 1)  %%节点在细胞膜上
        sigma(i,1) = sigma_m;
        epsilon(i,1) = epsilon_m;
        Ds(i,1) = Ds_p;
        c1(i,1) = c1_i0;
    elseif sqrt( (p(i,1) - 0 )^2 +  (p(i,2) - 0)^2 ) > (R + d_M) %%节点在细胞外
        sigma(i,1) = sigma_e;
        epsilon(i,1) = epsilon_e;
        Ds(i,1) = Ds_e;
        c1(i,1) = c1_e0; %胞外分子浓度
    elseif sqrt( (p(i,1) - 0 )^2 +  (p(i,2) - 0)^2 ) < R %%节点在细胞内
        sigma(i,1) = sigma_i;
        epsilon(i,1) = epsilon_i;   
        Ds(i,1) = Ds_i;
        c1(i,1) = c1_i0; %胞内分子浓度
    end
end

%%%%%%%%%%%%%为VC单元分配电阻值电容值%%%%%%%%%%%%%%%%%%%%%

Res_VC=zeros(size(Boun));Cap_VC=zeros(size(Boun));
for i = 1:size(Adj,1)
    for j = 1:size(Adj{i},2)
        sigma_av = (sigma(i) + sigma(Adj{i}(j)))/2;
        epsilon_av = (epsilon(i) + epsilon(Adj{i}(j)))/2;
        l_ij=sqrt((p(i,1)-p(Adj{i}(j),1))^2+(p(i,2)-p(Adj{i}(j),2))^2);
        Cap_VC(i,Adj{i}(j)) =(epsilon_av * (Boun(i,Adj{i}(j)) * 1e-6) * (VC_th * 1e-6))/(l_ij * 1e-6);
        Res_VC(i,Adj{i}(j)) = (sigma_av * (Boun(i,Adj{i}(j)) * 1e-6) * (VC_th * 1e-6))/(l_ij * 1e-6) ;

    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%节点导纳矩阵%%%%%%%%%%%%%%%%%%%%%%%
k1=find(roundn(p(:,2),-4)==500);
k2=find(roundn(abs(p(:,2)),-4)~=500); %%K1是上下边界 K2是内部
Res=zeros(size(k2,1)); 
Cap=zeros(size(k2,1));
Rus=zeros(size(k2,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%《推导正负号》%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(k2,1)
    for j=1:size(k2,1)
        if i==j    %%i=j 主对角线上，求自导纳
            Res(i,j) = sum(Res_VC(k2(i),:));
            Cap(i,j) = -sum(Cap_VC(k2(i),:)); %% https://ww2.mathworks.cn/help/matlab/math/solve-stiff-transistor-dae.html?searchHighlight=ode23t&s_tid=srchtitle_ode23t_2
        else       %%求互导纳 加-号
            Res(i,j) = -Res_VC(k2(i),k2(j));
            Cap(i,j) = Cap_VC(k2(i),k2(j));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%电源导纳矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(k2,1)           %判断Res1中第k2行第k1列元素是否为0,
    for j=1:size(k1,1)
        if Res_VC(k2(i),k1(j))~=0
            Rus(i)=Rus(i) - Res_VC(k2(i),k1(j));  %%排除膜上临近节点电流的影响
        end
    end
end
%%
%save(Cap.m);save(Res.m);save(Rus.m);save(cm_suoyin.m);
% save(Ds.m);save(sigma.m);save(epsilon.m);save(c1.m);
disp('完成！step_4');
toc