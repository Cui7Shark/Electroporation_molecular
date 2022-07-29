function dudt=transmem(t,u)
global A B C Mem_n_idex Res Rus  ;
dudt = zeros(A+B+C,1);
%%
% 模型
par.rcell = 8.55e-6; % 细胞半径 [m]
par.Box = [-500,-500;-500,500;500,-500;500,500]; %仿真区域 [μm]
par.dm = 5e-9; % 膜厚度 [m]
% 细胞电气参数
par.sigma_e = 1.58; % 胞外电导 [S/m]
par.sigma_i = 0.3; % 胞内电导 [S/m]
par.sigma_m = 9.5e-9; %膜电导 [S/m]

par.epsilon_e = 6.38e-10; % 外介电常数 [F/m]
par.epsilon_i = par.epsilon_e; % 内介电常数 [F/m]
par.epsilon_m = 4.43e-11; % 膜介电常数 [F/m]

par.Urest = -50e-3; % 膜静息电位 (Um = Vi - Ve) [V]

T = 295.15;    %绝对温度 [K]
k = 1.38065e-23;%玻尔兹曼常数 
q =2.46;
e=1.60e-19;   %单位电荷量 1.60e-19库伦

% 分子传输参数
%路西法黄 分子1
par.r1 = 0.61e-9;    % Radius of solute 1
par.l1 = 1.46e-9;    % Length of solute 1
par.z1 = -2;         % Valence of solute 1
par.D1_e = 4.77e-10; % 胞外溶质扩散系数 [m^2/s]
par.D1_i = 0.25*4.77e-10; % 胞内溶质扩散系数[m^2/s]
par.c1_i0 = 0;       % 胞内初始分子浓度 [mol/m^3]
par.c1_e0 = 1.0;     % 胞外初始分子浓度 [mol/m^3]

% 孔隙动力学参数
par.rstar = 0.65e-9; %疏水孔向亲水孔转变的孔半径 [m]
par.rmin = 1.0e-9;  %在Um = 0 时，孔径空间的局部极小值
par.max = 12.0e-9; % 最大孔径
par.Wstar_kT = 45; % 孔形成能量壁垒 Pore creation energy barrier [kT]
par.a = 1e9;       % 孔形成速率密度 Pore creation rate density [1/(m^2*s)]
par.alpha = 11;    % 不对称孔生成常数 Asymmetric pore creation constant [kT/V] 
par.beta = 18;     % 对称孔生成常数 Symmetric pore creation constant [kT/V^2]
par.Dp = 2e-13;    % 孔隙扩散系数 Diffusion coefficient in the pore radius space [m^2/s]
par.gamma = 2e-11; % 孔壁张力 Edge tension [N]
par.Gamma0 = 1e-5; % 膜表面张力Surface tension of nonelectroporated membrane [N/m]
par.Gamma1 = 2e-2; % 烃-水界面张力 Hydrocarbon-water interfacial tension [N/m]
par.Fmax = 6.9e-10;% 电力参数中有扩孔倾向的参数Parameter in the electrical force tending to expand the pore [N/V^2]
par.rh = 0.95e-9;  % 电力参数中有扩孔倾向的参数Parameter in the electrical force tending to expand the pore [m]
par.rt = 0.23e-9;  % 电力参数中有扩孔倾向的参数Parameter in the electrical force tending to expand the pore [m]
par.nu = 0.25;     % 孔隙相对入口长度 Relative length of the pore entrance
par.fprot = 0.5;   % 膜蛋白质占比 Areal fraction of proteins
par.taup = 4;      % 孔重新密封时间常数 Pore resealing time constant [s]
%%1ms   1kV/cm
Vapp= @(t) 0.*(t<0)+(1e8 * t).*(t<1e-6 & t>=0) + 100.*(t>=1e-6 & t<999e-6) + (-1e8.*t + 1e5).*(t>=999e-6 & t<1e-3)+ 0.*(t>=1e-3);
t;
%%
for i = 1:(A+B+C)
    [row,col]=find(Mem_n_idex(:,1:2)==i);                                  %%判断点i是否为膜上节点
    if size(col,1)==0 && i<=A                                              %%非膜上节点的电流

        dudt(i) = Res(i,:)*u(1:a) + Rus(i)*Vapp;                           %求电流
    elseif  size(col,1)==1 && i<=A
        Vm=@(u) u(Mem_n_idex(row,2))-u(Mem_n_idex(row,1));                                   %外电势-内电势
        vm=@(u) Vm(u)*X;                                                   %X=q*e/(k_Boltzmann*T);%%95.0725

        Gp=@(u) (sigma*pi*rm*rm/dm) * (exp(vm(u))-1)/((w0*exp(w0-yita*vm(u))-yita*vm(u))*exp(vm(u))/(w0-yita*vm(u))- ...
            ((w0*exp(w0+yita*vm(u))+yita*vm(u))/(w0+yita*vm(u))));
     
%%计算跨膜电流
        %细胞膜内节点
        if col==1 
            dudt(i) = Res(i,:)*u(1:A) - Gp(u)*u(row+a)*Vm(u)*Mem_n_idex(row,3)*Th*1e-12;      %1e-6*1e-6*;
        %细胞膜外-节点
        elseif col==2 
            dudt(i) = Res(i,:)*u(1:a) + Gp(u)*u(row+a)*Vm(u)*Mem_n_idex(row,3)*Th*1e-12;      %1e-6*1e-6;
        end
    end
%%计算孔隙密度
    if  i>A
        Vm1=@(u) (u(Mem_n_idex(i-a,2))-u(Mem_n_idex(i-a,1)));
        dudt(i) = alpha * (exp((Vm1(u)/Vep)^2) * (1-((u(i)/N0)*exp(-q*(Vm1(u)/Vep)^2))));
    end
end
