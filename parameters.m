% 模型
param.rcell = 8.55e-6; % 细胞半径 [m]
param.Box = [-500,-500;-500,500;500,-500;500,500]; %仿真区域 [μm]
param.dm = 5e-9; % 膜厚度 [m]
% 细胞电气参数
param.sigma_e = 1.58; % 胞外电导 [S/m]
param.sigma_i = 0.3; % 胞内电导 [S/m]
param.sigma_m = 9.5e-9; %膜电导 [S/m]

param.epsilon_e = 6.38e-10; % 外介电常数 [F/m]
param.epsilon_i = param.epsilon_e; % 内介电常数 [F/m]
param.epsilon_m = 4.43e-11; % 膜介电常数 [F/m]

param.Urest = -50e-3; % 膜静息电位 (Um = Vi - Ve) [V]

T = 295.15;    %绝对温度 [K]
k = 1.38065e-23;%玻尔兹曼常数 

% 分子传输参数
%路西法黄 分子1
param.r1 = 0.61e-9;    % Radius of solute 1
param.l1 = 1.46e-9;    % Length of solute 1
param.z1 = -2;         % Valence of solute 1
param.D1_e = 4.77e-10; % 胞外溶质扩散系数 [m^2/s]
param.D1_i = 0.25*4.77e-10; % 胞内溶质扩散系数[m^2/s]
param.c1_i0 = 0;       % 胞内初始分子浓度 [mol/m^3]
param.c1_e0 = 1.0;     % 胞外初始分子浓度 [mol/m^3]

% 孔隙动力学参数
param.rstar = 0.65e-9; %疏水孔向亲水孔转变的孔半径 [m]
param.rmin = 1.0e-9;  %在Um = 0 时，孔径空间的局部极小值
param.max = 12.0e-9; % 最大孔径
param.Wstar_kT = 45; % 孔形成能量壁垒 Pore creation energy barrier [kT]
param.a = 1e9;       % 孔形成速率密度 Pore creation rate density [1/(m^2*s)]
param.alpha = 11;    % 不对称孔生成常数 Asymmetric pore creation constant [kT/V] 
param.beta = 18;     % 对称孔生成常数 Symmetric pore creation constant [kT/V^2]
param.Dp = 2e-13;    % 孔隙扩散系数 Diffusion coefficient in the pore radius space [m^2/s]
param.gamma = 2e-11; % 孔壁张力 Edge tension [N]
param.Gamma0 = 1e-5; % 膜表面张力Surface tension of nonelectroporated membrane [N/m]
param.Gamma1 = 2e-2; % 烃-水界面张力 Hydrocarbon-water interfacial tension [N/m]
param.Fmax = 6.9e-10;% 电力参数中有扩孔倾向的参数Parameter in the electrical force tending to expand the pore [N/V^2]
param.rh = 0.95e-9;  % 电力参数中有扩孔倾向的参数Parameter in the electrical force tending to expand the pore [m]
param.rt = 0.23e-9;  % 电力参数中有扩孔倾向的参数Parameter in the electrical force tending to expand the pore [m]
param.nu = 0.25;     % 孔隙相对入口长度 Relative length of the pore entrance
param.fprot = 0.5;   % 膜蛋白质占比 Areal fraction of proteins
param.taup = 4;      % 孔重新密封时间常数 Pore resealing time constant [s]
