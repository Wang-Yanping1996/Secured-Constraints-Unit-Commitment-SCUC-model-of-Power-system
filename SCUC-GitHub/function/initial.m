mpc = feval(casename);

baseMVA = mpc.baseMVA;
bus = mpc.bus;
gen = mpc.gen;
branch = mpc.branch;
gencost = mpc.gencost;
%%
%一些常数
%Bus type
PQ=1; PV=2; REF=3; NONE=4; 
%Bus
BUS_I=1; BUS_TYPE=2; BUS_PD=3; BUS_QD=4; BUS_GS=5; BUS_BS=6; 
BUS_AREA=7; BUS_VM=8; BUS_VA=9; BUS_baseKV=10; BUS_zone=11; BUS_Vmax=12; BUS_Vmin=13;
%Gen
GEN_BUS=1; GEN_PG=2; GEN_QG=3; GEN_QMAX=4; GEN_QMIN=5; GEN_VG=6; GEN_MBASE=7; GEN_STATUS=8; 
GEN_PMAX=9; GEN_PMIN=10; GEN_Pc1=11; GEN_Pc2=12; GEN_Qc1min=13; GEN_Qc1max=14; GEN_Qc2min=15; GEN_Qc2max=16; 
GEN_ramp_agc=17; GEN_ramp_10=18; GEN_ramp_30=19; GEN_ramp_q=20; GEN_apf=21;
%Branch
F_BUS=1; T_BUS=2; BR_R=3; BR_X=4; BR_B=5; RATE_A=6; RATE_B=7; RATE_C=8;% standard notation (in input)
BR_RATIO=9; BR_ANGLE=10; BR_STATUS=11; BR_angmin=12; BR_angmax=13;% standard notation (in input)
BR_COEFF = 14; BR_MINDEX = 15;
%Gencost
GENCOST_TYPE=1; GENCOST_START=2; GENCOST_DOWN=3; GENCOST_N=4; GENCOST_C2=5; GENCOST_C1=6; GENCOST_C0=7;
%%
% --- convert bus numbering to internal bus numbering
i2e	= bus(:, BUS_I);
e2i = zeros(max(i2e), 1);
e2i(i2e) = [1:size(bus, 1)]';
bus(:, BUS_I)	= e2i( bus(:, BUS_I)	);
gen(:, GEN_BUS)	= e2i( gen(:, GEN_BUS)	);
branch(:, F_BUS)= e2i( branch(:, F_BUS)	);
branch(:, T_BUS)= e2i( branch(:, T_BUS)	);
branch_f_bus = branch(:, F_BUS);
branch_t_bus = branch(:, T_BUS);

%%
%一些用到的数组长度
n_gen = size(gen, 1);
n_bus = size(bus, 1);
n_branch = size(branch, 1);

%%
%机组组合需要引进的数据
GEN_UT=1; GEN_SU=2; GEN_SD=3; GEN_RU=4; GEN_RD=5; GEN_TCOLD=6;
MIN_UP=7; MIN_DOWN=8; INIT_UP=9; INIT_DOWN=10; COST_MAX=11; TIME_DELAY=12;
RU = mpc.SCUC_data(:, GEN_RU)/baseMVA;               %ramp-up 机组爬升限制
SU = mpc.SCUC_data(:, GEN_SU)/baseMVA;               %startup 机组开机限制
RD = mpc.SCUC_data(:, GEN_RD)/baseMVA;               %ramp-down 机组下坡限制
SD = mpc.SCUC_data(:, GEN_SD)/baseMVA;               %shutdown 机组关机限制

min_up = mpc.SCUC_data(:, MIN_UP);           %最小开机时间
min_down = mpc.SCUC_data(:, MIN_DOWN);       %最小停机时间
init_up = mpc.SCUC_data(:, INIT_UP);         %开始调度前开机时间
init_down = mpc.SCUC_data(:, INIT_DOWN);     %开始调度前停机时间
%开始调度前每个机组状态
init_state = zeros(n_bus, 1);
%开始调度前每个机组有功出力
init_gen_P = zeros(n_bus, 1);
for i = 1 : n_gen
    if (init_up(i) > 0 && init_down(i) == 0)    %初始时机组在运行
        init_state(gen(i,GEN_BUS), 1) = 1;
        init_gen_P(gen(i,GEN_BUS), 1) = (gen(i, GEN_PMAX) + 0*gen(i, GEN_PMIN))/2/baseMVA;
    elseif (init_up(i) == 0 && init_down(i) > 0)    %机组不在运行
        init_state(gen(i,GEN_BUS), 1) = 0;
    else
        error('开始调度前运行时间和停机时间必有一个为0，另一个为正');
    end
end
%机组开机费用 Cmax*(1-exp(-t*TIMEDELAY))
start_cost = (mpc.SCUC_data(:, COST_MAX)*ones(1,n_T+max(init_down))).*(1-exp(-mpc.SCUC_data(:,TIME_DELAY)*[1: n_T+max(init_down)]));
%%
%机组发电曲线 二次函数线性化
%均匀分成n_L段
P_interval = zeros(n_gen, n_L+1);
for i = 1: n_gen
    P_interval(i, :) = gen(i, GEN_PMIN): (gen(i, GEN_PMAX)-gen(i, GEN_PMIN))/n_L: gen(i, GEN_PMAX);
end
%出力最小值
A_gen = gencost(:, GENCOST_C2).*gen(:,GEN_PMIN).^2 + gencost(:, GENCOST_C1).*gen(:,GEN_PMIN) + gencost(:, GENCOST_C0);
%各段斜率
Fij = zeros(n_gen, n_L);
for i = 1: n_gen
    for l = 1: n_L
        Fij(i, l) = ((gencost(i, GENCOST_C2).*P_interval(i,l+1).^2 + gencost(i, GENCOST_C1).*P_interval(i,l+1) + gencost(i, GENCOST_C0)) - ...
                     (gencost(i, GENCOST_C2).*P_interval(i,l).^2 + gencost(i, GENCOST_C1).*P_interval(i,l) + gencost(i, GENCOST_C0)))/(P_interval(i,l+1)-P_interval(i,l));
    end
%     Fij(i, :) = 2*gencost(i,GENCOST_C2)*P_interval(i, 1: n_L) + gencost(i, GENCOST_C1);
end