%SCUC采用流潮流方程，但是采用了二阶锥松弛
clc;
clear;
close all;
clear all class;
warning off;

addpath('../example');
addpath('../function');

alltime = tic;
tic
%%
% 读入
% casename = input('Please enter case name : ', 's');
casename = 'case14mod_SCUC';
% casename = 'case30mod';
k_safe = 0.95;          %安全系数，用于留一定的裕度，针对潮流安全约束

n_T = 24;   % 时段数t，用于机组组合优化
n_L = 20;   % 发电机曲线 二次函数 分段线性化

% 初始化文件
initial;
PD = bus(:, BUS_PD)/baseMVA;
QD = bus(:, BUS_QD)/baseMVA;
% PD = PD*ones(1, n_T);
% QD = QD*ones(1, n_T);
% 24小时的负荷数据
Q_factor = QD/sum(QD);
P_factor = PD/sum(PD);
%P_sum = sum(PD)-sum(PD)/2*sin(pi/12*[0:n_T-1]+pi/3);
P_sum = mpc.PD'/baseMVA;
QD = Q_factor*sum(QD)*P_sum/sum(PD);
PD = P_factor*P_sum;
spnningReserve = 1.02*P_sum;
% gbus = gen(:, GEN_BUS);             % form generator vector

%%
%导纳矩阵计算
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);   % build admitance matrix
% G = real(Ybus);
% B = imag(Ybus);
% % G B 矩阵中都是导纳值的负数 
% g = -G;
% b = -B;

%%
% 创建决策变量
% 发电机出力 非发电机节点取0
gen_P = sdpvar(n_bus, n_T);
gen_Q = sdpvar(n_bus, n_T);
gen_P_upper = sdpvar(n_bus, n_T);   %发电机有功出力上界
% 各节点电压幅值 相角
% Vm = sdpvar(n_bus, 1);      %幅值
% Va = sdpvar(n_bus, 1);      %相角
% 松弛变量 代替
x_i = sdpvar(n_bus, n_T);                         %xi = Vi^2
% xij_1 = sdpvar(n_bus, n_bus, n_T, 'skew');           %xij1 = Vi*Vj*sin(Oi-Oj)  反对称矩阵 A = -A'
% xij_2 = sdpvar(n_bus, n_bus, n_T, 'symmetric');      %xij2 = Vi*Vj*cos(Oi-Oj)  对称矩阵   A = A'
xij_1 = sdpvar(n_branch, n_T);                    %xij_1 = Vi*Vj*sin(Oi-Oj)
xij_2 = sdpvar(n_branch, n_T);                    %xij_2 = Vi*Vj*cos(Oi-Oj)

% 各支路潮流
PF_D = sdpvar(n_branch, n_T);     %P Flow Direct 正向有功潮流 1->2
QF_D = sdpvar(n_branch, n_T);     %Q Flow Direct 正向无功潮流 1->2
PF_R = sdpvar(n_branch, n_T);     %P Flow Reverse 反向有功潮流 2->1
QF_R = sdpvar(n_branch, n_T);     %Q Flow Reverse 反向无功潮流 2->1

% 机组状态
u_state = binvar(n_bus, n_T);     %按母线数，非发电机节点取0
C = [];     %约束
% C = sdpvar(C)>=0;

assign(xij_1, 0);
assign(xij_2, 1);
assign(x_i, 1);

%% 
% 机组组合约束
%系统功率平衡约束
for t = 1: n_T
    C = [C,
        sum(gen_P(gen(:,GEN_BUS),t)) >= sum(PD(:,t)),
        ];
end
%%
% 旋转备用约束
for t = 1: n_T
    C = [C,
        sum(gen_P_upper(gen(:, GEN_BUS), t)) >= sum(spnningReserve(:, t))
        ];
end
%%
% 机组组合约束
% 爬坡限制
for t = 1: n_T
    if (t > 1)
    C = [C,
        %这个约束是照2006年文献
        %A Computationally Efficient Mixed-Integer Linear Formulation for the Thermal Unit Commitment Problem
        %写的，不知道Pmax*(1-u)这项有什么用
        %爬坡限制和启动限制 (ramp-up & startup)    (18)
        gen_P_upper(gen(:,GEN_BUS),t) <= gen_P(gen(:,GEN_BUS),t-1) + RU.*u_state(gen(:,GEN_BUS),t-1) + ...
                                         SU.*(u_state(gen(:,GEN_BUS),t)-u_state(gen(:,GEN_BUS),t-1)) + ...
                                         (gen(:, GEN_PMAX)/baseMVA).*(1-u_state(gen(:,GEN_BUS),t)),
        %下坡限制 (ramp-down)       (20)
        gen_P(gen(:,GEN_BUS),t-1) - gen_P(gen(:,GEN_BUS),t) <= RD.*u_state(gen(:,GEN_BUS),t) + ...
                                                               SD.*(u_state(gen(:,GEN_BUS),t-1)-u_state(gen(:,GEN_BUS),t)) + ...
                                                               (gen(:, GEN_PMAX)/baseMVA).*(1-u_state(gen(:,GEN_BUS),t-1)),
        ];
    end
    if (t < n_T)
        C = [C,
            %关机限制 (shutdown)    (19)
            gen_P_upper(gen(:,GEN_BUS),t) <= (gen(:, GEN_PMAX)/baseMVA).*u_state(gen(:,GEN_BUS),t+1) + ...
                                             SD.*(u_state(gen(:,GEN_BUS),t)-u_state(gen(:,GEN_BUS),t+1)),
                                             ];
    end
end
%%
% 机组组合约束
%最小启动时间限制
for i = 1: n_gen
    Gi = min(n_T, (min_up(i)-init_up(i))*init_state(gen(i,GEN_BUS)));
    %分三段 起始，中间，结尾
    %起始时考虑初始状态的影响，中间不需考虑太多，结尾保证开了以后没关过
    if (Gi >= 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS), [1: Gi])) == Gi,
            ];
    end
    for t = Gi+1: n_T-min_up(i)+1
        if (t > 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS),[t: t+min_up(i)-1])) >= min_up(i).*(u_state(gen(i,GEN_BUS),t)-u_state(gen(i,GEN_BUS),t-1)),
            ];
        elseif (t == 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS),[t: t+min_up(i)-1])) >= min_up(i).*(u_state(gen(i,GEN_BUS),t)-init_state(gen(i,GEN_BUS))),
            ];        
        else
        end
    end
    for t = n_T-min_up(i)+2: n_T
        if (t > 1)
            C = [C,
                sum(u_state(gen(i,GEN_BUS),[t: n_T])) >= (n_T-t+1).*(u_state(gen(i,GEN_BUS),t)-u_state(gen(i,GEN_BUS),t-1)),
                ];
        elseif (t == 1)
            C = [C,
                sum(u_state(gen(i,GEN_BUS),[t: n_T])) >= (n_T-t+1).*(u_state(gen(i,GEN_BUS),t)-init_state(gen(i,GEN_BUS))),
                ];            
        else
        end
    end
end
%%
% 机组组合约束
%最小关机时间限制 同最小开机时间限制
for i = 1: n_gen
    Li = min(n_T, (min_down(i)-init_down(i))*(1-init_state(gen(i,GEN_BUS))));
    if (Li >= 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS), [1: Li])) == 0,
            ];
    end
    for t = Li+1: n_T-min_down(i)+1
        if (t > 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS), [t: t+min_down(i)-1])) <= min_down(i).*(1-u_state(gen(i,GEN_BUS),t-1)+u_state(gen(i,GEN_BUS),t)),
            ];
        elseif (t == 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS), [t: t+min_down(i)-1])) <= min_down(i).*(1-init_state(gen(i,GEN_BUS))+u_state(gen(i,GEN_BUS),t)),
            ];
        else
        end
    end
    for t = n_T-min_down(i)+2: n_T
        if (t > 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS),[t: n_T])) <= (n_T-t+1).*(1-u_state(gen(i,GEN_BUS),t-1)+u_state(gen(i,GEN_BUS),t)),
            ];
        elseif (t == 1)
        C = [C,
            sum(u_state(gen(i,GEN_BUS),[t: n_T])) <= (n_T-t+1).*(1-init_state(gen(i,GEN_BUS))+u_state(gen(i,GEN_BUS),t)),
            ];
        else
        end            
    end
end

%%
New_Br_temp = 1: n_bus;
New_Br_temp(gen(:, GEN_BUS)) = [];
C = [C,
    gen_P(New_Br_temp, :) == 0,
    gen_Q(New_Br_temp, :) == 0,
    gen_P_upper(New_Br_temp, :) == 0,
    u_state(New_Br_temp, :) == 0      %考虑机组组合 非发电机节点取0
    ];  %非发电机节点有功无功为0
%%
% 关于引入的松弛变量x的约束 
for i = 1: n_branch
    m = branch(i, F_BUS);
    n = branch(i, T_BUS);
    for t = 1: n_T
    C = [C,
%         (2*xij_1(i,t)).^2 + (2*xij_2(i,t)).^2 + (x_i(m,t)-x_i(n,t)).^2 <= (x_i(m,t)+x_i(n,t)).^2
        rcone([xij_1(i,t); xij_2(i,t)], 0.5*x_i(m,t), x_i(n,t)),          
        ];
    end
end
%%
%潮流方程
% 支路潮流约束
for i = 1: n_branch
    f_bus = branch_f_bus(i);            % 支路i的起始母线  
    t_bus = branch_t_bus(i);            % 支路i的终端母线
    
    gff = real(Yf(i,branch_f_bus(i)));
    gft = real(Yf(i,branch_t_bus(i)));
    bff = imag(Yf(i,branch_f_bus(i)));
    bft = imag(Yf(i,branch_t_bus(i)));
    
    gtf = real(Yt(i,branch_f_bus(i)));
    gtt = real(Yt(i,branch_t_bus(i)));
    btf = imag(Yt(i,branch_f_bus(i)));
    btt = imag(Yt(i,branch_t_bus(i)));
    C = [C,
        PF_D(i,:) == x_i(f_bus,:)*gff+xij_2(i,:)*gft+xij_1(i,:)*bft,
        QF_D(i,:) == -x_i(f_bus,:)*bff+xij_1(i,:)*gft-xij_2(i,:)*bft,

        PF_R(i,:) == x_i(t_bus,:)*gtt+xij_2(i,:)*gtf-xij_1(i,:)*btf,
        QF_R(i,:) == -x_i(t_bus,:)*btt-xij_1(i,:)*gtf-xij_2(i,:)*btf
        ];  %各线路潮流方程
end

%%
%节点功率平衡约束
for i = 1: n_bus
    for t = 1: n_T
    C = [C,
        %转换成 n_T = 24 后，sum这么写应该是把行加起来，合成一行
        gen_P(i,t) == PD(i,t) + ...
                    sum(PF_D(branch(:, F_BUS) == i,t)) + ...
                    sum(PF_R(branch(:, T_BUS) == i,t)) + ...
                    x_i(i,t)*bus(i, BUS_GS)/baseMVA,
        gen_Q(i,t) == QD(i,t) + ...
                    sum(QF_D(branch(:, F_BUS) == i,t)) + ...
                    sum(QF_R(branch(:, T_BUS) == i,t)) - ...
                    x_i(i,t)*bus(i, BUS_BS)/baseMVA             %考虑Gs Bs影响
        ];      %节点功率平衡方程     %x_i = Vi^2
    end
end

%%
%各节点电压幅值约束
for t = 1: n_T
    C = [C,
        bus(:, BUS_Vmax).^2 >= x_i(bus(:, BUS_I), t) >= bus(:, BUS_Vmin).^2
        ];          %x_i = Vi^2
end
%%
%发电机有功出力约束
for  t = 1: n_T
    for i = 1: n_gen
        C = [C,
            gen_P_upper(gen(i, GEN_BUS),t) >= gen_P(gen(i, GEN_BUS),t) >= u_state(gen(i, GEN_BUS),t).*gen(i, GEN_PMIN)/baseMVA,
            u_state(gen(i, GEN_BUS),t).*gen(i, GEN_PMAX)/baseMVA >= gen_P_upper(gen(i, GEN_BUS),t) >= 0
            ];
    end
end
%%
%发电机无功出力约束
for t = 1: n_T
    for i = 1: n_gen
        C = [C,
            %考虑机组组合
            u_state(gen(i, GEN_BUS),t).*gen(i, GEN_QMAX)/baseMVA >= gen_Q(gen(i, GEN_BUS),t) >= u_state(gen(i, GEN_BUS),t).*gen(i, GEN_QMIN)/baseMVA
            ];
    end
end

%%
%支路潮流约束
% -Pmax <=P <= Pmax
for i = 1: n_branch
    if (branch(i, RATE_A) ~= 0)     %rateA为0则认为不需要添加安全约束
        C = [C,
            -k_safe*branch(i, RATE_A)/baseMVA <= PF_D(i,:) <= k_safe*branch(i, RATE_A)/baseMVA
            ];
    end
end

%%
%故障态约束

%%
%发电机成本函数，仅考虑2次函数  如果多次要重写
obj_value = sum(gencost(:, GENCOST_C2)'*(gen_P(gen(:, GEN_BUS),:)*baseMVA).^2) + ...
            sum(gencost(:, GENCOST_C1)'* gen_P(gen(:, GEN_BUS),:)*baseMVA) + ...
            sum(gencost(:, GENCOST_C0)'*u_state(gen(:, GEN_BUS),:));
        
%%     
ops = sdpsettings('solver','cplex','verbose',2,'usex0',1);       %使用初值 主要是电压取1
% ops.mosek.MSK_IPAR_MIO_CONIC_OUTER_APPROXIMATION = 'MSK_ON';
% ops.moesk.MSK_DPAR_OPTIMIZER_MAX_TIME = 300;
% ops.mosek.MSK_IPAR_MIO_HEURISTIC_LEVEL = 4;
% ops.mosek.MSK_DPAR_MIO_TOL_REL_GAP = 1e-6;
% ops.mosek.MSK_DPAR_MIO_REL_GAP_CONST = 1e-4;
%mosek 不需要yalmip转换二阶锥约束，cplex和gurobi需要 很耗内存
% ops = sdpsettings('verbose',2);
% ops.sedumi.eps = 1e-6;
% ops.gurobi.MIPGap = 1e-6;
%%
toc
%求解         
result = optimize(C, obj_value, ops);
toc(alltime)
if result.problem == 0 % problem =0 代表求解成功
%     value(x)
%     value(opf_value)   
else
%     tic
%     ops = sdpsettings('solver','cplex','verbose',2,'usex0',1);       %使用初值 主要是电压取1
%     result = optimize(C, opf_value, ops);
%     toc
    error('求解出错');
end  
plot([0: n_T], [init_gen_P(gen(:,GEN_BUS)) value(gen_P(gen(:,GEN_BUS),:))]);    %各机组出力

gen_P = value(gen_P(gen(:,GEN_BUS),:));
gen_Q = value(gen_Q(gen(:,GEN_BUS),:));
xij_1 = value(xij_1);
xij_2 = value(xij_2);
x_i = value(x_i);
PF_D = value(PF_D);
% QF_D = value(QF_D);
PF_R = value(PF_R);
% QF_R = value(QF_R);
u_state = value(u_state(gen(:,GEN_BUS),:));

