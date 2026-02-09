# HSTAR 程序完整执行流程分析

## 一、总体架构

HSTAR 是一个以全局变量驱动的有限元程序，所有状态通过 `Global.f90` 中的 `global_var` 模块共享。
程序的分支极多，由输入参数控制不同的分析路径。

### 核心模块依赖关系

```
Vartype.f90 (ink=int32, irk=float64)
    ↓
Global.f90  (global_var: 所有全局变量、数据结构定义)
    ↓
├── Elements.f90  (形函数、高斯积分、单元库)
├── Material.f90  (本构模型)
├── Stiff.f90     (刚度矩阵计算)
├── Residu.f90    (内力/残差计算)
├── Solver.f90    (PARDISO/PCG/Profile 求解)
├── Prescrib.f90  (边界条件)
├── Load.f90      (外荷载)
├── Temper.f90    (温度分析)
├── Output.f90    (后处理输出)
├── Level.f90     (Level Set 方法)
    ↓
Fem.f90 (主程序 + 所有分析流程子程序, ~18000行)
```

---

## 二、主程序入口 (Fem.f90:4-543)

```
PROGRAM FEM90
│
├─ 读取 inp 文件
│   ├─ text (标题)
│   ├─ restart, relis, sysrelis, ADINA, Uopt_R, gamamax
│   ├─ text
│   └─ probn (问题名称，用于拼接所有输入输出文件名)
│
├─ call global_data                    ← 读取所有输入文件，初始化全局数据
│
├─ [Uopt_R==1] 坐标修正 (modify_coord)
│
├─ 读取 runblks (实际运行的荷载块数)
│
├─ call material_set                   ← 设置材料参数
├─ call modf_element_lib               ← 修改单元库(接触、梁等特殊处理)
├─ call contact_point_to_point         ← 点对点接触初始化
├─ call link_concrete_and_steel        ← 钢筋混凝土连接
├─ call link_concrete_and_water_pipe   ← 冷却水管连接
│
├─ call stiff_interface_fluid_solid    ← 流固界面刚度
├─ call stiff_absorb_fluid             ← 流体吸收边界
├─ call stiff_absorb_solid             ← 固体吸收边界
├─ call stiff_ifs2006                  ← IFS2006 界面
│
├─ call output_read                    ← 读取输出控制参数
│
├─ 分配结果存储数组 (result_zero, result_first, result_second, tofor, stfor 等)
│
├─ [restart/=0] call resta_read_write(1)  ← 从重启文件恢复
│   └─ [restart==2] 仅后处理输出然后 stop
│
├─ 跳过已完成的荷载块 (rewind mainunit, 跳过 lblks 个块)
│
├─── 根据 Bparameter 分支 ────────────────────────
│    │
│    ├─ [-1,-2] 参数反分析验证
│    │   └─ call parameter_back_analysis_verify → stop
│    │
│    ├─ [1,2] & balgor<=1: 最优化反分析
│    │   └─ call parameter_back_analysis
│    │
│    ├─ [1,2] & balgor==2: 信赖域反分析
│    │   └─ call trust_region_back_analysis
│    │
│    ├─ [3] 刚体位移反分析
│    │   └─ call rigid_dis_back_analysis → stop
│    │
│    ├─ [4] 节点值反分析
│    │   └─ call nodal_value_back_analysis → stop
│    │
│    └─ [else] 正常分析
│        └─ call process_analysis          ← ★ 主分析入口
│
├─ [submodel<0] 子模型后处理
│
├─ call out_record_write
│
└─ 结束
```

---

## 三、process_analysis 核心流程 (Fem.f90:1633-1999)

这是正常分析的主入口，包含 **荷载块循环**。

```
subroutine process_analysis
│
├─ 初始化: ttime, 接触变量归零, call gpvar_initial
│
├─ call external_load_1                ← 预读荷载数据
│
└─── 荷载块循环 ═══════════════════════════════════════════
     for iblks = lblks+1 to runblks
     │
     ├─ 更新单元组出现状态: appear(igroup) = appear_process(igroup, iblks)
     │   └─ appear=-1 表示本块被"杀死"的组 (分级施工关键机制)
     │
     ├─ [gamamax/=0] 等效线性化处理
     │
     ├─ [meshc==1,2] 网格细化/粗化: 插值矩阵、结果传递
     │   └─ 涉及 trans, trans_c 数组的大量操作
     │
     ├─ call prescrib_set              ← 设置当前块的边界条件
     │
     ├─ [nbackdT==1] 反分析位移边界处理
     │
     ├─ call external_load_2           ← 读取当前块的荷载数据
     │
     ├─ call contact_pair_process      ← 接触对处理
     ├─ call boundt                    ← 温度边界条件
     │
     ├─ operation='SET'; call solve    ← 求解器初始化(分配存储、符号分解)
     │
     ├─ call modf_time_order           ← 修改时间积分阶次
     │
     ├─── 根据 type_problem 分支 ─────────────────────
     │    │
     │    ├─ 'Q' (准静态):
     │    │   ├─ [relis==1, block_stab==0] call STATIC_U_reli
     │    │   ├─ [relis==1, block_stab>=1] call static_rigid_reli
     │    │   ├─ [block_stab>=1, ebody==0] call static_rigid_1
     │    │   ├─ [mdofn==7, lmdofn(7)/=0] call static_U_P   ← U-P 耦合
     │    │   ├─ [mdofn==8, lmdofn(8)/=0] call static_U_Pw  ← U-Pw 耦合
     │    │   ├─ [nbackf/=0, nbackdT==0]  call back_analysis ← 反分析
     │    │   ├─ [nbackf/=0, nbackdT/=0]  call back_d_analysis
     │    │   └─ [else]                    call static_U      ← ★ 标准静力
     │    │
     │    ├─ 非 'Q','E','W' (动力):
     │    │   ├─ [type_solver=='EXPLICIT'] call explicit      ← 显式积分
     │    │   └─ [else]                    call time_dependent ← ★ Newmark隐式
     │    │
     │    ├─ 'E' (地震反应谱):
     │    │   └─ call response_spectrum
     │    │
     │    └─ 'W' (频率分析):
     │        └─ call frequency_analysis
     │
     ├─ call gpvar1_initial            ← 高斯点变量初始化
     │
     └─ [rmesh/=0] 自适应网格重划分
         └─ call mesh_refine
```

---

## 四、static_U 详细流程 (Fem.f90:3560-4246)

标准准静态分析，三层循环嵌套：**增量→步→迭代**

```
SUBROUTINE STATIC_U
│
├─ 读取 nincs (增量步数)
│
└─── 增量循环 (iincs = lincs+1 to nincs) ═══════════════════
     │
     ├─ 从 mainunit 读取控制参数:
     │   ├─ miter       最大迭代次数
     │   ├─ ditime      时间步长
     │   ├─ noutn       节点结果输出间隔
     │   ├─ noutf       全量输出间隔
     │   ├─ nstep       总步数
     │   ├─ inc_step    步长增量
     │   ├─ nresta      重启写入间隔
     │   ├─ cwater      水荷载标志
     │   └─ Qstatic     静力场标志
     │
     ├─ 读取收敛容差: toler_force, toler_var(1:mdofn)
     │
     ├─ [cwater/=0] 读取水压力系数 coef_water
     ├─ [Qstatic/=0] 读取静力场参数
     │
     └─── 步循环 (istep = inc_step to nstep, step inc_step) ════
          │
          ├─ 更新时间: ttime += ditime
          ├─ call dfact_time_curve(ttime)    ← 荷载时间曲线插值
          ├─ call modf_var_prescribed        ← 更新约束变量
          ├─ call saturation_judge           ← 饱和度判断
          ├─ call gravity                    ← 重力荷载
          ├─ [cwater/=0] call step_water_pressure  ← 水压力
          │
          ├─ call force_external             ← 组装外力向量 tofor
          │   ├─ 体力 (body force)
          │   ├─ 面力 (edge load)
          │   ├─ 集中力 (point load)
          │   ├─ 惯性力 (inertia force)
          │   └─ 温度荷载
          │
          ├─ call load_of_creep_and_temperature  ← 蠕变+温度荷载
          ├─ call creep_strain_of_rock_fill      ← 堆石体蠕变
          ├─ call wetting_strain_of_rock_fill    ← 湿化变形
          │
          └─── 迭代循环 (iiter = 1 to miter) Newton-Raphson ════
               │
               ├── 首次迭代初始化 ──────────────
               │   ├─ call local_stress          ← 计算局部应力
               │   ├─ call contact_state(0)      ← 接触状态初始化
               │   ├─ call predict               ← 位移预测
               │   └─ 材料更新:
               │       ├─ call strain_for_steel_bar
               │       ├─ call stran0_creep4     ← Burgers蠕变
               │       ├─ call effect_stres_modul_for_steel_beam
               │       └─ call stiffness_for_bolt_spring
               │
               ├── 刚度矩阵组装 ──────────────
               │   ├─ call stiff_u               ← 单元刚度矩阵计算
               │   │   └─ 遍历所有活跃组的单元，计算 element(ielem)%estif
               │   │       ├─ 弹性刚度 (dmatrix_* 系列)
               │   │       ├─ 几何刚度 (大变形)
               │   │       └─ 特殊单元 (梁、接触、弹簧)
               │   │
               │   ├─ [kthmat/=0] call htmatrx   ← 热传导矩阵
               │   │
               │   ├─ call estif_assemble         ← 组装到全局矩阵
               │   │   ├─ PARDISO: 组装到 CSR 格式 (sstore, iseq)
               │   │   ├─ PROFILE: 组装到天际线格式
               │   │   └─ JPCG: 保留单元刚度，不组装全局
               │   │
               │   └─ 对角线修正 (防止零主元)
               │
               ├── 残差力计算 ──────────────
               │   ├─ call gpvar2_initial         ← 高斯点变量初始化
               │   │
               │   ├─ [ninit/=0, kinit==2, 首步]  ← 初始应力处理
               │   │   ├─ call eload_initial_stress
               │   │   └─ call force_release
               │   │
               │   ├─ call eload_initialize       ← 单元荷载数组清零
               │   ├─ call residu_f               ← ★ 核心: 内力计算
               │   │   └─ 遍历所有单元:
               │   │       ├─ 提取节点位移 → 计算应变
               │   │       ├─ 调用本构模型 → 计算应力
               │   │       ├─ 积分内力: ∫ B^T σ dV
               │   │       └─ 累加到 element(ie)%eload
               │   │
               │   ├─ call eload_field            ← 场荷载
               │   ├─ call force_internal          ← stfor = Σ eload
               │   │   └─ 将单元荷载组装到全局内力向量
               │   │
               │   └─ [ngaps/=0] call ctfor_to_tofor  ← 接触力→总力
               │
               ├── 刚度矩阵分解 ──────────────
               │   └─ operation='FACTORIZE'; call solve
               │       ├─ PARDISO: 数值分解
               │       ├─ PROFILE: LDL^T 分解
               │       └─ JPCG: 不需要预分解
               │
               ├── 组装残差向量 ──────────────
               │   └─ rvector(ieq) = tofor(itotv) - stfor(itotv)
               │       └─ 处理约束条件 (trans 数组)
               │
               ├── 求解位移增量 ──────────────
               │   ├─ [type_nl==8] call bfgsr     ← 修正Newton (BFGS)
               │   └─ [else] operation='SOLVE'; call solve
               │       └─ 解 K·Δu = rvector
               │
               ├── 接触力求解 ──────────────
               │   └─ [ngaps/=0] call solve_ctt
               │
               ├── 变量更新 ──────────────
               │   ├─ call varupdate              ← 更新 result_zero += Δu
               │   ├─ call eload_initialize
               │   ├─ call residu_f               ← 重新计算内力
               │   ├─ call eload_field
               │   ├─ call reaction_prescribed    ← 支反力
               │   │
               │   ├── 收敛检查 ──────────────
               │   │   ├─ call conver_load        ← 力收敛: |R|/|F| < tol
               │   │   │   └─ 写入 chk 文件: "ratio for residu norm= ..."
               │   │   ├─ call conver_nodal_value ← 位移收敛: |Δu|/|u| < tol
               │   │   └─ nchek==0 → EXIT (收敛!)
               │   │
               │   └─ (不收敛则继续迭代)
               │
               └── 迭代结束后 ──────────────
                   ├─ call local_stress           ← 最终应力
                   ├─ call contact_state(1)       ← 更新接触状态
                   ├─ call state_and_stiff_2021   ← 状态变量更新
                   ├─ call gpvarupdate            ← 高斯点变量存储
                   │
                   ├── 大变形更新 ──────────────
                   │   ├─ [Blarge==1] call update_coord_blarge
                   │   └─ call modf_element_local_direction
                   │
                   └── 输出 ──────────────
                       ├─ [istep%noutn==0] call outputres      ← 节点结果
                       ├─ [istep%noutf==0] call out_full_write ← 全量输出
                       ├─ [istep%noutf==0] call OUT_GID_WRITE  ← GiD格式
                       └─ [istep%nresta==0] call resta_read_write(-1) ← 重启
```

---

## 五、time_dependent 详细流程 (Fem.f90:8470-9408)

动力分析 Newmark 隐式时间积分，与 static_U 的关键差别:

```
SUBROUTINE time_dependent
│
├─ 与 static_U 相同的三层循环 (iincs → istep → iiter)
│
├─── 关键差别 ─────────────────────────────
│
│   1. 质量矩阵:
│      call mcmatrx('U')                  ← 组装质量矩阵
│      call mcmatrx('W')                  ← 流体质量
│
│   2. 耦合矩阵:
│      call hmatrx('W')                   ← 热传导矩阵
│      call upwcouple                     ← U-P/U-W 耦合矩阵
│      call COUPLE_ASSEMBLE               ← 耦合项组装
│
│   3. 时间积分:
│      call predict                       ← Newmark 预测 (用 beeta1, beeta2)
│      call varupdate                     ← 更新 result_zero, result_first, result_second
│      (result_first = 速度, result_second = 加速度)
│
│   4. 吸收边界:
│      call assemble_absorb_solid
│      call eload_absorb_solid
│      call assemble_interface_fluid_solid
│
│   5. 地震荷载:
│      earthquake_curve → fachv 放大系数
│      call modf_inpwav                   ← MIF 输入波场
│
│   6. 其它:
│      call algort                        ← 算法控制 (应力更新算法)
│      call porepr                        ← 孔隙水压力
│      call heat_internal1                ← 热源计算
│      call liquifaction_judge            ← 液化判定
│
└─── 收敛检查与 static_U 相同
```

---

## 六、其它分析路径

### frequency_analysis (Fem.f90:9806-9938)
- 频率域分析，使用**复数**数组 (stforw, toforw)
- 组装刚度+质量+阻尼，求特征值
- 无时间推进，求稳态响应

### explicit (Fem.f90:9941-10112)
- 显式时间积分 (中心差分法)
- 使用**集中质量** rmid (对角矩阵)
- **不需要求解器**: a = (F_ext - F_int) / M
- 时间步长受 CFL 条件约束

### response_spectrum (Fem.f90:10118)
- 反应谱分析
- 模态叠加法

---

## 七、关键全局数据结构

### 自由度系统
```
mdofn = 最大每节点自由度数 (可达8+)
cdofn = 压缩后实际活跃的自由度数
lmdofn(i) = 第i个自由度的压缩编号 (0=不活跃)
nodfn(cdofn, npoin) = 节点→全局方程号映射
ntotv = 总自由度数

示例: 2D U-P 问题
  mdofn=8, lmdofn=[1,2,0,0,0,0,0,3]
  cdofn=3, 每节点3个自由度 (ux, uy, p)
```

### 单元数据
```
element(ielem)%field(ifield)%lnods_f → 节点编号
element(ielem)%estif(nevab,nevab)    → 单元刚度矩阵
element(ielem)%eload(nevab)          → 单元荷载向量
element(ielem)%egaus(ifield)%sigma   → 高斯点应力
element(ielem)%egaus(ifield)%epsilon → 高斯点应变
element(ielem)%egaus(ifield)%stran   → 应变增量
element(ielem)%egaus(ifield)%stres   → 应力增量
element(ielem)%egaus(ifield)%gpvar   → 高斯点状态变量
```

### 单元组
```
group(igroup)%index    → 单元类型编号 (3=T3, 5=Q4, 10=B8 等)
group(igroup)%matno    → 材料编号
group(igroup)%nelgroup → 组内单元数
group(igroup)%list(:)  → 单元编号列表
group(igroup)%fieldid  → 场类型 ('U', 'P', 'T', 'W')
appear(igroup)         → 当前块是否活跃 (分级施工)
```

### 接触系统
```
gaps(igaps)%npairs           → 接触对数
gaps(igaps)%pairnode(2,npairs) → 接触节点对
gaps(igaps)%state(:)         → 接触状态 (开/闭/滑移)
gaps(igaps)%ctforce(:,:)     → 接触力
gaps(igaps)%dxyz(:,:)        → 相对位移
gaps(igaps)%rot(:,:,:)       → 局部坐标旋转矩阵
gaps(igaps)%kxyz(:,:,:)      → 接触刚度

gapb(igapb)%cmatrix          → 接触柔度矩阵
```

---

## 八、输入文件体系

所有文件以 `probn` 为基名：

| 文件 | 内容 |
|------|------|
| `probn.glb` | 全局控制参数 (problem definition) |
| `probn.cor` | 节点坐标 |
| `probn.ele` | 单元连接关系 |
| `probn.pre` | 边界条件 (约束、预定位移) |
| `probn.mat` | 材料参数 |
| `probn.loa` | 荷载定义 (体力、面力、集中力) |
| `probn.tem` | 温度荷载 |
| `probn.man` | 主控文件 (增量步参数) |
| `probn.sol` | 求解器参数 |
| `probn.ctt` | 接触数据 |
| `probn.ifs` | 流固界面 |
| `probn.ftr` | 温度传递 |
| `probn.act` | 分级施工 (appear_process) |
| `probn.bar` | 钢筋数据 |
| `probn.bem` | 梁单元数据 |
| `probn.ini` | 初始条件 |

### man 文件格式 (控制增量步循环)
```
"Block 1"                          ← text
nincs                              ← 增量数
miter ditime noutn noutf nstep inc_step nresta cwater Qstatic
toler_force toler_var(1:mdofn)
... (重复 nincs 组)
"Block 2"
...
```

---

## 九、刚体分析 (详见 rigid-body-reliability.md)

```
static_rigid_1 (Fem.f90:6217-6424)
│
├─ 当 block_stab>=1 且 ebody==0 时调用
├─ 物理问题: 大坝/挡墙/边坡刚体块滑动稳定性
│
├─ 与 static_U 的关键差别:
│   1. 自由度: 3(ndimn-1)个刚体DOF (2D:tx,ty,θz; 3D:tx,ty,tz,θx,θy,θz)
│   2. 接触求解: solve_ctt_rigid (Solver.f90:3793-4745)
│   3. 状态更新: state_and_stiff_rigid_2021 (Solver.f90:5263-5432)
│   4. 安全系数: call safety_factor
│
└─ block_stab 值:
    0 = 标准可变形分析
    1 = 刚体块运动学 (接触面可变形)
    2 = 混合 (刚体+普通节点DOF)
```

---

## 十、可靠度分析 (详见 rigid-body-reliability.md)

```
STATIC_U_reli (Fem.f90:4708-5364)       ← relis==1, block_stab==0
static_rigid_reli (Fem.f90:5957-6214)   ← relis==1, block_stab>=1
│
├─ 方法: 一阶可靠度方法 (FORM)
├─ 搜索最可能失效点 (MPP)
├─ 随机变量: 摩擦角φ, 内聚力c (可扩展)
├─ 分布: 正态, 对数正态, 极值分布
│
├─ 核心子程序:
│   ├─ DANGLI: 分布变换 (原始→等效正态)
│   ├─ RI3: Rackwitz-Fiessler设计点更新
│   ├─ betaindex: 可靠度指标 β = Φ⁻¹(1-Pf)
│   └─ stab_rcandgy_reli: 极限状态函数 G = 抗力 - 荷载
│
└─ 系统可靠度 (sysrelis>0):
    ├─ 串联系统 (任一失效)
    ├─ 并联系统 (全部失效)
    └─ 窄域法 (一般系统)
```

---

## 十一、参数反演 (详见 back-analysis.md)

```
Bparameter 分支:
├─ [-1,-2] 验证: parameter_back_analysis_verify (Fem.f90:1561-1630)
│   └─ Monte Carlo随机采样验证
│
├─ [1,2] & balgor<=1: parameter_back_analysis (Fem.f90:668-1024)
│   ├─ Intel MKL DTRNLSP (Levenberg-Marquardt)
│   ├─ balgor=0: 有限差分Jacobian
│   └─ balgor=1: 解析灵敏度 (dudx, Fem.f90:4291-4390)
│
├─ [1,2] & balgor==2: trust_region_back_analysis (Fem.f90:1027-1276)
│   ├─ 信赖域+BFGS Hessian近似
│   └─ 双狗腿法子问题求解
│
├─ [3] rigid_dis_back_analysis (Fem.f90:1362-1463)
│   └─ 最小二乘分解: u_total = u_rigid + u_elastic
│
└─ [4] nodal_value_back_analysis (Fem.f90:1465-1559)
    └─ 从采样点观测反演网格节点值

nbackf/=0 (在static_U内部调用):
├─ nbackdT==0: back_analysis (Fem.f90:2368-2955)
│   └─ 边界位移监测，刚体/弹性分解
└─ nbackdT/=0: back_d_analysis (Fem.f90:2957-3400+)
    └─ 位移-温度耦合监测
```

---

## 十二、关键控制参数分支总结

| 参数 | 值 | 影响 |
|------|-----|------|
| `type_problem` | 'Q' | 准静态 → static_U 系列 |
| | 'F','D' | 动力 → time_dependent/explicit |
| | 'E' | 地震反应谱 |
| | 'W' | 频率分析 |
| `type_solver` | 'PARDISO' | Intel MKL PARDISO 直接法 |
| | 'JPCG' | 预条件共轭梯度 |
| | 'PROFILE' | 天际线直接法 |
| | 'EXPLICIT' | 显式(无求解器) |
| `type_load` | 'LOAD2' | 两步加载(力学+温度) |
| | 'ARCLENGTH' | 弧长法 |
| | 'DISCONTROL' | 位移控制 |
| `type_nl` | 8 | 修正Newton (BFGS) |
| | 其它 | 标准Newton-Raphson |
| `nlayer` | 2 | 双层分析 |
| `mdofn` | 7 | U-P 耦合 |
| | 8 | U-Pw 耦合 |
| `ngaps` | >0 | 有接触分析 |
| `meshc` | 1,2 | 网格细化/粗化 |
| `rmesh` | >0 | 自适应重划分 |
| `Blarge` | 1 | 大变形 (更新坐标) |
| `nbackf` | >0 | 边界反分析监测 |
| `Bparameter` | 1,2 | 参数反演 (balgor选算法) |
| | 3 | 刚体位移分解 |
| | 4 | 节点值反演 |
| | -1,-2 | 反演验证 (Monte Carlo) |
| `balgor` | 0 | 有限差分Jacobian |
| | 1 | 解析灵敏度 |
| | 2 | 信赖域优化 |
| `block_stab` | 0 | 标准可变形分析 |
| | 1 | 刚体块稳定性 |
| | 2 | 混合刚体+可变形 |
| `relis` | 1 | FORM可靠度分析 |
| `sysrelis` | >0 | 系统可靠度 (串联/并联/窄域) |
| `ninit,kinit` | /=0,==2 | 初始应力 |
| `restart` | 1 | 从断点恢复 |

---

## 十三、为什么"模块化重构"如此困难

1. **全局状态无处不在**: `result_zero`, `tofor`, `stfor`, `element(:)`, `gaps(:)` 等
   几百个全局变量在 Fem.f90 的各个子程序间隐式传递，没有明确的接口。

2. **分支嵌套极深**: `type_problem` × `type_solver` × `type_load` × `type_nl` ×
   `nlayer` × `mdofn` × `ngaps` × `meshc` × `Bparameter` × ... 的组合空间巨大。

3. **隐式依赖**: 例如 `residu_f` 会修改 `element(ie)%egaus%sigma`，
   而 `conver_load` 读取 `stfor`，`force_internal` 从 `element(ie)%eload` 组装 `stfor`。
   这些依赖关系没有通过参数传递体现。

4. **文件I/O贯穿始终**: `mainunit` 在循环内部被读取，文件位置（指针）是隐式状态。
   不同的 `restart` 值会改变 rewind 和 skip 的逻辑。

5. **力学过程本身复杂**: 本构模型有状态变量(`gpvar`)、接触有状态(`state`)、
   蠕变有历史(`stran0`)、分级施工有出现/消失(`appear`)。这些都是物理上必须的。
