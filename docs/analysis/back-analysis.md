# HSTAR 参数反演与边界反演分析

## 一、反演系统总览

### 1.1 控制参数

| 参数 | 值 | 含义 |
|------|-----|------|
| Bparameter | -2 | 反演验证 (随机采样, 增量位移) |
| | -1 | 反演验证 (随机采样, 总位移) |
| | 0 | 无反演 |
| | 1 | 参数反演 (总位移) |
| | 2 | 参数反演 (增量位移) |
| | 3 | 刚体位移分解反演 |
| | 4 | 节点值反演 |
| balgor | 0 | 有限差分Jacobian |
| | 1 | 解析Jacobian (dudx) |
| | 2 | 信赖域优化 (BFGS) |
| nbackf | >0 | 边界观测组数 |
| nbackdT | 0 | 仅位移反演 |
| | ≠0 | 位移-温度耦合反演 |

### 1.2 分支选择 (主程序, Fem.f90:331-363)

```
Bparameter == -1 或 -2:
    → parameter_back_analysis_verify  (验证后stop)

Bparameter == 1 或 2:
    if balgor <= 1:
        → parameter_back_analysis    (Levenberg-Marquardt)
    if balgor == 2:
        → trust_region_back_analysis (信赖域)

Bparameter == 3:
    → rigid_dis_back_analysis       (刚体位移分解, 后stop)

Bparameter == 4:
    → nodal_value_back_analysis     (节点值反演, 后stop)

nbackf /= 0 (在static_U/static_U_Pw内部):
    if nbackdT == 0:
        → back_analysis             (边界位移监测)
    else:
        → back_d_analysis           (位移-温度耦合监测)
```

---

## 二、参数反演 (Bparameter=1,2, balgor≤1)

### 2.1 子程序: parameter_back_analysis (Fem.f90:668-1024)

**数学问题:**
```
min ||F(x)||² = Σ(观测值ᵢ - 计算值ᵢ(x))²
```

其中 x = [E, c, φ, K, ...] 为待反演参数向量。

### 2.2 求解器: Intel MKL DTRNLSP

使用 Intel MKL 的非线性最小二乘信赖域求解器 (Levenberg-Marquardt型)。

### 2.3 算法流程

```
parameter_back_analysis(N=Npara, M=Mvalue)
│
├── 初始化:
│   ├── DTRNLSP_INIT(handle, N, M, xvalue, eps, iter1, iter2, rs)
│   ├── DTRNLSP_CHECK(handle, ...) → 验证参数
│   └── RCI_REQUEST = 0, SUCCESSFUL = 0
│
├── RCI循环 (Reverse Communication Interface):
│   │
│   ├── DTRNLSP_SOLVE(handle, fvec, fjac, RCI_REQUEST) → 获取请求
│   │
│   ├── RCI_REQUEST == 1: 计算目标函数
│   │   ├── 参数变换: xvalue → 物理参数
│   │   │   (mode_transform=0: xvalue*factor, mode_transform=1: factor/xvalue)
│   │   ├── call process_analysis()  ← 完整前向FEM分析
│   │   ├── 计算残差: fvec(i) = 观测值(i) - 计算值(i)
│   │   └── 处理耦合观测 (jvalue≠0时取差值)
│   │
│   ├── RCI_REQUEST == 2: 计算Jacobian
│   │   ├── balgor==0: 有限差分
│   │   │   └── 对每个参数 ±δ 扰动, 重跑分析
│   │   │       FJAC(:,i) = -(F(x+δeᵢ) - F(x-δeᵢ)) / (2δxᵢ)
│   │   └── balgor==1: 解析灵敏度
│   │       └── FJAC(:,i) = -Value_observ(:)%dudx(i)
│   │
│   └── RCI_REQUEST < 0: 收敛/终止
│
├── 后处理 (收敛后):
│   ├── 灵敏度矩阵: fjac2 = J^T · J
│   ├── 误差协方差: Cov = (J^T·J)⁻¹
│   ├── 参数不确定度: σ(j) = √(mean(errx(j,i)²))
│   └── 输出反演结果
│
└── 清理: DTRNLSP_DELETE, MKL_FREE_BUFFERS
```

### 2.4 收敛准则

```fortran
eps(1) = 1e-6     ! 残差绝对容差
eps(2) = 1e-6     ! Jacobian相对容差
eps(3) = 1e-8     ! 梯度范数容差
eps(4) = 1e-8     ! 步长容差
iter1              ! 外循环最大迭代
iter2              ! 内循环最大迭代 (线性求解器)
```

### 2.5 数据结构

**待反演参数 (para_back):**
```fortran
para_back(Npara):
  imat           ! 材料组编号
  name           ! 参数名 (E, C, PHI, PERM, ...)
  factor         ! 缩放系数
  factor_inc     ! 有限差分扰动增量
  mode_transform ! 0=线性 (x*factor), 1=倒数 (factor/x)
```

**观测点 (Value_observ):**
```fortran
Value_observ(Mvalue):
  iblks, iincs, istep  ! 时间位置
  ivalue_point          ! 空间点索引
  idofn                 ! 自由度 (1-mdofn)
  value_measure         ! 观测/实测值
  value_computation     ! FEM计算值
  dudx(Npara)           ! 灵敏度 ∂u/∂参数
  jvalue                ! 耦合观测索引 (差值计算)
  ic                    ! 激活标志
```

---

## 三、信赖域反演 (Bparameter=1,2, balgor==2)

### 3.1 子程序: trust_region_back_analysis (Fem.f90:1027-1276)

### 3.2 信赖域子问题

```
min_{s} f_k + g_k^T s + 0.5 s^T B_k s    s.t. ||s|| ≤ Δ_k
```

其中:
- g_k = 梯度: `g_k(i) = -2·Σ(fvec · J(:,i))`
- B_k = Hessian近似 (BFGS更新)
- Δ_k = 信赖域半径

### 3.3 算法流程

```
trust_region_back_analysis(N, M)
│
├── 初始化:
│   ├── B_k = I (单位矩阵)
│   ├── x0 = xvalue (初始点)
│   └── delta_k = delta0 (初始信赖域半径)
│
└── 主循环 (iter_tr = 1 to mtter):
    ├── 前向分析: call process_analysis
    ├── 计算残差: fvec = 观测 - 计算
    ├── 计算梯度: g_k(i) = -2·Σ(fvec·J(:,i))
    │
    ├── 收敛检查: if ||g_k|| < eps → 结束
    │
    ├── 求解子问题: call solve_dx  ← 双狗腿法
    │   ├── s1 = -α·g_k  (最速下降步, α = ||g||²/(g^T·B·g))
    │   ├── s2 = -B⁻¹·g  (Newton步, Householder求解)
    │   └── 选择:
    │       ├── ||s1|| ≥ Δ: S = -(Δ/||g||)·g  (截断最速下降)
    │       ├── ||s2|| ≤ Δ: S = s2               (完整Newton)
    │       └── 否则: S = s1 + λ(s2-s1), ||S||=Δ (狗腿插值)
    │
    ├── 试探步: x_trial = x0 + s_k
    ├── 前向分析: call process_analysis
    │
    ├── 计算约化比:
    │   ├── 实际约化: dk_star = f(x_k) - f(x_trial)
    │   ├── 预测约化: dk_pre = -g_k·s - 0.5·s^T·B·s
    │   └── 比值: r_k = dk_star / dk_pre
    │
    ├── 步长接受/拒绝:
    │   ├── r_k ≥ η1: 接受 (x_new = x_trial)
    │   └── r_k < η1: 拒绝 (x_new = x_k)
    │
    ├── 信赖域半径更新:
    │   ├── r_k < η1:       Δ = 0.5·(0 + γ1·Δ)          (缩小)
    │   ├── η1 ≤ r_k ≤ η2:  Δ = 0.5·(γ1·Δ + Δ)          (维持)
    │   └── r_k > η2:       Δ = min(γ2·Δ, Δ_max)         (扩大)
    │
    └── BFGS更新 (如果 r_k ≥ η1):
        └── B_new = B + (y·y^T)/(s^T·y) - (B·s·s^T·B)/(s^T·B·s)
            其中 y = g_new - g_old, s = x_new - x_old
```

### 3.4 信赖域参数

```fortran
eta1, eta2       ! 步长接受阈值 (0 < η1 < η2 < 1)
gama1, gama2     ! 缩放因子 (0 < γ1 < 1 < γ2)
eps              ! 梯度范数容差
delta0           ! 初始半径
deltab           ! 最大半径
mtter            ! 最大迭代次数
```

---

## 四、灵敏度分析 (dudx)

### 4.1 子程序: dudx (Fem.f90:4291-4390)

**方法:** 直接微分法

对每个参数 ivar:

```
1. 识别受影响的材料组
2. 计算灵敏度荷载:
   如 name=='E':
     stfor_bar = (K/E) · u    (刚度对E的导数 × 当前位移)
3. 求解灵敏度问题:
   K · δu = -∂K/∂p · u
   → δu = ∂u/∂p  (位移对参数的灵敏度)
4. 在观测点提取:
   dudx(ivar) = dot_product(rintf, result(nodfn))
```

### 4.2 支持的灵敏度参数

| 参数名 | 含义 | 灵敏度计算 |
|--------|------|-----------|
| E | 弹性模量 | ∂K/∂E 从单元刚度 |
| PERM | 渗透系数 | 流动方程对渗透的导数 |
| C, PHI | 强度参数 | (由有限差分计算) |

---

## 五、刚体位移分解反演 (Bparameter=3)

### 5.1 子程序: rigid_dis_back_analysis (Fem.f90:1362-1463)

**物理问题:** 将实测接触面位移分解为刚体运动分量和弹性变形分量:

```
u_total = u_rigid + u_elastic
```

### 5.2 刚体参数

```fortran
npara = 3*(ndimn-1)    ! 2D: 3 (tx, ty, θz)
                       ! 3D: 6 (tx, ty, tz, θx, θy, θz)
```

### 5.3 算法

```
对每个时间步:
1. 收集观测点位移: observ(nmbpoint)

2. 构建插值矩阵: interp(nmbpoint, npara)
   (刚体位移影响系数, 来自 gapb%npdisp)

3. 最小二乘求解:
   interp^T · interp · rgdis = interp^T · observ
   → call householder(interpt, observt, rgdis)

4. 分解:
   u_rigid(i) = interp(i,:) · rgdis     (刚体分量)
   u_elastic(i) = observ(i) - u_rigid(i) (弹性分量)

5. 约束处理:
   固定DOF → 对角线加大罚数: interpt(i,i) *= 1e10
```

---

## 六、节点值反演 (Bparameter=4)

### 6.1 子程序: nodal_value_back_analysis (Fem.f90:1465-1559)

**物理问题:** 从场采样/积分点观测反演网格节点值

```
采样点观测 → 网格节点值
```

### 6.2 算法

```
对每个DOF类型:
  对每个时间步:
    1. 构建插值矩阵: interp(nmbpoint, npara)
       (空间插值基函数)
    2. 最小二乘: nodvar = (I^T·I)⁻¹·I^T·obs
    3. 赋值: result_zero(itotv) = nodvar(jpoin)
```

---

## 七、边界反演 (nbackf≠0)

### 7.1 back_analysis (Fem.f90:2368-2955)

**物理问题:** 分级施工过程中监测边界位移，将实测位移分解为刚体和弹性分量。

**数据结构 (Global.f90:351-356):**
```fortran
backf(igapbf):
  groupb         ! 接触组索引
  mdism          ! 位移观测点数
  wstep          ! 权重/步参数
  listp(mdism)   ! 节点索引
  listdim(mdism) ! DOF索引
  ic(mdism)      ! 激活标志
  dism(mdism, nstep) ! 实测位移时间序列
  wtime(mdism)   ! 时间权重
```

**集成方式:** 在 static_U 或 static_U_Pw 的Newton-Raphson循环内部调用。

**位移分解 (Fem.f90:2843-2893):**
```fortran
对每个观测点 kpoin:
  dispoint1  = Σ(result_zero_e · rintf)      ! 弹性部分
  dispoint1g = Σ((result_zero - result_zero_e) · rintf) ! 刚体部分
  qi = |dispoint1g / dispoint1|               ! 比值
```

**误差度量:**
```fortran
qerr = Σ(观测 - 计算)²     ! 残差
qabs = Σ(计算)²             ! 绝对值
RMSE = √(qerr / qabs)       ! 均方根误差
```

### 7.2 back_d_analysis (Fem.f90:2957-3400+)

与 back_analysis 结构相同，但处理位移-温度 (或位移-孔压) 耦合场。

---

## 八、反演验证 (Bparameter=-1,-2)

### 8.1 子程序: parameter_back_analysis_verify (Fem.f90:1561-1630)

**目的:** 用随机参数样本验证反演方法的有效性

### 8.2 算法

```
1. 生成高斯随机样本:
   para_stoch(Npara, nstoch)
   (使用 Intel MKL VSL随机数生成器)

2. 对每个样本 istoch = 1 to nstoch:
   a. xvalue = para_stoch(:, istoch)
   b. 参数变换
   c. call process_analysis  (前向FEM)
   d. 收集结果: Value_vc(观测点, 时步, istoch)
      Bparameter==-1: 存总位移
      Bparameter==-2: 存增量位移

3. 输出: 观测格式的计算结果
```

---

## 九、输入文件格式

### 9.1 .btl 文件 (反分析控制)

```
<注释行>
Npara                              ! 待反演参数数
<注释>
ipar imat name factor factor_inc mode_transform   [重复Npara行]
<注释>
eps iter1 iter2 rs jac_eps         ! DTRNLSP收敛参数
<注释: "待反演参数初始值">
x1 x2 ... xNpara                  ! 初始猜测值
```

信赖域额外参数 (balgor==2):
```
eta1 eta2 gama1 gama2 eps eta01 eta02 delta0 deltab mtter
```

### 9.2 .obs 文件 (观测数据)

```
<注释>
Nblks_pb                           ! 分级施工块数
<每块:>
  <块描述>
  nincs_pb                         ! 增量数
  <每个增量:>
    dtime_pb nstep_pb observ_pb begin_day end_day
<观测值数据>
```

---

## 十、方法对照表

| 方面 | 参数反演(balgor=0,1) | 信赖域(balgor=2) | 刚体位移(Bpar=3) | 节点值(Bpar=4) |
|------|---------------------|-----------------|-----------------|---------------|
| 问题 | min\|\|残差\|\|² | min\|\|残差\|\|²+BFGS | u分解 | 场重构 |
| 求解器 | DTRNLSP (MKL) | 手写信赖域 | 直接最小二乘 | 直接最小二乘 |
| Jacobian | 有限差分/解析 | 解析+BFGS近似 | N/A (约束矩阵) | N/A (插值矩阵) |
| Hessian | DTRNLSP内部 | BFGS公式 | N/A | N/A |
| 收敛速 | 快 (近解处二次) | 线性 (取决于BFGS) | 即时 (一步LS) | 即时 (一步LS) |
| 鲁棒性 | 中等 (依赖梯度) | 高 (显式TR保护) | 高 | 高 |
| 输出 | 参数+不确定度(σ) | 参数+不确定度 | 刚体分量 | 节点值 |
| 用途 | 材料参数识别 | 材料参数识别(需鲁棒) | 施工监测 | 场值重建 |

---

## 十一、数据流总图

```
输入文件:
  probn.inp   ──┐
  probn.btl   ──┼──→ 初始化 (Bparameter, balgor, nbackf)
  probn.obs   ──┘
                 │
        ┌────────┼────────┬──────────┬─────────┐
        ↓        ↓        ↓          ↓         ↓
     Bpar>0   Bpar<0   nbackf≠0   Bpar==3  Bpar==4
        │        │        │          │         │
   参数反演   验证分析  边界监测  刚体分解  节点反演
        │        │        │          │         │
   ┌─ balgor     │        │          │         │
   ├─0: 有限差分  │        │          │         │
   ├─1: 解析灵敏  │        │          │         │
   └─2: 信赖域   │        │          │         │
        │        │        │          │         │
        └───┬────┴────────┴──────────┴─────────┘
            │
            ↓
      反复调用 process_analysis (前向FEM)
            │
            ├──→ result_zero (位移)
            ├──→ deltafi (增量位移)
            └──→ dudx (灵敏度, balgor≥1时)
            │
            ↓
        收敛检查 → 输出结果
```
