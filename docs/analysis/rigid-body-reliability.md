# HSTAR 刚体分析与可靠度分析

## 一、刚体分析 (Rigid Body Analysis)

### 1.1 核心控制变量

| 变量 | 位置 | 值 | 含义 |
|------|------|-----|------|
| block_stab | Global.f90:47 | 0 | 标准可变形分析 |
| | | 1 | 刚体块运动学 (3(ndimn-1)个DOF) |
| | | 2 | 混合: 刚体DOF + 普通节点DOF |
| ebody | Global.f90:48 | 0 | 单接触面追踪 |
| | | 1 | 多接触面追踪 |
| nrdof | Global.f90:445 | 3(2D)/6(3D) | 每个刚体块的刚体自由度数 |

### 1.2 分析路径选择 (process_analysis)

```
type_problem == 'Q' (准静态):
├─ relis==1, block_stab==0  → STATIC_U_reli      (可变形+可靠度)
├─ relis==1, block_stab>=1  → static_rigid_reli   (刚体+可靠度)
├─ block_stab>=1, ebody==0  → static_rigid_1      (刚体确定性)
└─ else                     → static_U            (标准静力)
```

### 1.3 static_rigid_1 (Fem.f90:6217-6424)

**物理问题:** 大坝/挡墙/边坡的刚体块滑动稳定性分析

**与 static_U 的关键差别:**

1. **自由度表示:** 不用每节点ndimn个DOF，而用3(ndimn-1)个刚体DOF
   - 2D: 3个 (tx, ty, θz — 两个平移+一个旋转)
   - 3D: 6个 (tx, ty, tz, θx, θy, θz)

2. **接触求解:** 使用 `solve_ctt_rigid` 而非 `solve_ctt`

3. **状态更新:** 使用 `state_and_stiff_rigid_2021` 包含:
   - Goodman接头非线性刚度
   - 薄层行为
   - 损伤演化追踪

4. **额外功能:**
   - Qstatic参数: 静力荷载模型
   - 安全系数计算: `call safety_factor`
   - 界面力输出: `call force_interface`

**算法流程:**

```
static_rigid_1
├── 读取增量步参数 (nincs, miter, ditime, nstep, ...)
├── 读取Qstatic参数 (如有)
└── 增量循环 (iincs)
    └── 步循环 (istep)
        ├── 初始化稳定分析 (iblks >= stab_matde)
        ├── 更新时间参数、边界条件、荷载
        └── 细分循环 (idiv = 1 to mdiv)
            └── 迭代循环 (iiter = 1 to miter)
                ├── call algort (全局系统组装)
                ├── if iiter==1: call predict
                ├── if ngaps/=0: call solve_ctt_rigid  ← 刚体接触求解
                ├── call varupdate
                ├── call conver_nodal_value (收敛检查)
                └── 收敛后:
                    ├── call state_and_stiff_rigid_2021  ← 刚体状态更新
                    └── 更新 dxyz0 (位移历史)
```

### 1.4 solve_ctt_rigid (Solver.f90:3793-4745)

**核心:** 求解带刚体运动学的接触问题

**DOF变换 (block_stab==1):**
```fortran
kkdimn = 3*(ndimn-1)           ! 刚体DOF数 (2D:3, 3D:6)
allocate(rdisp(3*(ndimn-1)))   ! 刚体位移
```

**迭代算法:**

```
1. 初始化接触力: ctfor = 0

2. 将节点位移转换到质心位移:
   result_zero(nodfn) → dislocal (接触坐标系)
   → resultx (刚体块坐标)

3. 组装接触刚度矩阵:
   对每个活跃接触对:
     global_stiff_bt += kxyz0 * fact (含旋转变换)
   PROFILE分解

4. 迭代求解 (标号20循环):
   a. 构建残差向量 (刚体位移增量 + 接触力变化 + 间隙)
   b. 求解: profile_ctt_rigid → resultm
   c. 提取刚体位移: rdisp = resultm(gapb%rldofs)
   d. 累积: gapb%rdisp_inc += rdisp
   e. 更新接触力: ctforce += kxyz * rdisp
   f. 反变换到节点力: ctfori = transpose(rot) * ctforl

   g. 收敛检查:
      dnorm = ||ctforcei - ctforce||
      tnorm = ||ctforce||
      if (dnorm < tol AND tnorm < tol): 收敛

5. 未收敛:
   a. 更新全局残差: rvector = ctfor + tofor - stfor
   b. 求解主FEM系统: call solve
   c. 更新 resultx
   d. 返回步骤4
```

### 1.5 state_and_stiff_rigid_2021 (Solver.f90:5263-5432)

**接触状态转换:**

| 状态 | 条件 | 行为 |
|------|------|------|
| 0 (开) | gap > 0 | 最小法向刚度 |
| 1 (闭) | gap ≤ 0, σn < ft | 法向刚度激活 |
| 2 (滑) | \|τ\| > c + σn·μ | 降低切向刚度 |

**Goodman接头刚度更新 (Solver.f90:5378-5417):**

```
R(i) = 1.0 - Rf · |τᵢ|/τmax
k(i,i) = kgdm(i) · (σn/Pa)^n · R(i)²

if block_stab==1: k(i+3,i+3) = k(i,i)  (扩展到旋转DOF)
if block_stab==0: k = 1/k (转为柔度)
```

### 1.6 刚体位移变量

```fortran
gapb(igapb)%rdisp_zero(:)      ! 基准位移 (通常=0)
gapb(igapb)%rdisp_inc(:)       ! 增量位移 (solve_ctt_rigid中累积)
gapb(igapb)%rdisp_inc0(:)      ! 上步增量
gapb(igapb)%rdisp_first(:)     ! 速度 (动力分析)
gapb(igapb)%rdisp_second(:)    ! 加速度 (动力: β·Δt²·a)
gapb(igapb)%mass_inertia(:)    ! 质量和惯量
    ! 2D: [mx, my, Iz]
    ! 3D: [mx, my, mz, Ix, Iy, Iz]
```

---

## 二、可靠度分析 (Reliability Analysis)

### 2.1 控制变量

```fortran
relis     ! 0=确定性分析, 1=FORM可靠度分析
sysrelis  ! 0=仅构件可靠度, >0=计算系统可靠度
```

从inp文件读取: `restart, relis, sysrelis, ADINA, Uopt_R, gamamax`

### 2.2 方法: 一阶可靠度方法 (FORM)

**基本思想:** 在标准正态空间中搜索最可能失效点 (MPP, Most Probable Point)

**可靠度指标:** β = 标准正态空间中原点到极限状态面的最短距离

### 2.3 STATIC_U_reli (Fem.f90:4708-5364)

**物理问题:** 可变形体结构的概率安全分析

**输入 (从stoc文件读取):**

```
nbeta     — 极限状态/失效模式数
nv        — 随机变量数
mkiter    — 搜索设计点最大迭代次数
ja(1:nv)  — 分布类型 (1=正态, 2=对数正态, 3=极值)
ee(1:nv)  — 均值
ss(1:nv)  — 标准差
cov(nv,nv) — 相关矩阵
```

**算法流程:**

```
STATIC_U_reli
├── 读取: nbeta, nv, mkiter
├── 读取: ja(:), ee(:), ss(:), cov(:,:)
├── 分配: betas(nbeta) (结果存储)
│
└── 对每个极限状态 ibeta = 1 to nbeta:
    ├── 初始化: SP=ss, EP=ee, XA=EP, YA=0, er=1e-3
    │
    └── 设计点搜索循环 (标号333):
        ├── iter = iter + 1
        │
        ├── call DANGLI(nv, EE, SS, XA, JA, SD, ED, SP, EP)
        │   └── 分布变换: 原始空间 → 等效正态空间
        │
        ├── 更新随机参数到接触属性:
        │   ├── gaps(igaps)%frict(ipairs) = xa(ivfri)  (摩擦角)
        │   └── gaps(igaps)%cohes(ipairs) = xa(ivcoh)  (内聚力)
        │
        ├── 执行完整FEM分析:
        │   └── 步循环: force_external → stiff_u → solve → varupdate
        │
        ├── 评估极限状态函数:
        │   └── call stab_rcandgy_reli(rc, xa, gy)
        │       └── gy = 抗力 - 荷载效应 (>0安全, <0失效)
        │
        ├── 更新设计点:
        │   └── call RI3(nv, GY, GA, RC, XA, YA, SD, ED, SP, EP, cov)
        │       └── Rackwitz-Fiessler变换 + 设计点更新
        │
        ├── 收敛检查: if |gy| > er AND iter < mkiter → 继续迭代
        │
        └── 计算可靠度指标:
            └── call betaindex(nv, rc, ep, sp, xa, cov, beta)
                └── β = μY^T·α / σY
```

### 2.4 static_rigid_reli (Fem.f90:5957-6214)

**物理问题:** 刚体块的概率滑动稳定性分析

与 STATIC_U_reli 结构相同，但:
- 使用 `solve_ctt_rigid` 代替完整FEM求解
- 计算更快 (不需要每次迭代组装/求解全局刚度)
- 随机变量: 接触面的摩擦系数和内聚力

### 2.5 极限状态函数

**stab_rcandgy_reli (Fem.f90:16707-16801):**

```
G(X) = Σ(抗滑力) - Σ(滑动力)

     = Σ(接触面摩擦力) - Σ(施加切向荷载)

G(X) > 0 → 安全
G(X) = 0 → 极限状态 (失效面)
G(X) < 0 → 失效
```

### 2.6 随机变量分布变换

**DANGLI (Fem.f90:5429-5480):**

| 分布 (ja) | 变换公式 |
|-----------|---------|
| 1 — 正态 | SP=ss, EP=ee (直接) |
| 2 — 对数正态 | VV=(ss/ee)², SP=xa·√ln(1+VV), EP=xa·(1+ln(ee)-ln(xa·√(1+VV))) |
| 3 — 极值 | AR=1.28255/ss, AK=ee-0.5772/AR, Q=exp(-AR·(xa-AK)), SP和EP由CDF变换 |

### 2.7 设计点更新

**RI3 (Fem.f90:5388-5427) — Rackwitz-Fiessler变换:**

```
1. 含相关性的灵敏度: rc(i) = Σ(xc(j)·cov(i,j)·SP(j))
2. 归一化: GG = √(Σrc²)
3. 方向余弦: ga(i) = -rc(i)/GG
4. 更新标准正态空间: Y1 = Σ(ya·ga) + GY/GG
5. 反变换: xa(i) = ya(i)·SP(i) + EP(i)
```

### 2.8 可靠度指标计算

**betaindex (Fem.f90:5367-5386):**

```
β = μY^T · α / σY

其中:
  μY = 等效正态空间的均值向量
  α  = 灵敏度因子 (方向余弦)
  σY = 等效正态空间的标准差
```

### 2.9 系统可靠度

**system_reliability (Fem.f90:5483-5651):**

当 sysrelis > 0 时:

1. **串联系统** (任一构件失效→系统失效):
   ```
   Pf_series = Π Pf_i
   β_series = Φ⁻¹(1 - Pf_series)
   ```

2. **并联系统** (所有构件失效→系统失效):
   ```
   relat_ave = 平均相关系数
   β_parallel = β_ave · √(n/(1 + relat_ave·(n-1)))
   ```

3. **窄域法** (一般系统):
   - 计算相关矩阵: relat(i,j) = ga(:,i) · ga(:,j)
   - 识别相关构件 (relat > 0.7)
   - 多元正态概率积分

**输出:**
- 串联系统失效概率和可靠度指标
- 并联系统可靠度指标
- PNET系统可靠度
- 概率上下界

### 2.10 可随机化的参数

基于代码分析 (Fem.f90:4827-4839):

```fortran
do igaps = 1, ngaps
    ivcoh = gaps(igaps)%ivcoh    ! 内聚力在随机变量数组中的索引
    ivfri = gaps(igaps)%ivfri    ! 摩擦角在随机变量数组中的索引
    if(ivfri>0) gaps(igaps)%frict(ipairs) = xa(ivfri)
    if(ivcoh>0) gaps(igaps)%cohes(ipairs) = xa(ivcoh)
end do
```

**主要随机变量:**
- 接触面摩擦系数 (φ)
- 接触面内聚力 (c)
- 可扩展到弹性模量、荷载等

---

## 三、分析组合总结

| 组合 | 问题类型 | 应用场景 |
|------|---------|---------|
| block_stab=0, relis=0 | 标准FEM | 可变形土石结构 |
| block_stab=1, ebody=0, relis=0 | 刚体块静力 | 大坝/挡墙滑动稳定 |
| block_stab=1, ebody=1, relis=0 | 多面刚体 | 节理岩体分析 |
| block_stab=1, relis=1 | 刚体+概率 | 大坝可靠度 (摩擦/内聚力不确定) |
| block_stab=0, relis=1 | 可变形+概率 | 边坡可靠度 (材料不确定) |
