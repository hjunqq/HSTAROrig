# HSTAR 本构模型分析

## 一、材料类型系统

### 1.1 类型概览 (Material.f90:1-200)

主容器: `solid_skeleton` 类型 (Material.f90:130-159)

| 材料类型 | 标识符 | 说明 |
|----------|--------|------|
| material_1 | ClassicalEP | 经典弹塑性 (MC/DP/VM/TC/MCC/DPC/MCJOINT) |
| material_2 | CamClay | 剑桥模型 |
| material_3 | SoilPZ | 广义塑性 (Pastor-Zienkiewicz) |
| material_4 | Concrete | 损伤混凝土模型 |
| material_5 | DuncanChang | 非线性弹性双曲线模型 |
| material_6 | Goodman | 接头/界面单元 |
| material_7 | Creep | 时间依赖变形 |
| material_8 | Elastic_Spring | 线性弹簧 |
| material_9 | Elastic_ep | 分段弹塑性 |
| material_10 | Plane_lowft | 低拉伸强度平面 |
| material_11 | Steel_ep | 钢材应力-应变曲线 |
| material_12 | Steel_sp | 钢材空间属性 |
| material_13 | SandPZ | 砂土广义塑性 |
| material_14 | ClayPZ | 粘土广义塑性 |
| material_15 | Scycl | 抗液化模型 |

---

## 二、弹性模型

### 2.1 线性各向同性弹性

**子程序:** `ecmat()` (Stiff.f90:4845-4936)

#### 1D
```
D = [E]
```

#### 2D 平面应力 (PS)
```
c = E/(1-ν²)
D = | c    cν   0         |
    | cν   c    0         |
    | 0    0    c(1-ν)/2  |
```

#### 2D 平面应变 (PE)
```
α = E(1-ν)/[(1+ν)(1-2ν)]
β = αν/(1-ν)
γ = α(1-2ν)/[2(1-ν)]
D = | α  β  0  β |
    | β  α  0  β |
    | 0  0  γ  0 |
    | β  β  0  α |
```

#### 3D
```
α = E(1-ν)/[(1+ν)(1-2ν)]
β = Eν/[(1+ν)(1-2ν)]
G = E/[2(1+ν)]
D(1:3,1:3) = α在对角线, β在非对角线
D(4:6,4:6) = G在对角线
```

### 2.2 Duncan-Chang 非线性弹性 (双曲线模型)

**类型定义:** `material_5` (Material.f90:48-54)

```fortran
type material_5
    character(2) model    ! 'EV', 'CR', 'EB'
    real(irk) Cohes, phi  ! 内聚力, 摩擦角
    real(irk) K, n        ! 体积模量参数
    real(irk) Rf          ! 破坏比 (0.7-0.9)
    real(irk) G, F, Vtf   ! EV模型参数
    real(irk) Kb, m, dphi ! EB模型参数
    real(irk) Pa, P0      ! 大气压, 参考围压
    real(irk) Nur, Kur    ! 卸载模量
    real(irk) phi_s, K_s  ! 湿化后参数
end type material_5
```

**应力依赖刚度:**
```
E_t = K·Pa·(σ₃/Pa)^n·[1 - Rf·(σ₁-σ₃)/((σ₁-σ₃)_f)]²
B = Kb·Pa·(σ₃/Pa)^m
```

**模型变体:**
- **EV:** 弹性体积变化模型
- **CR:** 帽模型变体
- **EB:** 弹塑性含剪胀

### 2.3 低拉伸强度平面 (PLANE_LOWFT)

**子程序:** `ecmat_lowft()` (Stiff.f90:4939-5036)

检查主拉应力，在开裂方向降低刚度:
```
if σ₁ ≥ ft₁ AND σ₂ ≥ ft₂:
    D → 1e-5·E (全面开裂)
elif σ₁ ≥ ft₁ XOR σ₂ ≥ ft₂:
    D → 单方向降低 (部分开裂)
else:
    D → 标准弹性 (未开裂)
```

### 2.4 板壳弹性 (P4)

**子程序:** `ecmat_p4()` (Stiff.f90:5206-5233)

```
膜部分 (1-3): 平面应力格式
弯曲部分 (4-6): D_b = E·t³/[12(1-ν²)]
横向剪切 (7-8): D_s = k₀·E/[2(1+ν)]  (k₀=1.2)
```

---

## 三、弹塑性模型

### 3.1 经典弹塑性框架

**类型定义:** `material_1` (Material.f90:16-20)

```fortran
type material_1
    character(20) criteria    ! MC, DP, VM, TC, MCC, DPC, MCJOINT
    real(irk) sigma0          ! 初始屈服应力/内聚力
    real(irk) hardening       ! 硬化模量 dσ/dεp
    real(irk) frict_angle     ! 内摩擦角 (度)
    real(irk) dilan_angle     ! 剪胀角 (度)
    real(irk) ft              ! 拉伸强度 (拉伸截断)
    real(irk) sigmat          ! 拉伸应力限 (MCJOINT)
    integer(ink) csigma0, cfrict, cdilan, cft, csigmat  ! 曲线索引
end type material_1
```

### 3.2 Mohr-Coulomb (MC)

**位置:** `Stiff.f90:7206-7212` (屈服面评估), `Stiff.f90:5746-6066` (流动法则)

**屈服面:**
```
f = σ_mean·sin(φ) + √(3J₂)·[cos(θ) - sin(θ)·sin(φ)/√3] - c·cos(φ) = 0
```

其中:
- φ: 内摩擦角
- θ: Lode角 (-π/6 ≤ θ ≤ π/6)
- c: 内聚力 (sigma0)
- J₂: 第二偏应力不变量

**流动法则:** 非关联 (用剪胀角ψ替代φ)
```
dεp = dλ · ∂Q/∂σ    其中Q用ψ角
```

**代码:**
```fortran
SNPHI = SIN(FRICT)
EQSTR = SMEAN*SNPHI + STEFF*(COS(THETA) - SIN(THETA)*SNPHI/ROOT3)
COHES = UNIAX + EPSTN*HARDS    ! 线性硬化
YVALU = COHES*COS(FRICT)
```

### 3.3 Drucker-Prager (DP)

**位置:** `Stiff.f90:7234-7240`

**屈服面:**
```
f = [6sin(φ)/(√3(3-sin(φ)))]·σ_mean + √(3J₂) - [6cos(φ)/(√3(3-sin(φ)))]·c = 0
```

**特点:**
- 比MC更光滑 (π平面无角点)
- 更好的数值稳定性
- 与MC的φ, c有相似的物理意义

**代码:**
```fortran
SNPHI = SIN(FRICT)
EQSTR = 6.0*SMEAN*SNPHI/(ROOT3*(3.0-SNPHI)) + STEFF
YVALU = 6.0*COHES*COS(FRICT)/(ROOT3*(3.0-SNPHI))
```

### 3.4 von Mises (VM)

**位置:** `Stiff.f90:7201-7205`

**屈服面:**
```
f = √(3J₂) - σy = 0
```

不依赖平均应力，适用于金属等材料。

```fortran
EQSTR = ROOT3 * STEFF
YVALU = UNIAX + EPSTN*HARDS
```

### 3.5 Tresca (TC)

**屈服面:**
```
f = 2·cos(θ)·√J₂ - σy = 0
```

### 3.6 Mohr-Coulomb + 拉伸截断 (MCC)

**修正屈服面:**
```
f = [MC准则] + α₀·σ₁ = 0

α₀ = (c·cos(φ) - ft) / ft - 0.5 - 0.5·sin(φ)
```

防止不合理的拉伸发展，适用于混凝土和岩石。

### 3.7 Drucker-Prager + 拉伸截断 (DPC)

类似MCC，但基于DP准则。

### 3.8 Mohr-Coulomb接头 (MCJOINT)

**位置:** `Stiff.f90:5643-5671` (INVART), `Stiff.f90:5815-5868` (FLOWFQ)

在局部接头坐标系中:
```
f = √(3J₂) + tan(φ)·σ_mean - c = 0
```

使用旋转矩阵转换到局部 (法向, 切向) 坐标系。

---

## 四、损伤模型

### 4.1 混凝土损伤模型

**类型定义:** `material_4` (Material.f90:31-40)

```fortran
type material_4
    real(irk) A, B, C, D    ! 损伤参数
    real(irk) Fc             ! 抗压强度
    real(irk) Gf             ! 断裂能
    real(irk) h              ! 特征网格尺寸
    real(irk) Ct             ! 拉/压强度比 (通常0.1)
    integer(ink) icr         ! 混凝土模型类型 (1-6)

    ! 拉伸损伤参数
    real(irk) at, bt, alfat
    real(irk) t1, t2, t3, t4, ft0, eft

    ! 压缩损伤参数
    real(irk) ac, bc, alfac
    real(irk) c1, c2, c3, c4, fc0, efc

    ! 其他
    real(irk) bb, dt, et0
end type material_4
```

### 4.2 icr变体

**icr=1: 一般脆性弹性断裂** (Residu.f90:3223)
- 简单拉伸/压缩损伤演化

**icr=2: 光滑损伤函数** (Residu.f90:3312)
- 多项式损伤演化 (Material.f90:674-714)
- 分离的拉伸和压缩分支
- 多项式系数 c1-c4 由拟合确定:
```
σ = E·ε·D(ε)    D为应变的多项式函数
```

**icr=3,5,6: 指数损伤** (Material.f90:660-672)
- 拉伸损伤:
```
ft(D) = (c·Ct·Fc)·exp(-ft·h·εp/Gf)
bb = 3/[ε₀·(2Gf·E/(h·ft²) - 1)]
```

### 4.3 损伤实现

**子程序位置:**
- `concrete_1()` — Residu.f90:3223 (类型1)
- `concrete_2()` — Residu.f90:3312 (类型2)
- `concrete_3()` — Residu.f90:3458 (类型3)
- `concrete_5()` — Residu.f90:3525 (类型5)
- `concrete_6()`, `concrete_6x()` — Residu.f90:3589-3811 (类型6)

**核心算法:**
```
1. 计算应变和有效应力
2. 评估损伤准则
3. 更新损伤参数 (0 ≤ D ≤ 1)
4. 降低刚度: D_tangent = (1-D)²·D_elastic
5. 更新应力: σ = (1-D)·σ_elastic
```

---

## 五、蠕变/粘弹性模型

### 5.1 蠕变类型系统

**类型定义:** `material_7` (Material.f90:68-73)

```fortran
type material_7
    integer(ink) nr           ! 组分数
    integer(ink) cvstrain     ! 应变控制标志
    real(irk) a, b            ! 幂律参数
    real(irk) alfas, bs, cs, ds  ! 岩石蠕变参数
    real(irk) m1, m2, m3      ! 应力依赖指数
    real(irk) Ek, etak, etam  ! Burgers模型参数
    real(irk),pointer:: c(:), d(:), f(:), k(:)  ! Kelvin链参数
end type material_7
```

**控制标志 icreep:**
- 0: 无蠕变
- 1: 幂律蠕变
- 2: Burgers模型 (多组分粘弹性链)
- 3: 堆石体蠕变
- 4: Burgers变体

### 5.2 幂律蠕变 (icreep=1)

**时间依赖模量演化:**
```
E(t) = E₀·(1 - a·t^b)
```

### 5.3 堆石体蠕变 (icreep=3)

**子程序:** `creep_strain_of_rock_fill()` (Residu.f90:~5650+)

**nr=3 (简化):**
```
εv(t) = bs·(p/Pa)
εq(t) = ds·(γ/(1-γ))
```

**nr=7 (完整):**
```
εv(t) = bs·(p/Pa)^m₁ + cs·(q/Pa)^m₂
εq(t) = ds·(γ/(1-γ))^m₃
```

蠕变增量:
```
Δε = alfas·εv·(1 - t_elapsed/t_total)·Δt
```

### 5.4 Burgers模型 (icreep=2,4)

**icreep=4 (简化Burgers):**
```
参数: Ek (弹性模量), etak (主粘度), etam (次粘度)
```

**icreep=2 (含Kelvin-Voigt链):**

多组分 (nr个弹簧/阻尼器):
```
ε_elastic = σ/Ek
ε_primary = (σ/ηk)·(1 - exp(-k(i)·t))     ! Kelvin-Voigt响应
ε_secondary = (σ/ηm)·t                      ! 阻尼器响应 (无界)

总应变: ε = ε_elastic + ε_primary + ε_secondary
```

**实现 (Residu.f90:5098-5126):**
```fortran
do ir = 1, nr    ! 对每个Kelvin组分
    cir = c(ir) + d(ir)·ptime^(-f(ir))
    qn += cir·(1 - exp(-k(ir)·Δt·dt/2))
end do
E_eff = 1 + qn·E    ! 修正模量
```

---

## 六、广义塑性模型

### 6.1 SoilPZ (Pastor-Zienkiewicz)

**类型定义:** `material_3` (Material.f90:26-29)

```fortran
type material_3
    integer(ink) ntest
    real(irk) d(24)    ! 24个参数
end type material_3
```

**核心概念:** 双屈服面 (加载面和卸载面)

**参数 d() 数组:**

| 索引 | 参数 | 含义 |
|------|------|------|
| 1 | PHIG | 加载摩擦角 |
| 2 | PHIF | 卸载摩擦角 |
| 3,4 | SIN(PHIG), SIN(PHIF) | 正弦值 |
| 5 | XMFC/XMGC | 应力比因子 |
| 6,7 | ALFAF, ALFAG | 屈服面形状参数 |
| 8 | PCUT | 压力截断 |
| 9,10 | HEV0, HES0 | 体积/剪切模量参数 |
| 12 | ICELS | 弹性控制 (0=双变, 1=K不变, 2=G不变, 3=双常) |
| 13,14 | BETA0, BETA1 | 指数形状参数 |
| 15 | H0 | 初始硬化参数 |
| 16 | GAMDM | 剪胀参数 |
| 20 | HU0 | 卸载体积模量 |
| 21 | GAMHU | 卸载硬化 |
| 22,23 | XMGC, XMFC | 临界状态参数 |
| 24 | EXPF | 指数因子 |

**核心子程序 (Stiff.f90):**
- `IVRMDL`: 应力不变量计算 (p, q, J₂, J₃, θ)
- `FNVMDL`: 加载/卸载方向判断
- `FAVMDL`: 屈服面法向量计算
- `DEPMDL` (Stiff.f90:7581+): 弹塑性刚度矩阵
- `RINMDL`: 回映更新塑性应变

**应力依赖属性:**
```
K(p) = HEV0·(p/PCUT)^β₀·(指数调整)
G(p) = HES0·(p/PCUT)^β₁
```

### 6.2 SandPZ (砂土广义塑性)

**类型定义:** `material_13` (Material.f90:106-111)

```fortran
type material_13
    integer(ink) ntest, humidification, pztype
    integer(ink),allocatable:: bline(:), eline(:)
    real(irk),allocatable:: sigmad(:)
    real(irk) d(24)
end type material_13
```

**扩展功能:**
- pztype: 模型类型 (11=与Duncan-Chang耦合)
- humidification: 湿度修正
  - ±2: 时间依赖湿化
  - 3: 循环/幅值依赖
- 额外参数 d(17-20): cohes, phi, p0, pa (湿化参数)

### 6.3 ClayPZ (粘土广义塑性)

**类型定义:** `material_14` (Material.f90:113-116)

```fortran
type material_14
    integer(ink) ntest
    real(irk) d(11)    ! 11个参数
end type material_14
```

**参数:**
| 索引 | 参数 | 含义 |
|------|------|------|
| 1 | Kevo | 体积模量 |
| 2 | Keso | 剪切模量 |
| 3 | Mg | 临界状态线斜率 |
| 4 | alfag | 表面形状 |
| 5 | H0 | 硬化 |
| 6 | expf | 指数因子 |
| 7,8 | beta0, beta1 | 形状参数 |
| 9 | nu | 泊松比 |
| 10 | pc0 | 初始固结压力 |
| 11 | icels | 弹性控制标志 |

### 6.4 Cam-Clay

**类型定义:** `material_2` (Material.f90:22-24)

```fortran
type material_2
    real(irk) Pc        ! 预固结压力
    real(irk) lamda      ! 塑性压缩指数
    real(irk) Mg         ! 临界状态线斜率 (p-q空间)
    real(irk) Mf         ! 破坏斜率 (常=Mg)
    real(irk) D0, D1     ! 剪胀参数
    real(irk) gamma      ! 硬化参数
end type material_2
```

---

## 七、应力不变量与屈服评估

### 7.1 INVART — 不变量计算

**位置:** `Stiff.f90:5629-5743`

**输出:**
- THETA: Lode角 (弧度, -π/6 ≤ θ ≤ π/6)
- STEFF: 有效偏应力 √J₂ 或 √(3J₂)
- SMEAN: 平均正应力 σ_mean = (σ₁+σ₂+σ₃)/3
- varj2: J₂ 第二偏应力不变量
- VARJ3: J₃ 第三不变量
- sint3: sin(3θ)

**3D公式:**
```
J₂ = [(σ₁-σ₂)² + (σ₂-σ₃)² + (σ₁-σ₃)²]/6 + τ₁₂² + τ₂₃² + τ₃₁²

J₃ = σ'₁·σ'₂·σ'₃ + 2τ₁₂·τ₂₃·τ₃₁ - σ'₁·τ²₂₃ - σ'₂·τ²₃₁ - σ'₃·τ²₁₂

STEFF = √J₂

sin(3θ) = -2.598·J₃/(J₂·STEFF)    [截断到 [-1, 1]]
THETA = (1/3)·asin(sint3)
```

### 7.2 YIELDS — 屈服准则评估

**位置:** `Stiff.f90:7137-7270`

计算等效应力(eqstr)和屈服值(yvalu):
- eqstr > yvalu → 屈服 (需要回映)
- eqstr ≤ yvalu → 弹性

### 7.3 FLOWFQ — 流动法则

**位置:** `Stiff.f90:5746-6066`

计算塑性流动方向:
- AVECT: 加载法向量 (∂f/∂σ)
- AVECQ: 卸载法向量 (∂g/∂σ, 塑性势)

**A向量分解:**
```
n = cons1·A₁ + cons2·A₂ + cons3·A₃

A₁ = ∂I₁/∂σ = [1,1,1,0,0,0]ᵀ         (体积方向)
A₂ = ∂√(3J₂)/∂σ = 偏应力分量缩放      (偏差方向)
A₃ = ∂θ/∂σ                             (Lode角导数)
```

---

## 八、应力更新算法

### 8.1 后退Euler积分

**位置:** `Residu.f90:2291-2498` (`bkwd_euler`)

```
步骤1: 弹性试应力
    σᵗʳⁱᵃˡ = σₙ + Dₑ·Δε

步骤2: 屈服检查
    f(σᵗʳⁱᵃˡ) ≤ 0 → 弹性 (直接返回)
    f(σᵗʳⁱᵃˡ) > 0 → 需要回映

步骤3: Newton-Raphson回映 (最多20次迭代)
    初始化: Δλ = 0

    迭代:
        n = ∂f/∂σ|σc            (流动方向)
        σc = σᵗʳⁱᵃˡ - Δλ·(Dₑ·n)   (回映到屈服面)
        εpc = εpn + Δλ·Q(n)       (更新等效塑性应变)

        fc = f(σc, εpc)
        if |fc| < 1e-5: 收敛

        Δλ_new = Δλ + fc/(H + nᵀ·Dₑ·n)   (Newton更新)

步骤4: 更新状态变量
    σₙ₊₁ = σc
    εpₙ₊₁ = εpc
    D_tangent = Dₑ - (Dₑ·n)⊗(nᵀ·Dₑ)/(H + nᵀ·Dₑ·n)   (一致切线刚度)
```

### 8.2 塑性乘子计算

**calcdlan1** (Residu.f90:4156-4178) — 一阶近似:
```
Δλ = f / (nᵀ·Dₑ·n + H)
```

**calcdlan2** (Residu.f90:4180-4213) — Newton-Raphson精化:
```
解: Q·x = r  [Householder分解]
Δλ_new = Δλ + Δλ_correction
```

### 8.3 法向量导数

**dadsig** (Stiff.f90:6600-7045):
- 计算法向量关于应力的导数
- 处理角点奇异条件
- ABTHE < 29° 时特殊处理

---

## 九、总结对照表

| 模型 | 标识 | 关键参数 | 塑性 | 损伤 | 蠕变 |
|------|------|----------|------|------|------|
| 线性弹性 | ELASTIC | E, ν | 无 | 无 | 无 |
| Duncan-Chang | DUNCANCHANG | K,n,φ,Rf | MC类 | 无 | 可选 |
| 经典EP | CLASSICALEP | σ₀,H,φ,ψ | MC/DP/VM/TC | 无 | 无 |
| 混凝土 | CONCRETE | Fc,Ct,Gf,h | 无 | icr=1-6 | 无 |
| Goodman接头 | GOODMAN | K₁,n,φ,kzz | MC | FCM可选 | 无 |
| Cam-Clay | CAMCLAY | Pc,λ,Mg | 帽 | 无 | 无 |
| SoilPZ | SOILPZ | 24参数 | 广义塑性 | 无 | 无 |
| SandPZ | SANDPZ | 20参数 | 广义塑性 | 无 | 无 |
| ClayPZ | CLAYPZ | 11参数 | 广义塑性 | 无 | 无 |
| 弹簧 | ELASTIC_SPRING | a,b,c,d,l0 | 无 | 无 | 无 |
| 幂律蠕变 | icreep=1 | a, b | 弹性 | 无 | 有 |
| Burgers | icreep=2,4 | Ek,ηk,ηm | 弹性 | 无 | 有 |
| 堆石蠕变 | icreep=3 | bs,cs,ds,m | 弹性 | 无 | 有 |

---

## 十、关键子程序索引

| 子程序 | 文件 | 行号 | 功能 |
|--------|------|------|------|
| ecmat | Stiff.f90 | 4845-4936 | 各向同性弹性矩阵 |
| ecmat_lowft | Stiff.f90 | 4939-5036 | 低拉伸弹性矩阵 |
| ecmat_p4 | Stiff.f90 | 5206-5233 | 板壳弹性矩阵 |
| INVART | Stiff.f90 | 5629-5743 | 应力不变量 |
| FLOWFQ | Stiff.f90 | 5746-6066 | 流动法则/法向量 |
| dadsig | Stiff.f90 | 6600-7045 | 法向量导数 |
| YIELDS | Stiff.f90 | 7137-7270 | 屈服准则评估 |
| DEPMDL | Stiff.f90 | 7581+ | SoilPZ弹塑性矩阵 |
| PKPN | Stiff.f90 | 6067-6200 | Goodman接头刚度 |
| bkwd_euler | Residu.f90 | 2291-2498 | 后退Euler应力更新 |
| calcdlan1/2 | Residu.f90 | 4156-4213 | 塑性乘子 |
| concrete_1~6 | Residu.f90 | 3223-3811 | 混凝土损伤模型 |
| material_set | Material.f90 | 210-1004 | 材料参数读取/初始化 |
