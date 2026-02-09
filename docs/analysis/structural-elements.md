# HSTAR 结构单元分析

## 一、单元类型索引映射

**位置:** `Elements.f90:128-153` (`kinddefine` 子程序)

| 索引 | 类型 | 名称 | 节点数 | 维度 | 场数 | 说明 |
|------|------|------|--------|------|------|------|
| 1 | 线性2节点 | l2 | 2 | 1 | 1 | 杆/桁架单元 |
| 2 | 线性3节点 | l3 | 3 | 1 | 1 | 高阶线单元 |
| 3 | 三角形3节点 | t3 | 3 | 2 | 1 | 平面应变/应力 |
| 4 | 三角形6节点 | t6 | 6 | 2 | 1 | 高阶三角形 |
| 5 | 四边形4节点 | q4 | 4 | 2 | 1 | 双线性四边形 |
| 6 | 四边形8节点 | q8 | 8 | 2 | 1 | 二次Serendipity |
| 7 | 四面体4节点 | h4 | 4 | 3 | 1 | 线性四面体 |
| 8 | 四面体10节点 | h10 | 10 | 3 | 1 | 二次四面体 |
| 9 | 六面体8节点 | b8 | 8 | 3 | 1 | 三线性六面体 |
| 10 | 六面体20节点 | b20 | 20 | 3 | 1 | 二次六面体 |
| 11 | T6-C3耦合 | t6c3 | 6 | 2 | 2 | U-P耦合三角形 |
| 12 | Q8-C4耦合 | q8c4 | 8 | 2 | 2 | U-P耦合四边形 |
| 13 | H10-C4耦合 | h10c4 | 10 | 3 | 2 | U-P耦合四面体 |
| 14 | B20-C8耦合 | b20c8 | 20 | 3 | 2 | U-P耦合六面体 |
| 15 | T3-C3耦合 | t3c3 | 3 | 2 | 2 | U-P耦合T3 |
| 16 | Q4-C4耦合 | q4c4 | 4 | 2 | 2 | U-P耦合Q4 |
| 17 | H4-C4耦合 | h4c4 | 4 | 3 | 2 | U-P耦合H4 |
| 18 | B8-C8耦合 | b8c8 | 8 | 3 | 2 | U-P耦合B8 |
| 19 | L2-C2耦合 | l2c2 | 2 | 1 | 2 | U-P耦合线单元 |
| 20 | **梁2节点** | b2 | 2 | 1 | 1 | 2D/3D梁单元 |
| 21 | **梁-接触** | b2c2 | 2 | 1 | 2 | 梁+接触耦合 |
| 22 | **板4节点** | p4 | 4 | 2 | 1 | DKT型板壳单元 |
| 23 | **棱柱6节点** | pr6 | 6 | 3 | 1 | 三棱柱/楔形体 |
| 24 | **棱柱耦合** | pr6c6 | 6 | 3 | 2 | 棱柱U-P耦合 |
| 25 | **钢筋/弹簧** | steel | 2 | 1 | 1 | 钢筋粘结/弹簧(2006) |
| 26 | **薄膜** | thin_film | 4 | 2 | 1 | 接触薄膜层(2023) |

---

## 二、梁单元 (索引20, 21)

### 2.1 单元定义

**位置:** `Elements.f90:886-924`

```
b2: nnode=2, ndimn=1, nrfields=1, nr_intrules=0 (解析积分)
b2c2: nnode=2, ndimn=1, nrfields=2 (位移场+接触场)
```

### 2.2 自由度

**2D梁:** 每节点3个自由度 (ux, uy, θz)，单元共6个自由度
**3D梁:** 每节点6个自由度 (ux, uy, uz, θx, θy, θz)，单元共12个自由度

### 2.3 截面属性

**定义:** `Material.f90:170-175` (`geometry_property` 类型)

```fortran
type geometry_property
    real(irk) Aera      ! 截面面积
    real(irk) J         ! 扭转惯性矩
    real(irk) Iy        ! y方向截面二阶矩
    real(irk) Iz        ! z方向截面二阶矩
    real(irk),pointer:: rotlg(:,:)    ! 局部-全局旋转矩阵
    integer(ink) ipd
    integer(ink),pointer::point_direct(:)
end type geometry_property
```

### 2.4 刚度矩阵公式

**位置:** `Stiff.f90:450-565`

**理论:** Euler-Bernoulli 经典梁理论

#### 2D梁局部刚度矩阵 (6×6)

```
K_local = | EA/L    0           0         -EA/L    0           0        |
          | 0       12EI/L³    -6EI/L²    0      -12EI/L³    -6EI/L²   |
          | 0      -6EI/L²      4EI/L     0        6EI/L²     2EI/L    |
          |-EA/L    0           0          EA/L    0           0        |
          | 0      -12EI/L³     6EI/L²    0       12EI/L³     6EI/L²   |
          | 0      -6EI/L²      2EI/L     0        6EI/L²     4EI/L    |
```

其中: `iea = E*Aera/L`, `iiy = E*Iy/L`

**代码 (Stiff.f90:469-478):**
```fortran
estifm(1,1)= iea
estifm(1,4)=-iea; estifm(4,4)=iea
estifm(2,2)= 12*iiy/dl**2
estifm(2,3)=-6*iiy/dl
estifm(2,5)=-12*iiy/dl**2
estifm(2,6)=-6*iiy/dl
estifm(3,3)= 4*iiy
estifm(3,5)= 6*iiy/dl
estifm(3,6)= 2*iiy
estifm(5,5)= 12*iiy/dl**2
estifm(5,6)= 6*iiy/dl
estifm(6,6)= 4*iiy
```

#### 3D梁局部刚度矩阵 (12×12)

**代码 (Stiff.f90:493-516):**

在2D基础上增加:
- **扭转刚度:** `itj = G*J/L` (剪切模量G × 扭转惯性矩J / 长度)
- **双向弯曲:** y-z两个平面分别有Iy, Iz对应的弯曲刚度
- `estifm(4,4) = itj; estifm(4,10) = -itj; estifm(10,10) = itj`

#### 坐标变换

**代码 (Stiff.f90:540-542):**
```fortran
estif = estifm .x. trot           ! K_local × T
estifm = estif
estif = transpose(trot) .x. estifm  ! T^T × K_local × T
```

变换矩阵 `trot` 由单元方向余弦构成:
- 2D: 6×6 分块对角矩阵 (3个2×2旋转块 + 角度不变)
- 3D: 12×12 分块对角矩阵 (4个3×3旋转块)

### 2.5 质量矩阵

**位置:** `Stiff.f90:2635-2755` (`beam_mass` 子程序)

#### 集中质量 (type_mass=0)

```fortran
estifm(1,1) = 0.5*ρA*L     ! 节点1 轴向
estifm(2,1) = 0.5*ρA*L     ! 节点1 横向
estifm(4,1) = 0.5*ρA*L     ! 节点2 轴向
estifm(5,1) = 0.5*ρA*L     ! 节点2 横向
```

#### 一致质量 (type_mass=1)

2D一致质量矩阵 (Stiff.f90:2684-2697):
```
M = ρAL/420 × | 140   0     0    70   0     0   |
               | 0     156  -22L  0    54    13L |
               | 0    -22L   4L²  0   -13L  -3L² |
               | 70    0     0    140  0     0   |
               | 0     54   -13L  0    156   22L |
               | 0     13L  -3L²  0    22L   4L² |
```

3D一致质量矩阵 (Stiff.f90:2698-2730) 包含 Timoshenko 剪切修正:
```
estifm(2,2) = 13/35 + 6Iz/(5A·L²)    ! 包含剪切变形效应
estifm(3,3) = 13/35 + 6Iy/(5A·L²)
estifm(4,4) = twist/(3A)               ! 扭转惯量
```

### 2.6 应力恢复

**位置:** `Residu.f90:4723-4798`

- 2D: `nstre = 3` → (轴力N, 剪力Vy, 弯矩Mz)
- 3D: `nstre = 6` → (轴力N, 剪力Vy, Vz, 扭矩Mx, 弯矩My, Mz)
- 积分点: `ngaus = 1` (解析积分)

### 2.7 预应力钢筋 (STEEL_SP)

**代码 (Stiff.f90:479-486, 517-524):**

当 `material == 'STEEL_SP'` 时，刚度从 `element(ielem)%field(1)%kdiag(idimn)` 读取对角值，用于预计算或外部分析的刚度注入。

---

## 三、板壳单元 (索引22)

### 3.1 单元定义

**位置:** `Elements.f90:926-951`

```
p4: nnode=4, ndimn=2, nrfields=1, ngaus=16 (4×4高斯积分)
```

### 3.2 自由度

每节点5个自由度: u, v (面内), w (法向), θx, θy (转角)
单元总自由度: 4×5 = 20

### 3.3 刚度公式

**位置:** `Stiff.f90:1853-1977` (`stif_p4_st` 子程序)

**理论:** DKT型板壳单元 (薄板Kirchhoff理论)

#### 膜刚度 (面内, 8×8)

**代码 (Stiff.f90:1909-1932):**

面内位移 (u, v) 对应2D弹性刚度。4节点×2DOF = 8个自由度。

#### 弯曲刚度 (面外, 12×12)

**代码 (Stiff.f90:1867-1906):**

弯曲自由度 (w, θx, θy)，取决于:
```
D_bend = E·t³ / [12(1-ν²)]
```

单元尺寸:
```fortran
a = 0.5·√(Σ(x₂ᵢ-x₁ᵢ)²)    ! x方向半长
b = 0.5·√(Σ(x₄ᵢ-x₁ᵢ)²)    ! y方向半长
```

#### 组合刚度 (20×20 → 24×24局部)

**代码 (Stiff.f90:1935-1961):**

将膜和弯曲刚度组合:
- 行1-2, 7-8, 13-14, 19-20: 膜自由度 (u,v)
- 行3-6, 9-12, 15-18, 21-24: 弯曲自由度 (w, θx, θy)

#### 坐标旋转

**代码 (Stiff.f90:1964-1971):**
```fortran
estif = rotation^T × local_stif × rotation
```

### 3.4 质量矩阵

**位置:** `Stiff.f90:2758-2901` (`plate_mass` 子程序)

- 16点高斯积分 (与刚度相同)
- 支持集中质量 (type_mass=0) 和一致质量 (type_mass=1)

### 3.5 本构矩阵

**位置:** `Stiff.f90:5206-5233` (`ecmat_p4`)

```
膜部分 (1-3): 平面应力格式
弯曲部分 (4-6): D[4,4] = t³·E/[12(1-ν²)], ...
横向剪切 (7-8): D[7,7] = k₀·E/[2(1+ν)]  (k₀=1.2 剪切修正系数)
```

---

## 四、钢筋/弹簧单元 (索引25)

### 4.1 单元定义

**位置:** `Elements.f90:232-249`

```
steel: nnode=2, ndimn=1
```

### 4.2 刚度计算

**位置:** `Stiff.f90:570-670`

#### 粘结-滑移关系

**代码 (Stiff.f90:617):**
```fortran
call steel_bond_slip_relation(ikindks, abs(dgap1), ks(1,1),
                              tao, ftx, fcx, dx, strabar, coefMpa)
ks(1,1) = ks(1,1) * aera    ! 缩放面积
```

#### 横向刚度

```fortran
if (doubsig==2) then
    if (ndimn>=2) ks(2,2) = ktan1    ! 横向1
    if (ndimn==3) ks(3,3) = ktan2    ! 横向2
endif
```

#### 特殊功能

- `icpspring` 数组标识弹簧附着端 (Stiff.f90:577-583)
- `prot` 每节点旋转矩阵用于局部坐标系
- `pstrain` 应变历史用于硬化模型
- 间隙追踪 (开/闭状态)

---

## 五、薄膜接触单元 (索引26)

**位置:** `Elements.f90:953-975`

```
thin_film: nnode=4, ndimn=2, ngaus=4 (2×2高斯积分)
```

2023年新增，用于接触界面薄层模型。
使用 `ecmat_thin_film()` (Stiff.f90:5237-5253) 和 `gbmat_thin_film()` 计算。

本构为简化膜模型 (平面应力):
```
D = E/(1-ν²) × | 1   ν   0        |
                | ν   1   0        |
                | 0   0   (1-ν)/2  |
```

---

## 六、棱柱单元 (索引23, 24)

**位置:** `Elements.f90:978-1044`

```
pr6:   nnode=6, ndimn=3, nrfields=1, ngaus=6 (刚度), ngaus=6 (质量)
pr6c6: nnode=6, ndimn=3, nrfields=2, ncouple=1 (U-P耦合)
```

用于3D轴对称问题转化为3D几何形状时的楔形体单元。

---

## 七、弹性弹簧单元 (索引1 + ELASTIC_SPRING材料)

**位置:** `Stiff.f90:327-351`

当 nnode==2 且材料为 'ELASTIC_SPRING' 时:

```fortran
l0 = props(matno)%mechanical%solid%Elastic_Spring%l0
! 间隙依赖刚度:
dmatx = (b + 2*c*|dgap| + 3*d*dgap²) * l0 / thick
```

状态:
- 'OPEN': dgap > 1e-8 → dmatx = 0
- 'CLOSE': dgap ≤ 1e-8 → dmatx = 计算值

---

## 八、坐标系统

### 局部坐标系存储

**全局变量 (Global.f90:135):**
```fortran
real(irk),allocatable::prot(:,:,:)   ! prot(ndimn, ndimn, npoin) 每节点旋转矩阵
```

**单元存储:**
```fortran
real(irk),pointer::rotation(:,:)     ! rotation(ndimn, ndimn) 方向余弦
```

### 梁单元方向余弦

局部x轴: 沿梁轴线方向
```fortran
x_local = (coord(:,node2) - coord(:,node1)) / length
```
垂直方向由正交化构建完整坐标基。

---

## 九、总结对照表

| 属性 | 梁(20) | 板(22) | 钢筋(25) | 薄膜(26) |
|------|--------|--------|----------|----------|
| 节点数 | 2 | 4 | 2 | 4 |
| 每节点DOF | 3(2D)/6(3D) | 5 | 1-3 | 5 |
| 总DOF | 6(2D)/12(3D) | 20 | 2-6 | 20 |
| 积分 | 解析 | 4×4高斯 | 解析 | 2×2高斯 |
| 理论 | Euler-Bernoulli | DKT薄板 | 粘结-滑移 | 接触薄膜 |
| 截面属性 | A, Iy, Iz, J | 厚度t | 面积A | 厚度t |
