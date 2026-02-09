# FEM Theory-Code Mapping (有限元原理-代码映射)

Use this skill when the user wants to understand the relationship between finite element theory and the HSTAR code implementation.

## Instructions

You are an expert in finite element analysis (FEM) and the HSTAR Fortran codebase. When the user asks about FEM concepts, help them understand:

1. **The theoretical principle** - Mathematical formulation and physical meaning
2. **The code implementation** - Where and how it's implemented in HSTAR
3. **The connection** - How the code translates the theory

## HSTAR Code Architecture Overview

### Core Modules and Their Roles

| Module | File | FEM Concept | Description |
|--------|------|-------------|-------------|
| `variable_types` | Vartype.f90 | 基础数据类型 | 定义精度(irk/ink)和基本类型 |
| `elements` | Elements.f90 | 单元库 | 形函数、高斯积分、单元类型定义 |
| `global_var` | Global.f90 | 全局变量 | 节点、单元、自由度、求解器参数 |
| `materials` | Material.f90 | 本构模型 | 材料属性、应力-应变关系 |
| `stiffness_matrix` | Stiff.f90 | 刚度矩阵 | 单元刚度、总体刚度组装 |
| `internal_force` | Residu.f90 | 内力计算 | 残差向量、内力向量 |
| `solver` | Solver.f90 | 线性求解 | PARDISO、PCG、Profile求解器 |
| `prescribed` | Prescrib.f90 | 边界条件 | 位移约束、力边界条件 |
| `applied_load` | Load.f90 | 荷载 | 外力向量、荷载增量 |
| `temperature` | Temper.f90 | 温度场 | 热传导、热应力 |
| `output` | Output.f90 | 后处理 | 结果输出、GiD格式 |

### Key FEM Concepts -> Code Mapping

#### 1. Shape Functions (形函数)
- **Theory**: N_i(ξ,η,ζ) - 插值函数，满足δ性质
- **Code**: `Elements.f90` - 各单元定义子程序
  - `l2_define` - 2节点线单元
  - `t3_define` - 3节点三角形
  - `q4_define` - 4节点四边形
  - `h4_define` - 4节点四面体
  - `b8_define` - 8节点六面体
- **Data**: `elkn(index)%ggaus(order_int)%shape` 存储形函数值

#### 2. Gauss Integration (高斯积分)
- **Theory**: ∫f(x)dx ≈ Σ w_i·f(x_i)
- **Code**: `Elements.f90` - `gauss_global` type
  - `weigp(:)` - 权重
  - `posgp(:,:)` - 积分点坐标
  - `ngaus` - 积分点数
- **Usage**: `DO igaus = 1, ngaus` 循环遍历积分点

#### 3. Strain-Displacement Matrix (B矩阵)
- **Theory**: ε = B·u, B = ∂N/∂x
- **Code**: `Stiff.f90` - `bmatx(:,:)` 数组
  - `cartd(:,:)` - 笛卡尔坐标导数 ∂N/∂x
  - 通过Jacobian变换: cartd = J^(-1) · deriv
- **Storage**: `element(ielem)%field(1)%bmatx` (部分单元)

#### 4. Constitutive Matrix (D矩阵/本构矩阵)
- **Theory**: σ = D·ε (弹性) 或 dσ = D_ep·dε (弹塑性)
- **Code**: `Material.f90` - 各种本构模型
  - `ELASTIC_ISOTROPIC` - 各向同性弹性
  - `DUNCANCHANG` - Duncan-Chang模型
  - `MOHRCOLOMB` - Mohr-Coulomb模型
  - `DRUCKER` - Drucker-Prager模型
- **Data**: `dmatx(:,:)` - 局部变量存储D矩阵

#### 5. Element Stiffness (单元刚度矩阵)
- **Theory**: K_e = ∫ B^T·D·B dV = Σ w_i·|J|·B^T·D·B
- **Code**: `Stiff.f90` - `SUBROUTINE STIFF_U`
  - `estif(:,:)` - 单元刚度矩阵
  - 关键代码: `estif = estif + djacb*weigp*matmul(transpose(bmatx),dbmat)`

#### 6. Global Assembly (总体刚度组装)
- **Theory**: K = Σ L_e^T·K_e·L_e (定位矩阵)
- **Code**: `Stiff.f90` / `Solver.f90`
  - `ldofs(:)` - 自由度编号(定位向量)
  - CSR格式存储: `iseq(:)`, `global_stiff1(:)`
- **Solver**: PARDISO (Intel MKL稀疏求解器)

#### 7. Internal Force (内力向量)
- **Theory**: F_int = ∫ B^T·σ dV
- **Code**: `Residu.f90` - `SUBROUTINE RESIDU_F`
  - `stfor(:)` - 内力向量
  - `tofor(:)` - 总力向量

#### 8. Newton-Raphson Iteration (牛顿迭代)
- **Theory**: K·Δu = F_ext - F_int, 迭代至收敛
- **Code**: `Fem.f90` - 主程序循环
  - `iiter` - 迭代次数
  - `deltafi(:)` - 位移增量
  - `toler_force` - 收敛容差

#### 9. Time Integration (时间积分)
- **Theory**: Newmark-β法, θ法
- **Code**: `Fem.f90` / `Global.f90`
  - `beeta1, beeta2` - Newmark参数
  - `theta1` - θ参数
  - `ditime` - 时间步长

### Supported Element Types (支持的单元类型)

| Index | Name | Description |
|-------|------|-------------|
| 1 | L2 | 2-node line (杆单元) |
| 2 | L3 | 3-node line (二次杆) |
| 3 | T3 | 3-node triangle (三角形) |
| 4 | T6 | 6-node triangle (二次三角形) |
| 5 | Q4 | 4-node quad (四边形) |
| 6 | Q8 | 8-node quad (二次四边形) |
| 7 | H4 | 4-node tetrahedron (四面体) |
| 8 | H10 | 10-node tetrahedron |
| 9 | B8 | 8-node hexahedron (六面体) |
| 10 | B20 | 20-node hexahedron |
| 20 | B2 | 2-node beam (梁单元) |
| 22 | P4 | Point-to-point contact |

### Material Models (本构模型)

- **ELASTIC_ISOTROPIC** - 线弹性
- **MOHRCOLOMB** - Mohr-Coulomb (岩土)
- **DRUCKER** - Drucker-Prager
- **DUNCANCHANG** - Duncan-Chang (非线性弹性)
- **GOODMAN** - 节理单元
- **CONTACT** - 接触单元
- **CREEP** - 蠕变模型 (icreep=1,2,3,4)

## How to Use This Skill

When asked about FEM concepts, follow this pattern:

1. **Explain the theory** with equations
2. **Locate the code** - specific file:line
3. **Show key variables** and their meanings
4. **Trace the data flow** through the program

Example queries:
- "如何计算单元刚度矩阵?" → 解释K_e公式 + Stiff.f90代码
- "Mohr-Coulomb模型在哪里实现?" → Material.f90相关子程序
- "荷载是如何施加的?" → Load.f90 + 增量步骤
- "如何处理边界条件?" → Prescrib.f90 + 自由度约束
