# HSTAR 重构策略 V2 — 基于深度代码分析的重新设计

## 零、为什么 V1 方案不可行

### V1 方案的核心假设

原始 12 周路线图（`migration-roadmap.md`）假设:
1. HSTAR 可以被分解为标准的 6 层架构
2. Newton-Raphson 循环可以统一为一个 `newton_driver`
3. 单元和材料可以通过 Factory 模式逐个替换
4. 12 周内可以完成首批物理迁移

### 深度分析揭示的现实

通过对 Fem.f90（~18000行）、Global.f90（~5000行）、Stiff.f90、Residu.f90、Solver.f90、
Material.f90、Elements.f90 的全面分析，发现:

**1. 不存在"一个" Newton-Raphson 循环**

实际存在 **10+ 个不同的分析变体**，每个都有自己的三层循环（增量→步→迭代）:

| 变体 | 行号 | 关键差异 |
|------|------|----------|
| `static_U` | 3560-4246 | 标准准静态 |
| `static_U_P` | — | U-P 耦合，额外 7 自由度 |
| `static_U_Pw` | — | U-Pw 耦合，额外 8 自由度 |
| `static_rigid_1` | 6217-6424 | 刚体 DOF 替代节点 DOF |
| `STATIC_U_reli` | 4708-5364 | FORM 可靠度包裹完整 FEM |
| `static_rigid_reli` | 5957-6214 | 刚体+可靠度 |
| `time_dependent` | 8470-9408 | Newmark 隐式，含质量/阻尼/耦合 |
| `explicit` | 9941-10112 | 中心差分，集中质量，无求解器 |
| `frequency_analysis` | 9806-9938 | 特征值，复数数组 |
| `back_analysis` | 2368-2955 | 边界监测，嵌入 N-R 内 |
| `back_d_analysis` | 2957-3400 | 位移-温度耦合监测 |

这些变体共享约 70% 的代码（刚度组装、残差计算、收敛检查），但在关键处有不同的:
- 自由度定义（节点 DOF vs 刚体 DOF vs 耦合 DOF）
- 矩阵类型（刚度 vs 刚度+质量+阻尼 vs 复数矩阵）
- 求解路径（直接法 vs 无求解器 vs 特征值）
- 外部包裹（可靠度搜索循环 vs 反演优化循环）

**"统一 newton_driver"** 需要处理所有这些变体的组合，比新写一个求解器更复杂。

**2. 接触不是"单元"，而是一个独立的求解子系统**

接触力学拥有:
- 独立的数据结构: `gap_group` (~30 个指针数组), `gap_block_group`
- 独立的状态机: open/contact/sliding 三态转换
- 独立的求解器: `solve_ctt`, `solve_ctt_rigid` (Solver.f90:2171-4745)
- 独立的柔度矩阵构建: `forAdirect` (Fem.f90:6724-6984)
- 三种本构模型: Goodman (Janbu 刚度), FCM (5 种软化律), MCJOINT
- 接触力→全局力的传递: `ctfor_to_tofor` (在 9 处被调用)

接触系统与主求解器之间是**双向耦合**：主求解器的位移结果影响接触状态，
接触力又参与全局平衡方程。这不是一个可以用 Factory 模式替换的组件。

**3. 可靠度和反演是"元分析"**

- FORM 可靠度: 在设计点搜索循环内反复调用完整的 `static_U` 或 `static_rigid_1`
- 参数反演: Intel MKL DTRNLSP 的 RCI 循环内反复调用 `process_analysis()`
- 信赖域反演: 每次迭代需要完整的前向分析

这意味着**任何对 static_U 的修改都会影响可靠度和反演**。
它们不是独立的层，而是嵌套的调用层次。

**4. 全局状态是"隐式参数传递系统"**

Global.f90 中的 ~300 个全局变量不仅仅是"共享数据"。它们构成了子程序之间的**隐式接口**:
- `residu_f` 修改 `element(ie)%egaus%sigma`
- `force_internal` 读取 `element(ie)%eload` 组装 `stfor`
- `conver_load` 读取 `stfor` 和 `tofor` 计算收敛
- `gpvarupdate` 保存高斯点状态变量到下一步

这些依赖关系**没有通过子程序参数传递体现**。如果把它们改为显式参数，
每个子程序的参数列表会膨胀到 20-50 个参数，代码可读性反而会下降。

**5. 文件 I/O 是隐式状态机**

`mainunit` 在循环内部被顺序读取，文件指针位置是隐式状态。
不同的 `restart` 值改变 rewind 和 skip 逻辑。
这意味着输入数据的读取顺序与分析循环结构**紧密耦合**。

---

## 一、重新定义目标

### 不做什么

1. **不做完整的架构替换** — 现有代码可以工作，而且经过了实际工程验证
2. **不追求"纯净"的分层设计** — HSTAR 的物理复杂性决定了组件之间必然有深度耦合
3. **不设定固定时间表** — 12 周对这个规模的程序是不现实的
4. **不同时改多个子系统** — 任何修改都需要回归验证

### 要做什么

1. **可维护性**: 把 18000 行的 Fem.f90 分解为可理解的模块
2. **可测试性**: 建立覆盖所有分析路径的自动化测试
3. **可扩展性**: 让添加新材料/新单元类型不需要修改 5 个文件
4. **渐进性**: 每一步改动都可以验证，可以回退

---

## 二、核心策略: 手术刀式重构 (Surgical Refactoring)

### 原则

**不是推倒重建，而是精准手术。**

就像外科手术一样:
1. 先做全面体检（已完成: 6 份分析文档）
2. 制定手术方案（本文档）
3. 每次只切一个部位
4. 每次手术后验证病人（程序）还在正常工作
5. 恢复后再做下一次手术

### 策略分解

```
Phase 0: 测试基础设施（安全网）
    ↓
Phase 1: 文件级分解（Fem.f90 → 多文件，逻辑不变）
    ↓
Phase 2: 输入数据结构化（消除 mainunit 隐式状态）
    ↓
Phase 3: 本构模型提取（Material/Stiff/Residu 三角）
    ↓
Phase 4: 单元库现代化（Elements/Stiff/Residu 三角）
    ↓
Phase 5: 求解器封装（PARDISO/PROFILE/JPCG 统一接口）
    ↓
Phase 6: 分析流程模板化（共性/差异分离）
```

---

## 三、Phase 0: 测试基础设施

### 目标

建立覆盖所有主要分析路径的回归测试集，这是所有后续重构的**安全网**。

### 现状

已有 4 个基线算例 (C001-C004: static, fix, freq, dynamic)。
但这远远不够 — 需要覆盖关键分支组合。

### 需要覆盖的分析路径

| 编号 | 路径 | 控制参数 | 优先级 |
|------|------|----------|--------|
| T01 | 标准静力 | type_problem='Q', 无接触 | P0 (已有 C001) |
| T02 | 固定边界 | type_problem='Q', 约束 | P0 (已有 C002) |
| T03 | 频率分析 | type_problem='W' | P0 (已有 C003) |
| T04 | 动力分析 | type_problem≠'Q','E','W' | P0 (已有 C004) |
| T05 | 接触静力 | ngaps>0, block_stab=0 | P1 |
| T06 | 刚体稳定 | block_stab≥1, ebody=0 | P1 |
| T07 | MC 弹塑性 | material_1, criteria='MC' | P1 |
| T08 | Duncan-Chang | material_5 | P1 |
| T09 | 混凝土损伤 | material_4, icr=2 or 3 | P2 |
| T10 | 蠕变分析 | icreep≠0 | P2 |
| T11 | 大变形 | Blarge=1 | P2 |
| T12 | 分级施工 | appear_process 变化 | P2 |
| T13 | 弧长法 | type_load='ARCLENGTH' | P3 |
| T14 | 显式动力 | type_solver='EXPLICIT' | P3 |
| T15 | U-P 耦合 | mdofn=7 | P3 |

### 测试工具增强

```
migration/baseline/
├── cases.csv                  # 测试算例清单
├── extract_results.py         # 结果提取 (已有)
├── compare_results.py         # 结果对比 (已有)
├── tolerances.yaml            # 容差配置 (已有)
├── run_baseline.bat           # 一键运行 (已有)
├── results.csv                # 黄金基线 (已有)
│
├── test_cases/                # 新增: 分类测试算例
│   ├── static/                # T01-T02
│   ├── dynamic/               # T04, T14
│   ├── frequency/             # T03
│   ├── contact/               # T05-T06
│   ├── material/              # T07-T10
│   └── advanced/              # T11-T15
│
└── run_regression.bat         # 新增: 全量回归测试
```

### 验收标准

- P0 (4 例): 立即可用
- P1 (4 例): Phase 1 开始前完成
- P2 (4 例): Phase 3 开始前完成
- P3 (3 例): Phase 5 开始前完成

---

## 四、Phase 1: 文件级分解

### 目标

将 Fem.f90 从单个 18000 行文件分解为 ~15 个文件，**不改变任何逻辑**。

### 原则

1. 纯机械操作: move subroutine → new file, add `use global_var`
2. 不修改接口、不修改逻辑、不修改变量名
3. 每分解一个子程序，立即跑回归测试
4. 使用 Fortran `include` 或 `submodule` 机制

### 分解方案

```
HSTAR/
├── Fem.f90              # 瘦身后: 主程序 + process_analysis (≤800行)
├── Fem_static.f90       # static_U, static_U_P, static_U_Pw (~2000行)
├── Fem_dynamic.f90      # time_dependent, explicit (~2000行)
├── Fem_frequency.f90    # frequency_analysis, response_spectrum (~400行)
├── Fem_rigid.f90        # static_rigid_1 (~600行)
├── Fem_reliability.f90  # STATIC_U_reli, static_rigid_reli (~1500行)
├── Fem_backanalysis.f90 # parameter_back_analysis, trust_region 等 (~2000行)
├── Fem_contact.f90      # contact_state, contact_pair_process 等 (~1000行)
├── Fem_init.f90         # global_data, material_set 等初始化 (~2000行)
├── Fem_output.f90       # outputres, out_full_write 等 (~1500行)
├── Fem_mesh.f90         # mesh_refine, modify_coord 等 (~1000行)
├── Fem_util.f90         # 小工具子程序集合 (~1500行)
├── Global.f90           # 不变
├── Elements.f90         # 不变
├── Material.f90         # 不变
├── Stiff.f90            # 不变
├── Residu.f90           # 不变
├── Solver.f90           # 不变
├── Prescrib.f90         # 不变
├── Load.f90             # 不变
├── Temper.f90           # 不变
├── Output.f90           # 不变
└── Level.f90            # 不变
```

### 关键技术点

**方法 A: `include` 方式** (最安全)
```fortran
! Fem.f90 (修改后)
module fem_procedures
    use global_var
    implicit none
contains
    include 'Fem_static.f90'
    include 'Fem_dynamic.f90'
    ! ...
end module
```
优点: 编译器看到完全相同的代码，零风险
缺点: 不是真正的模块化

**方法 B: `submodule` 方式** (推荐)
```fortran
! Fem.f90 — 父模块，声明接口
module fem_module
    use global_var
    implicit none
    interface
        module subroutine static_U()
        end subroutine
        module subroutine time_dependent()
        end subroutine
        ! ...
    end interface
end module

! Fem_static.f90 — 子模块，包含实现
submodule(fem_module) fem_static
contains
    module subroutine static_U()
        ! ... 原始代码不变 ...
    end subroutine
end submodule
```
优点: 真正的编译分离，修改子模块不重编译父模块
缺点: 需要 ifx 支持 submodule (已支持)

### 验收标准

- 所有回归测试 (P0+P1) 通过
- Fem.f90 从 ~18000 行降至 <1000 行
- 每个分解文件 <2500 行

---

## 五、Phase 2: 输入数据结构化

### 问题

当前输入通过 `mainunit` 顺序读取，文件指针位置是隐式状态:
```fortran
! Fem.f90 中的典型模式
read(mainunit) nincs
do iincs = 1, nincs
    read(mainunit) miter, ditime, noutn, noutf, nstep, inc_step, nresta, cwater, Qstatic
    read(mainunit) toler_force, toler_var(1:mdofn)
    ! ... 用 miter, ditime 等控制这个增量步 ...
end do
```

问题:
1. 文件读取和计算逻辑混合在一起
2. 重启时需要 rewind + skip 来恢复文件位置
3. 无法独立测试某个增量步的行为

### 解决方案: 预读取 + 结构化存储

**Step 2.1: 定义分析控制数据结构**

```fortran
! analysis_data_mod.f90
module analysis_data_mod
    use precision_mod
    implicit none

    type :: increment_control
        integer(IP) :: miter         ! 最大迭代次数
        real(WP)    :: ditime        ! 时间步长
        integer(IP) :: noutn, noutf  ! 输出间隔
        integer(IP) :: nstep         ! 总步数
        integer(IP) :: inc_step      ! 步长增量
        integer(IP) :: nresta        ! 重启间隔
        integer(IP) :: cwater        ! 水荷载标志
        integer(IP) :: Qstatic       ! 静力场标志
        real(WP)    :: toler_force   ! 力收敛容差
        real(WP), allocatable :: toler_var(:)  ! 变量收敛容差
        real(WP)    :: coef_water    ! 水压系数 (cwater≠0)
    end type

    type :: load_block_data
        character(80) :: description
        integer(IP) :: nincs          ! 增量步数
        type(increment_control), allocatable :: increments(:)
    end type

    type :: analysis_control
        integer(IP) :: nblks          ! 总块数
        integer(IP) :: runblks        ! 实际运行块数
        type(load_block_data), allocatable :: blocks(:)
    end type
end module
```

**Step 2.2: 实现预读取**

```fortran
! 在 global_data 中，读取完所有控制参数后:
! 原来: read(mainunit) 在循环内
! 新的: call read_analysis_control(mainunit, analysis_ctrl)
!       然后循环内直接访问 analysis_ctrl%blocks(iblks)%increments(iincs)
```

**Step 2.3: 渐进替换**

不是一次性替换所有 `read(mainunit)`，而是:
1. 先实现预读取，把数据存入结构体
2. 在原来 `read(mainunit)` 处，改为从结构体读取
3. 保留原始 `read` 语句作为注释，方便对照
4. 逐步移除注释

### 验收标准

- 所有 `mainunit` 的顺序读取被预读取替代
- 重启逻辑简化为"从结构体中找到正确的 block/increment 索引"
- 回归测试全部通过

---

## 六、Phase 3: 本构模型提取

### 问题

当前本构模型的代码分散在三个文件中:

```
Material.f90: 数据类型定义 + 参数读取
    ↕
Stiff.f90: D矩阵计算 (ecmat, INVART, YIELDS, FLOWFQ, DEPMDL, PKPN)
    ↕
Residu.f90: 应力更新 (bkwd_euler, concrete_*, calcdlan*)
```

添加一个新材料模型需要修改**所有三个文件**，且需要了解它们之间的调用关系。

### 解决方案: 每个模型一个模块

```
next/src/materials/
├── material_interface_mod.f90    # 抽象接口
├── elastic_mod.f90               # 线性弹性
├── mohr_coulomb_mod.f90          # MC 弹塑性
├── drucker_prager_mod.f90        # DP 弹塑性
├── von_mises_mod.f90             # VM 弹塑性
├── duncan_chang_mod.f90          # 非线性弹性
├── concrete_damage_mod.f90       # 混凝土损伤 (icr=1-6)
├── goodman_joint_mod.f90         # 接头模型
├── cam_clay_mod.f90              # 剑桥模型
├── soil_pz_mod.f90               # Pastor-Zienkiewicz
├── creep_mod.f90                 # 蠕变 (幂律/Burgers/堆石)
└── material_factory_mod.f90      # 工厂函数
```

### 抽象接口

```fortran
module material_interface_mod
    use precision_mod
    implicit none

    ! 材料状态 (高斯点级别)
    type :: material_state
        real(WP) :: stress(6)         ! Voigt 应力
        real(WP) :: strain(6)         ! Voigt 应变
        real(WP) :: plastic_strain    ! 等效塑性应变
        real(WP), allocatable :: internal_vars(:)  ! 内变量
    end type

    ! 抽象材料接口
    type, abstract :: material_base
    contains
        procedure(calc_dmatrix_if), deferred :: calc_dmatrix
        procedure(stress_update_if), deferred :: stress_update
    end type

    abstract interface
        ! 计算本构矩阵 D
        subroutine calc_dmatrix_if(self, state, ndimn, dmatrix, ierr)
            import :: material_base, material_state, WP, IP
            class(material_base), intent(in) :: self
            type(material_state), intent(in) :: state
            integer(IP), intent(in) :: ndimn
            real(WP), intent(out) :: dmatrix(:,:)
            integer(IP), intent(out) :: ierr
        end subroutine

        ! 应力更新 (给定应变增量，更新应力和内变量)
        subroutine stress_update_if(self, state, dstrain, ndimn, ierr)
            import :: material_base, material_state, WP, IP
            class(material_base), intent(inout) :: self
            type(material_state), intent(inout) :: state
            real(WP), intent(in) :: dstrain(:)
            integer(IP), intent(in) :: ndimn
            integer(IP), intent(out) :: ierr
        end subroutine
    end interface
end module
```

### 迁移方法: 适配器模式

不是立刻替换 Stiff.f90 中的 `ecmat`、`YIELDS` 等，而是:

```fortran
! mohr_coulomb_mod.f90
module mohr_coulomb_mod
    use material_interface_mod
    implicit none

    type, extends(material_base) :: mohr_coulomb_material
        real(WP) :: sigma0, hardening, frict_angle, dilan_angle, ft
    contains
        procedure :: calc_dmatrix => mc_calc_dmatrix
        procedure :: stress_update => mc_stress_update
    end type

contains
    subroutine mc_calc_dmatrix(self, state, ndimn, dmatrix, ierr)
        ! ... 从 Stiff.f90 中的 ecmat + YIELDS + FLOWFQ 提取 ...
    end subroutine

    subroutine mc_stress_update(self, state, dstrain, ndimn, ierr)
        ! ... 从 Residu.f90 中的 bkwd_euler 提取 ...
    end subroutine
end module
```

**关键: 新模块和旧代码并行运行**

```fortran
! 在 Stiff.f90 的 ecmat 中添加:
#ifdef USE_NEW_MATERIALS
    call new_material%calc_dmatrix(state, ndimn, dmatrix, ierr)
#else
    ! ... 原始代码 ...
#endif
```

用预处理器控制切换，确保随时可以回退。

### 迁移顺序

1. 线性弹性 (最简单，验证框架)
2. Mohr-Coulomb (最常用的弹塑性)
3. Duncan-Chang (非线性弹性)
4. Goodman 接头 (接触系统依赖)
5. 混凝土损伤 (6 个变体)
6. 蠕变模型
7. PZ 系列 (最复杂)

### 验收标准

每迁移一个模型:
- 新旧代码在相同输入下产生 bit-identical 结果
- 新模块有独立单元测试
- 回归测试全部通过

---

## 七、Phase 4: 单元库现代化

### 问题

当前单元库的形函数在 Elements.f90，刚度在 Stiff.f90，内力在 Residu.f90。
26 种单元类型用 `group(igroup)%index` 索引 (1-26) 分支。

### 解决方案

类似 Phase 3，每种单元一个模块:

```
next/src/elements/
├── element_interface_mod.f90     # 抽象接口
├── tri3_mod.f90                  # T3 三角形
├── tri6_mod.f90                  # T6 三角形
├── quad4_mod.f90                 # Q4 四边形
├── quad8_mod.f90                 # Q8 四边形
├── hex8_mod.f90                  # B8 六面体
├── hex20_mod.f90                 # B20 六面体
├── tet4_mod.f90                  # H4 四面体
├── tet10_mod.f90                 # H10 四面体
├── beam2d_mod.f90                # 二维梁
├── beam3d_mod.f90                # 三维梁
├── plate_mod.f90                 # 板壳
├── contact_element_mod.f90       # 接触单元 (薄层等)
└── element_factory_mod.f90       # 工厂
```

### 抽象接口

```fortran
type, abstract :: element_base
    integer(IP) :: nnodes   ! 节点数
    integer(IP) :: ndofn    ! 每节点自由度
    integer(IP) :: ngauss   ! 高斯点数
contains
    procedure(shape_if), deferred :: shape_functions
    procedure(bmatrix_if), deferred :: b_matrix
    procedure(stiffness_if), deferred :: element_stiffness
    procedure(internal_force_if), deferred :: element_internal_force
end type
```

### 迁移顺序

1. Q4 (2D, 最常用)
2. B8 (3D, 最常用)
3. T3, T6
4. H4, H10
5. Q8, B20
6. 梁/板 (特殊刚度矩阵)
7. 接触单元 (依赖接触子系统)

---

## 八、Phase 5: 求解器封装

### 问题

当前 `solve` 子程序通过全局变量 `operation` 和 `type_solver` 控制:
```fortran
operation='SET';       call solve   ! 初始化
operation='FACTORIZE'; call solve   ! 分解
operation='SOLVE';     call solve   ! 求解
```

`type_solver` 分支: PARDISO / PROFILE / JPCG / EXPLICIT

### 解决方案

```fortran
module solver_interface_mod
    implicit none

    type, abstract :: solver_base
    contains
        procedure(setup_if), deferred :: setup
        procedure(factorize_if), deferred :: factorize
        procedure(solve_if), deferred :: solve
        procedure(cleanup_if), deferred :: cleanup
    end type
end module

! pardiso_solver_mod.f90
type, extends(solver_base) :: pardiso_solver
    ! PARDISO handle, parameters
contains
    procedure :: setup => pardiso_setup
    procedure :: factorize => pardiso_factorize
    procedure :: solve => pardiso_solve
    procedure :: cleanup => pardiso_cleanup
end type
```

这一步相对独立，因为求解器的接口比较清晰。

---

## 九、Phase 6: 分析流程模板化

### 目标

这是最难的一步，也是最后一步。前面 5 个 Phase 都是在不改变主流程的前提下提取组件。
Phase 6 才开始触碰 `static_U`、`time_dependent` 等核心流程。

### 方法: 模板方法模式 (Template Method)

识别所有分析变体的**共性**和**差异点**:

```fortran
! analysis_template_mod.f90
module analysis_template_mod
    implicit none

    type, abstract :: analysis_base
    contains
        ! 模板方法 (共性流程)
        procedure :: run_analysis

        ! 钩子方法 (差异点，子类重载)
        procedure(init_step_if), deferred :: init_step
        procedure(assemble_if), deferred :: assemble_system
        procedure(predict_if), deferred :: predict
        procedure(update_if), deferred :: update_variables
        procedure(converge_if), deferred :: check_convergence
    end type

contains
    subroutine run_analysis(self)
        class(analysis_base), intent(inout) :: self
        ! 增量循环
        do iincs = 1, nincs
            ! 步循环
            do istep = ...
                call self%init_step()
                ! 迭代循环
                do iiter = 1, miter
                    call self%assemble_system()
                    call self%predict()
                    call solve()
                    call self%update_variables()
                    if (self%check_convergence()) exit
                end do
            end do
        end do
    end subroutine
end module
```

子类:
- `static_analysis`: 标准静力
- `dynamic_analysis`: Newmark 隐式 (额外: 质量/阻尼/加速度)
- `rigid_analysis`: 刚体 (额外: 刚体 DOF, solve_ctt_rigid)
- `explicit_analysis`: 显式 (无全局求解)

### 注意

Phase 6 的实施取决于 Phase 1-5 的完成程度。只有当组件（材料、单元、求解器）
都已经被提取和封装后，分析流程的模板化才有意义。

**如果 Phase 1-5 做得好，Phase 6 可能变成自然而然的结果而非强制的重构。**

---

## 十、接触系统特别策略

接触力学是 HSTAR 中最复杂的子系统之一，需要特殊处理。

### 阶段分解

**阶段 A (Phase 1 期间): 文件分离**
- `Fem_contact.f90`: contact_state, contact_pair_process, forAdirect
- Solver.f90 中的 solve_ctt, solve_ctt_rigid 暂不动

**阶段 B (Phase 3 之后): 接触本构提取**
- Goodman 接头 → goodman_joint_mod.f90
- FCM 损伤 → fcm_damage_mod.f90
- MCJOINT → 已在 MC 材料中

**阶段 C (Phase 5 之后): 接触求解器封装**
- solve_ctt → contact_solver_mod.f90
- solve_ctt_rigid → rigid_contact_solver_mod.f90
- 柔度矩阵构建 → compliance_matrix_mod.f90

**阶段 D (Phase 6 期间): 状态机重构**
- 将隐式的 state/state0/statei 转换为显式的状态机类型
- 接触力传递 (ctfor_to_tofor) 规范化

---

## 十一、可靠度和反演特别策略

### 关键约束

可靠度和反演都是**元分析**——它们包裹完整的 FEM 分析流程。
这意味着它们不能在 FEM 流程内部被拆解，而必须在外层处理。

### 处理方式

```
Phase 1: 分离为独立文件 (Fem_reliability.f90, Fem_backanalysis.f90)
Phase 6: 将它们实现为分析模板的外层包装:

reliability_analysis
  └── for each design point:
        └── static_analysis.run_analysis()  (或 rigid_analysis)

parameter_inversion
  └── for each DTRNLSP iteration:
        └── call analysis_template.run_analysis()
```

### 不需要重写的部分

- FORM 算法 (DANGLI, RI3, betaindex): 数学算法，独立性好，保持不变
- MKL DTRNLSP 接口: 黑盒调用，保持不变
- 信赖域 BFGS: 纯数学算法，保持不变

### 需要适配的部分

- 参数→材料属性的映射 (para_back → gaps%frict, material%E 等)
- 观测值提取 (Value_observ → result_zero)
- 灵敏度计算 (dudx → 需要访问刚度矩阵)

---

## 十二、已完成工作的定位

### Core 层 (next/src/core/)

已实现的 6 个模块 (`precision_mod`, `error_codes_mod`, `error_handling_mod`,
`tensor_mod`, `sparse_matrix_mod`, `gauss_quadrature_mod`) **保留并继续使用**。

它们在新策略中的位置:
- Phase 3-4 的新材料/单元模块将依赖这些 Core 模块
- Phase 5 的新求解器接口将使用 `sparse_matrix_mod`
- `gauss_quadrature_mod` 的数值必须与 legacy `Elements.f90:getgauss` 保持一致

### 基线测试 (migration/baseline/)

已有的 4 个测试算例和对比工具**保留并扩展**。
Phase 0 的工作就是在此基础上增加更多覆盖。

---

## 十三、实施路线图

```
                    Phase 0         Phase 1          Phase 2
                    测试基础设施     文件级分解        输入结构化
时间 →  ──────────┼───────────────┼────────────────┼──────────
已完成             P0(4例)已有      Core层已有
                   │               │
待做               P1(4例)新增      Fem.f90分解      mainunit
                   P2(4例)          为~15个文件       预读取
                   P3(3例)          纯机械操作         结构化


                    Phase 3          Phase 4          Phase 5
                    本构模型提取      单元库现代化      求解器封装
时间 →  ──────────┼────────────────┼────────────────┼──────────
                   Elastic           Q4               PARDISO
                   MC                B8               PROFILE
                   Duncan-Chang      T3/T6            JPCG
                   Goodman           H4/H10
                   Concrete          Q8/B20
                   Creep             Beam/Plate
                   PZ系列            Contact


                    Phase 6
                    分析流程模板化
时间 →  ──────────┼──────────
                   static_analysis
                   dynamic_analysis
                   rigid_analysis
                   explicit_analysis
                   (可靠度和反演作为外层包装)
```

### 每个 Phase 的估计工作量

| Phase | 复杂度 | 风险 | 说明 |
|-------|--------|------|------|
| 0 | 低 | 低 | 收集/创建测试算例 |
| 1 | 低 | 低 | 纯机械分解，不改逻辑 |
| 2 | 中 | 中 | 需要仔细保持读取顺序 |
| 3 | 高 | 中 | 本构模型逻辑复杂，但可逐个迁移 |
| 4 | 高 | 中 | 类似 Phase 3 |
| 5 | 中 | 低 | 求解器接口清晰 |
| 6 | 极高 | 高 | 触碰核心流程，必须前面 5 步都完成 |

---

## 十四、关键风险和缓解

### 风险 1: 全局变量依赖图不完整

**缓解**: Phase 1 (文件分解) 会自然暴露出每个子程序依赖哪些全局变量。
可以在分解时记录依赖表。

### 风险 2: 本构模型的隐式耦合

例如 `bkwd_euler` 内部调用 `YIELDS` → `FLOWFQ` → `INVART`，
这些函数通过参数传递很密，但也通过全局变量传递部分状态。

**缓解**: Phase 3 中，先不拆开这个调用链，而是整体迁移（把 INVART + YIELDS +
FLOWFQ + bkwd_euler 作为一个包迁入新模块）。

### 风险 3: 接触-主求解器耦合

接触力和位移是双向耦合的，不能简单地把接触"移出"主循环。

**缓解**: 接触系统的重构推迟到 Phase 5-6，在那之前只做文件分离。

### 风险 4: 数值精度漂移

本构模型的状态更新是数值积分过程，对浮点运算顺序敏感。
重构后即使逻辑相同，编译器优化可能导致微小数值差异。

**缓解**:
- 使用 `/fp:strict` 编译选项
- 容差设置为合理范围（位移 1e-10，应力 1e-8）
- 逐个模型对比，确认差异在可接受范围内

---

## 十五、与 V1 架构设计的关系

V1 的 `architecture-design.md` 定义的**目标架构**仍然有效。
本文档改变的是**到达目标的路径**:

| 方面 | V1 方案 | V2 方案 |
|------|---------|---------|
| 方向 | 从新架构出发，逐步迁入 | 从现有代码出发，逐步提取 |
| 第一步 | 实现新的 Core 层 | 分解 Fem.f90（不改逻辑） |
| 风险控制 | 新旧并行运行 + 对比 | 每一步都是原地重构 + 回归测试 |
| 全局变量 | 一开始就消除 | 先分文件，再逐步消除 |
| 时间估计 | 12 周 | 不设固定期限，按里程碑推进 |
| 前提条件 | 假设代码可分层 | 承认物理耦合，只要求可测试 |

**最终目标不变**: 一个分层的、可测试的、可扩展的 FEM 框架。
**路径改变**: 从"自顶向下设计+迁入"改为"自底向上提取+封装"。
