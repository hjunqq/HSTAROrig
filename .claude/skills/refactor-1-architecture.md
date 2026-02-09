# HSTAR Refactoring Guide - Part 1: Architecture Overview
# HSTAR重构指南 - 第一部分：架构总览

Use this skill when planning or discussing the overall refactoring strategy for HSTAR.

## Current Code Problems Analysis (现有代码问题分析)

### 1. God Module Anti-Pattern (上帝模块反模式)
```
global_var 模块被 9 个核心模块直接依赖
├── Material.f90
├── Stiff.f90
├── Residu.f90
├── Solver.f90
├── Prescrib.f90
├── Load.f90
├── Temper.f90
├── Output.f90
└── Fem.f90
```

**Problem**: Global.f90 contains 300+ global variables, making it impossible to:
- Test modules independently
- Understand data dependencies
- Refactor safely

### 2. Tight Coupling (紧耦合)
```fortran
! Current: Everything depends on everything
module materials
    use global_var  ! <-- 引入数百个不需要的变量
```

### 3. No Abstraction (缺乏抽象)
```fortran
! Current: Material types are hardcoded strings
if (material=='MOHRCOLUMB') then
    ! implementation
elseif (material=='DRUCKER') then
    ! different implementation
endif
```

### 4. I/O Embedded in Logic (输入输出嵌入逻辑)
```fortran
! Current: Direct file operations scattered everywhere
read(mainunit,*)text
read(mainunit,*)nincs
```

### 5. Poor Naming (命名不规范)
```fortran
! Current: Cryptic variable names
integer i0, ic, f0, e, nu
real estif(:,:), dmatx(:,:)
```

---

## Target Architecture (目标架构)

```
HSTAR-Modern/
├── src/
│   ├── core/                    # 核心数学和数据结构
│   │   ├── precision_mod.f90    # 精度定义
│   │   ├── tensor_mod.f90       # 张量运算
│   │   ├── sparse_matrix_mod.f90 # 稀疏矩阵
│   │   └── vector_mod.f90       # 向量操作
│   │
│   ├── mesh/                    # 网格模块
│   │   ├── node_mod.f90         # 节点定义
│   │   ├── element_mod.f90      # 单元基类
│   │   ├── mesh_mod.f90         # 网格容器
│   │   └── mesh_io_mod.f90      # 网格输入输出
│   │
│   ├── elements/                # 单元库
│   │   ├── element_interface.f90 # 单元抽象接口
│   │   ├── line2_mod.f90        # 2节点线单元
│   │   ├── tri3_mod.f90         # 3节点三角形
│   │   ├── quad4_mod.f90        # 4节点四边形
│   │   ├── tet4_mod.f90         # 4节点四面体
│   │   ├── hex8_mod.f90         # 8节点六面体
│   │   └── ...
│   │
│   ├── materials/               # 材料模块
│   │   ├── material_interface.f90 # 材料抽象接口
│   │   ├── elastic_mod.f90      # 弹性材料
│   │   ├── mohr_coulomb_mod.f90 # Mohr-Coulomb
│   │   ├── drucker_prager_mod.f90 # Drucker-Prager
│   │   ├── duncan_chang_mod.f90 # Duncan-Chang
│   │   └── ...
│   │
│   ├── boundary/                # 边界条件
│   │   ├── bc_interface.f90     # 边界条件接口
│   │   ├── dirichlet_bc_mod.f90 # 位移边界
│   │   ├── neumann_bc_mod.f90   # 力边界
│   │   └── contact_mod.f90      # 接触边界
│   │
│   ├── solver/                  # 求解器
│   │   ├── solver_interface.f90 # 求解器接口
│   │   ├── direct_solver_mod.f90 # 直接求解器
│   │   ├── iterative_solver_mod.f90 # 迭代求解器
│   │   └── newton_raphson_mod.f90 # NR迭代
│   │
│   ├── analysis/                # 分析类型
│   │   ├── static_analysis_mod.f90
│   │   ├── dynamic_analysis_mod.f90
│   │   ├── thermal_analysis_mod.f90
│   │   └── coupled_analysis_mod.f90
│   │
│   ├── io/                      # 输入输出
│   │   ├── input_parser_mod.f90 # 输入解析
│   │   ├── output_writer_mod.f90 # 输出写入
│   │   └── gid_interface_mod.f90 # GiD接口
│   │
│   └── app/                     # 应用层
│       └── fem_app_mod.f90      # 主程序封装
│
├── tests/                       # 测试
│   ├── unit/                    # 单元测试
│   ├── integration/             # 集成测试
│   └── benchmarks/              # 基准测试
│
├── examples/                    # 示例问题
├── docs/                        # 文档
└── CMakeLists.txt              # CMake构建

```

---

## Dependency Diagram (依赖关系图)

### Before (重构前)
```
         ┌──────────────┐
         │  global_var  │  (God Module)
         └──────┬───────┘
                │
    ┌───────────┼───────────┐
    │     ┌─────┼─────┐     │
    ▼     ▼     ▼     ▼     ▼
┌──────┐┌────┐┌────┐┌────┐┌──────┐
│Stiff ││Resi││Solv││Load││Output│
└──────┘└────┘└────┘└────┘└──────┘
    │           │           │
    └───────────┴───────────┘
         Circular deps!
```

### After (重构后)
```
┌─────────────────────────────────────────┐
│              Application                │
└──────────────────┬──────────────────────┘
                   │
┌──────────────────┼──────────────────────┐
│      ┌───────────┴───────────┐          │
│      ▼           ▼           ▼          │
│ ┌─────────┐ ┌─────────┐ ┌─────────┐     │
│ │Analysis │ │ Solver  │ │   I/O   │     │
│ └────┬────┘ └────┬────┘ └────┬────┘     │
│      │           │           │          │
│      └───────────┴───────────┘          │
│                  │                      │
│      ┌───────────┴───────────┐          │
│      ▼           ▼           ▼          │
│ ┌─────────┐ ┌─────────┐ ┌─────────┐     │
│ │Elements │ │Materials│ │Boundary │     │
│ └────┬────┘ └────┬────┘ └────┬────┘     │
│      │           │           │          │
│      └───────────┴───────────┘          │
│                  │                      │
│                  ▼                      │
│            ┌─────────┐                  │
│            │  Mesh   │                  │
│            └────┬────┘                  │
│                 │                       │
│                 ▼                       │
│            ┌─────────┐                  │
│            │  Core   │                  │
│            └─────────┘                  │
└─────────────────────────────────────────┘
```

---

## Key Refactoring Principles (核心重构原则)

### 1. Single Responsibility (单一职责)
每个模块只负责一件事

### 2. Dependency Inversion (依赖倒置)
高层模块不依赖低层模块，都依赖抽象接口

### 3. Open-Closed (开放封闭)
对扩展开放，对修改封闭 - 通过抽象接口添加新材料/单元

### 4. Interface Segregation (接口隔离)
接口要小而专一，不要大而全

### 5. Explicit Dependencies (显式依赖)
所有依赖通过参数传递，不使用全局变量
