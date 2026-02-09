# HSTAR Refactoring Overview (重构总览)

Use this skill to get a quick overview of the HSTAR refactoring plan and current code issues.

## Current Code Issues Summary (现有代码问题汇总)

| 问题 | 严重程度 | 影响 | 解决方案 |
|------|---------|------|---------|
| 上帝模块 (Global.f90) | ★★★★★ | 无法单独测试模块 | 拆分为独立模块 |
| 300+ 全局变量 | ★★★★★ | 难以追踪数据流 | 封装到类型中 |
| 紧耦合 | ★★★★☆ | 修改一处影响全局 | 依赖注入 |
| 无抽象接口 | ★★★★☆ | 难以扩展新材料/单元 | 抽象类+工厂模式 |
| I/O嵌入逻辑 | ★★★☆☆ | 难以测试、难以更换格式 | 分离I/O层 |
| 变量命名不规范 | ★★★☆☆ | 难以理解代码 | 遵循命名规范 |
| 无测试框架 | ★★★☆☆ | 无法验证正确性 | 添加单元测试 |
| 老旧Fortran语法 | ★★☆☆☆ | 可维护性差 | 使用现代Fortran特性 |

## Refactoring Phases (重构阶段)

```
┌─────────────────────────────────────────────────────────────────┐
│                     REFACTORING ROADMAP                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Phase 1: Foundation (基础层)          [2-3 weeks]              │
│  ════════════════════════════════════                          │
│  ├── precision_mod.f90     精度定义                             │
│  ├── tensor_mod.f90        张量运算                             │
│  ├── sparse_matrix_mod.f90 稀疏矩阵                             │
│  └── gauss_quadrature_mod  高斯积分                             │
│                                                                 │
│  Phase 2: Mesh Layer (网格层)          [2 weeks]                │
│  ════════════════════════════════════                          │
│  ├── node_mod.f90          节点定义                             │
│  ├── element_interface.f90 单元抽象接口                         │
│  ├── mesh_mod.f90          网格容器                             │
│  └── mesh_io_mod.f90       网格I/O                             │
│                                                                 │
│  Phase 3: Material Layer (材料层)      [3 weeks]                │
│  ════════════════════════════════════                          │
│  ├── material_interface.f90 材料抽象接口                        │
│  ├── elastic_mod.f90        弹性材料                            │
│  ├── mohr_coulomb_mod.f90   Mohr-Coulomb                       │
│  ├── drucker_prager_mod.f90 Drucker-Prager                     │
│  └── state_manager_mod.f90  状态管理                            │
│                                                                 │
│  Phase 4: Solver Layer (求解器层)      [2 weeks]                │
│  ════════════════════════════════════                          │
│  ├── solver_interface.f90   求解器接口                          │
│  ├── pardiso_solver_mod.f90 PARDISO封装                        │
│  └── newton_raphson_mod.f90 NR迭代控制                         │
│                                                                 │
│  Phase 5: Integration (集成)           [2 weeks]                │
│  ════════════════════════════════════                          │
│  ├── fem_engine_mod.f90     FEM引擎                            │
│  ├── static_analysis_mod.f90 静力分析                          │
│  └── output_writer_mod.f90  输出系统                            │
│                                                                 │
│  Phase 6: Migration (迁移)             [4+ weeks]               │
│  ════════════════════════════════════                          │
│  ├── 创建适配层桥接新旧代码                                      │
│  ├── 逐步迁移单元类型                                           │
│  ├── 逐步迁移材料模型                                           │
│  └── 验证与基准测试                                             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Architecture Comparison (架构对比)

### Before (重构前)
```
                    ┌──────────────────┐
                    │    global_var    │  ← 300+ 全局变量
                    │   (God Module)   │
                    └────────┬─────────┘
           ┌─────────────────┼─────────────────┐
           │                 │                 │
      ┌────▼────┐       ┌────▼────┐       ┌────▼────┐
      │ Stiff   │       │ Residu  │       │ Solver  │
      │ .f90    │       │ .f90    │       │ .f90    │
      └────┬────┘       └────┬────┘       └────┬────┘
           │                 │                 │
           │    ┌────────────┴────────────┐    │
           │    │                         │    │
      ┌────▼────▼───┐               ┌─────▼────▼───┐
      │  Elements   │               │   Material   │
      │   .f90      │               │    .f90      │
      └─────────────┘               └──────────────┘

    问题: 循环依赖、全局状态、难以测试
```

### After (重构后)
```
    ┌─────────────────────────────────────────────────┐
    │               Application Layer                 │
    │  ┌─────────────┐  ┌─────────────────────────┐  │
    │  │ FEM Engine  │  │ Analysis Controllers    │  │
    │  └──────┬──────┘  └────────────┬────────────┘  │
    └─────────┼──────────────────────┼───────────────┘
              │                      │
    ┌─────────▼──────────────────────▼───────────────┐
    │               Service Layer                     │
    │  ┌─────────┐  ┌─────────┐  ┌─────────────────┐ │
    │  │ Solver  │  │ Assembly│  │ I/O Manager     │ │
    │  └────┬────┘  └────┬────┘  └────────┬────────┘ │
    └───────┼────────────┼────────────────┼──────────┘
            │            │                │
    ┌───────▼────────────▼────────────────▼──────────┐
    │               Domain Layer                      │
    │  ┌─────────┐  ┌─────────┐  ┌─────────────────┐ │
    │  │Elements │  │Materials│  │ Boundary Conds  │ │
    │  │(Factory)│  │(Factory)│  │                 │ │
    │  └────┬────┘  └────┬────┘  └────────┬────────┘ │
    └───────┼────────────┼────────────────┼──────────┘
            │            │                │
    ┌───────▼────────────▼────────────────▼──────────┐
    │               Core Layer                        │
    │  ┌─────────┐  ┌─────────┐  ┌─────────────────┐ │
    │  │  Mesh   │  │ Tensor  │  │ Sparse Matrix   │ │
    │  └─────────┘  └─────────┘  └─────────────────┘ │
    └─────────────────────────────────────────────────┘

    优势: 清晰分层、单向依赖、可独立测试
```

## Key Design Patterns (关键设计模式)

| 模式 | 应用场景 | 示例 |
|------|---------|------|
| **Factory** | 创建单元/材料对象 | `create_element('Q4')` |
| **Strategy** | 不同求解器算法 | `solver_base` 抽象类 |
| **Template Method** | 单元刚度计算流程 | `calc_stiffness()` |
| **Observer** | 收敛监控 | `convergence_monitor` |
| **RAII** | 资源管理 | `final` 析构函数 |
| **Adapter** | 新旧代码桥接 | `legacy_adapter_mod` |

## Recommended Reading Order (建议阅读顺序)

1. **refactor-1-architecture.md** - 了解目标架构
2. **refactor-5-best-practices.md** - 学习现代Fortran
3. **refactor-2-modules.md** - 核心模块设计
4. **refactor-3-io-testing.md** - 测试框架
5. **refactor-4-migration.md** - 实施迁移

## Quick Start Commands (快速开始命令)

```bash
# 创建新项目结构
mkdir -p HSTAR-Modern/src/{core,mesh,elements,materials,solver,analysis,io}
mkdir -p HSTAR-Modern/tests/{unit,integration,benchmarks}

# 初始化CMake
cd HSTAR-Modern
cmake -B build -S .

# 运行测试
cmake --build build
ctest --test-dir build
```

## FAQ (常见问题)

**Q: 应该从哪里开始重构?**
A: 从Phase 1的core模块开始，它们没有依赖，可以独立开发和测试。

**Q: 如何保持旧代码可用?**
A: 使用Adapter模式，新代码可以调用旧代码，旧代码也可以调用新代码。

**Q: 如何验证重构正确性?**
A: 每个阶段都要编写单元测试，并与原代码的数值结果进行基准对比。

**Q: 重构需要多长时间?**
A: 完整重构预计需要3-4个月。可以分阶段进行，每个阶段都产出可用的代码。

**Q: 可以只重构部分代码吗?**
A: 可以。建议优先重构材料模块（最常需要扩展）和I/O模块（最需要灵活性）。
