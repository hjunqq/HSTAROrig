# HSTAR Finite Element Analysis Program

## Overview

HSTAR is a comprehensive Fortran-based finite element analysis program developed by M. Pastor, Tonchun Li, and P. Mira (1997). It supports structural, thermal, and coupled analyses with advanced material models.

## Project Structure

```
HSTAR/
├── Fem.f90          # Main program entry point
├── Global.f90       # Global variables and data structures
├── Elements.f90     # Element library (shape functions, integration)
├── Material.f90     # Constitutive models
├── Stiff.f90        # Stiffness matrix calculation
├── Residu.f90       # Internal force (residual) calculation
├── Solver.f90       # Linear solvers (PARDISO, PCG, etc.)
├── Prescrib.f90     # Boundary conditions
├── Load.f90         # External loads
├── Temper.f90       # Temperature analysis
├── Output.f90       # Post-processing and output
├── Level.f90        # Level set method
├── Vartype.f90      # Variable type definitions
├── Array.f90        # Array utilities
├── gidpost.F90      # GiD post-processing interface
└── OPT.F90          # Optimization routines
```

## Key Capabilities

- **Element Types**: Line (L2,L3), Triangle (T3,T6), Quad (Q4,Q8), Tet (H4,H10), Hex (B8,B20), Beam, Contact
- **Material Models**: Elastic, Mohr-Coulomb, Drucker-Prager, Duncan-Chang, Goodman joint, Creep
- **Analysis Types**: Static, Dynamic (Newmark), Thermal, Seepage, Coupled U-P
- **Solvers**: PARDISO (Intel MKL), JPCG, Profile, PBCG
- **Features**: Staged construction, Contact analysis, Level set

## Design Documents (设计文档)

For refactoring this program, refer to the architecture design documents:

- `docs/design/architecture-design.md` - 架构设计文档 Part 1
  - 设计目标与原则
  - 系统分层架构
  - 核心数据类型定义
  - 模块规范与接口定义
  - 数据流设计

- `docs/design/architecture-design-part2.md` - 架构设计文档 Part 2
  - 输入输出格式规范
  - 错误处理策略
  - 扩展机制
  - 文件组织结构
  - 命名规范

## Build Requirements

- Intel Fortran Compiler (ifort)
- Intel MKL (for PARDISO solver)
- Visual Studio project files included

## Usage Notes

- Input file: `inp` (contains problem name and control parameters)
- Main control: `probn.msh` (mesh), `probn.mtr` (materials), `probn.lod` (loads)
- Output: GiD format supported via gidpost library

## Skills Available

The following Claude skills are available for understanding this codebase:

### Understanding the Code (理解代码)
1. **fem-theory.md** - FEM theory to code mapping (有限元原理-代码映射)
2. **code-navigator.md** - Code navigation and data structures (代码导航)
3. **fem-explainer.md** - Detailed explanations of FEM concepts (概念解读)

### Refactoring Guide (重构指南)
4. **refactor-overview.md** - Quick overview of refactoring plan (重构总览)
5. **refactor-1-architecture.md** - Target architecture and design principles (目标架构)
6. **refactor-2-modules.md** - Core module design patterns (核心模块设计)
7. **refactor-3-io-testing.md** - I/O system and testing framework (输入输出与测试)
8. **refactor-4-migration.md** - Migration strategy and steps (迁移策略)
9. **refactor-5-best-practices.md** - Modern Fortran best practices (最佳实践)

## Common Tasks

When working with this code:
- Use `Grep` to find subroutines: `subroutine NAME`
- Key entry point is `PROGRAM FEM90` in Fem.f90
- Variables are typically in Global.f90 module `global_var`
- Material models in Material.f90 follow pattern `dmatrix_*`
