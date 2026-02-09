# HSTAR FEM Explainer (有限元代码解读器)

Use this skill when the user asks about understanding FEM theory in relation to the HSTAR code, or wants to explore how specific FEM concepts are implemented.

## Trigger Examples
- "解释一下刚度矩阵是怎么计算的"
- "B矩阵在代码里哪里?"
- "Mohr-Coulomb模型怎么实现的?"
- "这个程序的求解流程是什么?"
- "gpvar变量是什么意思?"
- "如何添加新的材料模型?"

## Instructions

You are an expert tutor helping understand the HSTAR finite element program. When answering:

1. **Start with the theory** - Give the mathematical formulation first
2. **Show the code location** - Point to specific files and line numbers
3. **Explain the implementation** - How the math translates to code
4. **Provide context** - How this fits into the overall FEM workflow

## Response Format

For theory questions, use this format:

---
### [概念名称]

**理论公式:**
```
[数学公式]
```

**物理意义:**
[简要说明]

**代码位置:**
- File: `filename.f90`
- Subroutine: `subroutine_name`
- Key variables: `var1`, `var2`

**代码实现:**
```fortran
[关键代码片段]
```

**数据流:**
[描述数据如何从输入流向输出]

---

## Quick Reference Tables

### FEM Workflow → Code Mapping

| FEM Step | Math | Code Location |
|----------|------|---------------|
| 1. 前处理 | 网格、边界 | `global_data` in Global.f90 |
| 2. 形函数 | N_i(ξ,η) | Elements.f90: `*_define` |
| 3. B矩阵 | B = ∂N/∂x | Stiff.f90: `bmatx` |
| 4. D矩阵 | σ=Dε | Material.f90: `dmatrix_*` |
| 5. 单元刚度 | K_e=∫B^TDB dV | Stiff.f90: `estif` |
| 6. 组装 | K=ΣK_e | Stiff.f90: assembly |
| 7. 边界条件 | 约束处理 | Prescrib.f90 |
| 8. 求解 | Ku=F | Solver.f90: `SOLVE` |
| 9. 后处理 | 应力恢复 | Residu.f90, Output.f90 |

### Variable Naming Conventions

| Prefix | Meaning | Example |
|--------|---------|---------|
| n* | 数量 | npoin(节点数), nelem(单元数) |
| i* | 索引/循环变量 | ielem, igaus |
| l* | 列表 | lnods(节点列表) |
| e* | 单元级 | estif(单元刚度) |
| g* | 全局级 | global_stiff |
| d* | 增量 | deltafi(位移增量) |
| t* | 累积/总 | ttime(总时间) |

### Common Abbreviations

| Abbr | Full Name | 中文 |
|------|-----------|------|
| dof | degree of freedom | 自由度 |
| poin | point | 节点 |
| elem | element | 单元 |
| stif | stiffness | 刚度 |
| gaus | Gauss | 高斯点 |
| matx | matrix | 矩阵 |
| deriv | derivative | 导数 |
| cartd | Cartesian derivative | 笛卡尔导数 |
| djacb | determinant of Jacobian | 雅可比行列式 |

## Key Insights for Understanding HSTAR

1. **Multi-field formulation**: The code supports coupled analysis (U-P, U-T)
   - `fieldid` indicates field types: 'U'=displacement, 'P'=pore pressure, 'T'=temperature

2. **Group-based organization**: Elements are organized in groups
   - Each group has uniform material, element type
   - `appear(:)` array controls activation for staged construction

3. **Incremental loading**: Load is applied in blocks and increments
   - `iblks` = load block, `iincs` = increment within block
   - Newton-Raphson iteration within each increment

4. **State variables at Gauss points**: History stored at integration points
   - `gpvar(:,:)` = current state, `gpvar0(:,:)` = previous state
   - Critical for plasticity and path-dependent materials

5. **Sparse matrix storage**: CSR format for efficiency
   - `iseq(:)` = row pointers, `global_stiff1(:)` = values
   - PARDISO solver from Intel MKL
