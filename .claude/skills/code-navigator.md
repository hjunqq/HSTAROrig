# HSTAR Code Navigator (代码导航助手)

Use this skill when the user wants to navigate, understand, or debug specific parts of the HSTAR finite element code.

## Instructions

Help users navigate the HSTAR codebase by:
1. Finding where specific functionality is implemented
2. Explaining variable meanings and data structures
3. Tracing execution flow for specific operations
4. Identifying connections between modules

## Key Data Structures

### Global Variables (Global.f90)

```fortran
! Mesh data
npoin          ! 节点总数
nelem          ! 单元总数
ndimn          ! 问题维度 (2D/3D)
ngroup         ! 单元组数
coord(:,:)     ! 节点坐标 coord(ndimn, npoin)

! Degrees of freedom
ntotv          ! 总自由度数
mdofn          ! 最大节点自由度
nodfn(:,:)     ! 节点自由度编号 nodfn(mdofn, npoin)

! Solution vectors
deltafi(:)     ! 位移增量向量
result_zero(:) ! 累积位移
stfor(:)       ! 内力向量
tofor(:)       ! 总力向量 (外力-内力)

! Time stepping
iblks          ! 当前荷载块
iincs          ! 当前增量步
iiter          ! 当前迭代次数
ditime         ! 时间步长
ttime          ! 累积时间
```

### Element Library (Elements.f90)

```fortran
type element_lib
    name           ! 单元名称
    index          ! 单元类型索引
    group          ! 所属组号
    matno          ! 材料号
    nstre          ! 应力分量数
    estif(:,:)     ! 单元刚度矩阵
    field(:)       ! 场信息 (位移场、压力场等)
    egaus(:)       ! 高斯点信息
end type

type element_field
    lnods_f(:)     ! 单元节点编号
    ldofs_f(:)     ! 单元自由度编号
    gpvar(:,:)     ! 高斯点变量 (应力、应变等)
    gpvar0(:,:)    ! 上一步高斯点变量
    state(:)       ! 材料状态 (弹性/塑性)
end type
```

### Group Structure (Global.f90)

```fortran
type group_of_elements
    name           ! 组名
    index          ! 单元类型索引
    matno          ! 材料号
    nelgroup       ! 组内单元数
    list(:)        ! 单元列表
    nstre          ! 应力分量数
    fieldid        ! 场标识 ('U', 'P', 'T'...)
end type
group(:)           ! 单元组数组
appear(:)          ! 组激活状态
```

### Material Properties (Material.f90)

```fortran
type property_solid
    material       ! 材料类型名称
    e              ! 弹性模量
    nu             ! 泊松比
    density        ! 密度
    c              ! 粘聚力
    phi            ! 摩擦角
    icreep         ! 蠕变类型
end type
props(:)           ! 材料属性数组
```

## Program Flow

### Main Program (Fem.f90)

```
FEM90 (主程序)
├── global_data          ! 读取网格、材料数据
├── material_set         ! 设置材料参数
├── modf_element_lib     ! 修改单元库
├── contact_point_to_point  ! 接触对设置
│
└── DO iblks = lblks+1, runblks    ! 荷载块循环
    ├── external_load_control       ! 读取荷载控制
    │
    └── DO iincs = 1, nincs         ! 增量步循环
        ├── external_load           ! 计算外力向量
        │
        └── DO iiter = 1, miter     ! Newton迭代
            ├── stiff_u             ! 计算刚度矩阵
            ├── stiff_assemble      ! 组装总刚度
            ├── solve               ! 求解线性方程
            ├── residu_f            ! 计算残差
            │
            └── 收敛检查
                ├── 收敛 → 下一增量步
                └── 不收敛 → 继续迭代
```

### Stiffness Calculation (Stiff.f90)

```
STIFF_U
├── DO igroup = 1, ngroup
│   └── DO ielgroup = 1, nelgroup
│       ├── 获取单元信息
│       ├── DO igaus = 1, ngaus
│       │   ├── 计算形函数导数 cartd
│       │   ├── 计算B矩阵 bmatx
│       │   ├── 计算D矩阵 dmatx (调用材料子程序)
│       │   ├── 计算 dbmat = D * B
│       │   └── estif += w * |J| * B^T * D * B
│       │
│       └── 存储单元刚度 element(ielem)%estif
```

### Solver Options (Solver.f90)

```
SOLVE
├── PARDISO    ! Intel MKL直接求解器 (推荐)
├── JPCG       ! Jacobi预条件共轭梯度
├── PROFILE    ! Profile法直接求解
├── PBCG       ! 预条件双共轭梯度
└── SSORPBCG   ! SSOR预条件BCG
```

## Common Debugging Points

### Check Stiffness Matrix
- `Stiff.f90:52` - 单元循环开始
- `element(ielem)%estif` - 单元刚度
- `global_stiff1(:)` - 总刚度(CSR格式)

### Check Internal Force
- `Residu.f90` - RESIDU_F子程序
- `stfor(:)` - 内力向量
- `gpvar(:,:)` - 高斯点应力

### Check Convergence
- `Fem.f90` - 主迭代循环
- `toler_force` - 力容差
- `err` - 当前误差

### Check Material Response
- `Material.f90` - 各本构模型子程序
- `dmatx(:,:)` - 当前切线刚度
- `state(:)` - 材料状态

## File Quick Reference

| 功能 | 文件 | 关键子程序 |
|------|------|-----------|
| 主程序 | Fem.f90 | PROGRAM FEM90 |
| 全局数据 | Global.f90 | global_data |
| 单元库 | Elements.f90 | kinddefine, *_define |
| 材料 | Material.f90 | material_set, dmatrix_* |
| 刚度 | Stiff.f90 | STIFF_U |
| 内力 | Residu.f90 | RESIDU_F |
| 求解 | Solver.f90 | SOLVE, MAIN_PARDISO |
| 边界 | Prescrib.f90 | prescrib_set |
| 荷载 | Load.f90 | external_load |
| 输出 | Output.f90 | out_full_write, OUT_GID_WRITE |
| 温度 | Temper.f90 | temperature analysis |
| Level Set | Level.f90 | level set method |
