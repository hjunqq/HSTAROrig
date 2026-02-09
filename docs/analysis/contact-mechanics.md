# HSTAR 接触力学系统分析

## 一、接触数据结构

### 1.1 gap_group 类型 (Global.f90:429-442)

接触对组的主要数据结构:

```fortran
type gap_group
    ! 控制标志
    integer(ink) npairs         ! 接触对数量
    integer(ink) ngroupt        ! 顶面组数
    integer(ink) stateix        ! 接触状态索引
    integer(ink) frict_less     ! 无摩擦标志
    integer(ink) goodman        ! Goodman模型标志
    integer(ink) thin_layer     ! 薄层标志
    integer(ink) ivcoh, ivfri   ! 内聚力/摩擦标志

    ! 节点与配对
    integer(ink),pointer:: pairnode(:,:)  ! 接触节点对 (2*nnode, npairs)
    integer(ink),pointer:: state(:)       ! 当前接触状态
    integer(ink),pointer:: state0(:)      ! 上步接触状态
    integer(ink),pointer:: statei(:)      ! 迭代中间状态
    integer(ink),pointer:: paire(:)       ! 所属单元
    integer(ink),pointer:: group(:)       ! 所属组
    integer(ink),pointer:: pair_process(:)! 处理标志

    ! 几何
    real(irk),pointer:: rot(:,:,:)     ! 局部旋转矩阵 (3,3,npairs)
    real(irk),pointer:: aera(:)        ! 接触面积
    real(irk),pointer:: gap0(:,:)      ! 初始间隙 (ndimn,npairs)
    real(irk),pointer:: gap(:,:)       ! 当前间隙

    ! 力和位移
    real(irk),pointer:: ctforce(:,:)   ! 当前接触力 (ndimn,npairs)
    real(irk),pointer:: ctforce0(:,:)  ! 上步接触力
    real(irk),pointer:: ctforcei(:,:)  ! 迭代接触力i
    real(irk),pointer:: ctforcej(:,:)  ! 迭代接触力j
    real(irk),pointer:: dxyz(:,:)      ! 位移增量 (ndimn,npairs)
    real(irk),pointer:: dxyz0(:,:)     ! 上步位移增量
    real(irk),pointer:: dxyzi(:,:)     ! 迭代位移增量

    ! 材料参数
    real(irk),pointer:: ft(:)          ! 拉伸极限
    real(irk),pointer:: cohes(:)       ! 内聚力
    real(irk),pointer:: frict(:)       ! 摩擦角
    real(irk),pointer:: ft0(:)         ! 初始拉伸极限
    real(irk),pointer:: cohes0(:)      ! 初始内聚力
    real(irk),pointer:: frict0(:)      ! 初始摩擦角

    ! 刚度
    real(irk),pointer:: kgroup0(:,:)   ! 初始组刚度 (3*(ndimn-1), 3*(ndimn-1))
    real(irk),pointer:: kgroup1(:,:)   ! 当前组刚度
    real(irk),pointer:: kxyz0(:,:,:)   ! 局部刚度张量 (ndimn,ndimn,npairs)
    real(irk),pointer:: kxyz(:,:,:)    ! 当前局部刚度

    ! 损伤
    real(irk),pointer:: damage0(:)     ! 上步损伤变量
    real(irk),pointer:: damage(:)      ! 当前损伤变量
    real(irk),pointer:: sigmad(:)      ! 损伤应力
    real(irk),pointer:: Gf0(:), Gf(:) ! 断裂能

    ! 标量参数
    real(irk) gapi                     ! 初始间隙
    real(irk) Rf, n, Pa, e            ! Janbu模型参数
    real(irk) miu                      ! 摩擦系数
    real(irk) thick                    ! 界面厚度
    real(irk) Ke                       ! 刚度系数
    real(irk) gamaw                    ! 水的重度
end type gap_group
```

### 1.2 gap_block_group 类型 (Global.f90:444-457)

刚体块接触结构:

```fortran
type gap_block_group
    integer(ink) npgblock        ! 刚体块潜在接触点数
    integer(ink) ngroupb         ! 底面组数 (可变形基础)
    integer(ink) ngroupt         ! 顶面组数 (刚体块表面)
    integer(ink) ntotv_bt        ! 块接触总自由度
    integer(ink) nrdof           ! 刚体自由度数 (0=可变形, >0=刚体)
    integer(ink) npblock         ! 块节点数
    integer(ink) eblock          ! 外部块标志 (2017)

    integer(ink),pointer:: nodeblock(:)       ! 块节点编号
    integer(ink),pointer:: nppt(:)            ! 节点-点映射
    integer(ink),pointer:: nodegblock(:)      ! 块组节点
    integer(ink),pointer:: nodegblock_ipairs(:) ! →接触对索引
    integer(ink),pointer:: nodegblock_igaps(:)  ! →间隙组索引
    integer(ink),pointer:: nodegblock_onetwo(:) ! →表面1/2标识
    integer(ink),pointer:: ldofs(:)           ! 局部自由度列表
    integer(ink),pointer:: listrele(:)        ! 释放/刚体自由度

    real(irk),pointer:: cmatrix(:,:)    ! 柔度矩阵 (ntotv_bt, ntotv_bt)
    real(irk),pointer:: rdisp_zero(:)   ! 刚体参考位移
    real(irk),pointer:: npdisp(:,:,:)   ! 节点位移矩阵 (ndimn, npblock, nrdof)
    real(irk),pointer:: rdisp_inc(:)    ! 刚体位移增量
    real(irk),pointer:: center(:)       ! 质心坐标
    real(irk),pointer:: mass_inertia(:) ! 质量和惯量 (3*(ndimn-1))
    real(irk),pointer:: rstiff(:,:)     ! 刚体刚度
    real(irk),pointer:: rdisp_first(:)  ! Newmark速度
    real(irk),pointer:: rdisp_second(:) ! Newmark加速度
end type gap_block_group
```

### 1.3 全局控制变量 (Global.f90:46, 156-162)

```fortran
integer(ink) ngaps        ! 接触间隙组总数
integer(ink) ngapb        ! 刚体块接触组总数
integer(ink) npbt         ! 块接触节点总数
integer(ink) ntotvbt      ! 块接触总自由度
integer(ink) contactpe    ! 接触点评估方式: 1=单点, 2=多点
integer(ink) nonsbt       ! 非对称块刚度标志
integer(ink) xlwsol       ! 扩展摩擦模型标志
integer(ink) mpairs       ! 最大接触对数
integer(ink) miter_state  ! 接触状态最大迭代次数
integer(ink) block_stab   ! 块稳定方法: 0=无, 1=完全刚体, 2=点对
```

---

## 二、接触状态机

### 2.1 主入口

**位置:** `Fem.f90:12711-12997` (`contact_state` 子程序)

调用方式:
- `call contact_state(0)` — 初始化 (每步首次迭代)
- `call contact_state(1)` — 更新 (迭代结束后)

### 2.2 状态定义

| 状态值 | 含义 |
|--------|------|
| 0 | 不活跃/未初始化 |
| 1 | 闭合/受压接触 |
| 2-3 | 摩擦滑移 |
| ≥4 | 特殊状态 (拉伸、损伤等) |

字符表示:
- `'open'` — 间隙张开，无接触力
- `'contact'` — 法向接触，存在法向力
- `'sliding'` — 摩擦滑移接触

### 2.3 状态转换逻辑

**代码 (Fem.f90:12920-12945):**

```
前状态='contact':
    if σ_mean > ft0 (拉应力超过极限):
        → 'open'
    else:
        → 'contact' (保持)

前状态='open':
    if gap_change < ε AND current_gap < natural_gap:
        → 'contact' (闭合)
    else:
        → 'open' (保持)
```

### 2.4 间隙初始化 (Fem.f90:12791-12824)

- **igap0=1** (均匀间隙): 所有高斯点统一赋值 gap0
- **igap0=2, gap_kind=1** (空间变化): 用形函数插值节点间隙值到高斯点
- **igap0=99** (由位移计算): 从节点相对运动计算间隙

### 2.5 收敛检查 (Fem.f90:12977-12995)

```fortran
iiii = 0  ! 已收敛单元计数
jjjj = 0  ! 总单元计数
do igroup...
    do ielgroup...
        jjjj = jjjj + 1
        if (element(ielem)%field(1)%icok == 1) iiii = iiii + 1
        ! icok=1 表示状态未改变: state1 == state
    end do
end do
if (iiii == jjjj) iccontact = 1  ! 所有单元收敛
```

---

## 三、接触力求解

### 3.1 求解流程 (Solver.f90:2171-3789)

`solve_ctt` 子程序的主要步骤:

1. **初始化** (2171-2200): 分配位移和接触力向量
2. **位移恢复** (2242-2290): 计算刚体运动产生的位移跳跃
3. **组装** (2300-2380): 构建柔度矩阵/刚度矩阵
4. **求解** (2382-2752): 使用 PROFILE 或 PARDISO 求解块接触系统
5. **力提取** (2495-2752): 从位移解提取接触力

### 3.2 接触刚度组装 (Solver.f90:2321-2355)

```fortran
do igapb = 1, ngapb
    do i0 = 1, gapb(igapb)%npgblock
        ipairs = gapb(igapb)%nodegblock_ipairs(i0)
        igaps = gapb(igapb)%nodegblock_igaps(i0)
        istate = gaps(igaps)%state(ipairs)

        if (istate==0) cycle  ! 不活跃跳过
        if (istate==1 .or. istate>=4 .or. ...) then
            ! 组装 kxyz 到全局块系统
            ! K_bt(ieq,jeq) += kxyz(idimn,jdimn,ipairs) * fact
        endif
    end do
end do
```

### 3.3 右端向量组装 (Solver.f90:2422-2476)

残差方程:
```
R_contact = -dxyz - K_contact · (F_contact - F_contact_old)
```

其中:
- `dxyz`: 闭合间隙所需的位移差
- `K_contact`: 接触刚度矩阵 (kxyz)
- `F_contact`, `F_contact_old`: 当前和上步接触力

### 3.4 力更新 (Solver.f90:2495-2752)

求解 `K_bt · u_bt = R_bt` 后:

```fortran
! 提取刚体位移
rdisp = result_bt(对应自由度)

! 更新接触力: F_new = F_old + K · u
gaps(igaps)%ctforce(idimn, ipairs) = &
    gaps(igaps)%ctforce0(idimn, ipairs) + &
    dot_product(gaps(igaps)%kxyz(:, idimn, ipairs), rdisp)
```

---

## 四、柔度矩阵 (A矩阵) 构建

### 4.1 forAdirect 子程序 (Fem.f90:6724-6984)

#### 步骤1: 质心计算 (6744-6772)

```fortran
center = Σ(djacb · gpcod · density) / Σ(djacb · density)
```

对所有底面组的高斯点积分。

#### 步骤2: 位移-刚体DOF关系 (6827-6863)

对每个块节点，构建位移矩阵:

```
u_node = npdisp · u_rigid
```

2D (3个刚体DOF: tx, ty, θz):
```
npdisp = | 1  0  -Δy |    Δx = x_node - x_center
         | 0  1   Δx |    Δy = y_node - y_center
```

3D (6个刚体DOF: tx, ty, tz, θx, θy, θz):
```
npdisp = | 1  0  0   0   Δz  -Δy |
         | 0  1  0  -Δz  0    Δx |
         | 0  0  1   Δy -Δx   0  |
```

#### 步骤3: 柔度矩阵组装 (6879-6959)

对每个接触点:

```fortran
rot = gaps(igaps)%rot(:,:,ipair)  ! 局部旋转矩阵
dist = coord(:, contact_point) - center

! 全局位移到局部接触坐标的变换
disli = rot @ disgi * coef  ! coef = ±1 (表面1/2)

! 填充柔度矩阵
C(contact_DOF, rigid_DOF) = disli
C(rigid_DOF, contact_DOF) = disli  ! 对称
```

#### 步骤4: 动力效应 (6969-6976)

```fortran
! 加入惯性项
C(rigid_DOF, rigid_DOF) -= mass_inertia * (1 + damp*θ₁*Δt)
```

#### 柔度矩阵结构

```
C = | C_pp  C_pr |    C_pp: 接触点→接触点
    | C_rp  C_rr |    C_pr: 接触点→刚体DOF
                      C_rp: 刚体DOF→接触点
                      C_rr: 刚体DOF→刚体DOF (惯性/阻尼)
```

满足: `u_contact = C · F_contact`

---

## 五、接触本构模型

### 5.1 Goodman接头模型

**材料定义:** `Material.f90:56-66` (`material_6`)

**刚度计算:** `Stiff.f90:6067-6200` (`PKPN` 子程序)

#### Janbu型应力依赖刚度

```
K_n = Kzz                    (法向刚度 — 常数)
K_s = K1·γw·(|σn|/Pa)^n·[1 - Rf·|τ|/τmax]²   (切向刚度)
```

其中:
- K1: 参考切向刚度
- γw: 水重度缩放
- Pa: 大气压参考值
- n: 应力依赖指数
- Rf: 破坏比 (0-1)
- τmax = c + σn·tan(φ) (Mohr-Coulomb包络线)

#### 开裂判断

```fortran
if (σn >= ft .or. |σn| < 0.01) then
    K_n = Pa (极小值)
    K_s = Pa
endif
```

### 5.2 FCM损伤模型 (虚拟裂缝模型)

**位置:** `Stiff.f90:6204-6290, 6311-6367`

**材料定义:** `Material.f90:42-46` (`material_fcm`)

```fortran
type material_fcm
    real(irk) w0, w1, w2   ! 裂缝张开位移参数
    real(irk) ft, ft1, Gf  ! 拉伸强度, 软化参数, 断裂能
    real(irk),pointer:: kns0(:)  ! 法向/切向初始刚度
    integer(ink) xlwmodel   ! 软化模型类型
end type material_fcm
```

#### 软化模型 (xlwmodel)

**1 — 线性软化:**
```
σ(w) = ft · (1 - w/w0)         当 w ≤ w0
σ(w) = 0.01                     当 w > w0
```

**2 — 双线性软化 (Petersson):**
```
σ(w) = ft + (ft1-ft)·w/w1      当 w ≤ w1
σ(w) = ft1·(w0-w)/(w0-w1)      当 w1 < w ≤ w0
```

**3 — 指数软化:**
```
σ(w) = ft · [(1+(w/w0)³)·exp(-5.64w/w0) - (w/w0)·7.106×10⁻³]
```

**5 — 多线性 (杜加吉等):**
```
σ = ft                           当 w ≤ w1
σ = ft + (ft1-ft)·(w-w1)/(w2-w1)  当 w1 < w ≤ w2
σ = ft1·(w0-w)/(w0-w2)           当 w2 < w < w0
σ = 0.01                          当 w ≥ w0
```

**损伤参数:**
```
D = 1 - σ(w) / (gap · kns0)
```

### 5.3 MCJOINT (Mohr-Coulomb接头)

**位置:** `Residu.f90:4018-4114`

在局部接触坐标系中:
```
f = τ + σn·tan(φ) - c ≤ 0
```

超过屈服面时:
```fortran
τ_allowed = -tan(φ)·σn + c
if (τ_allowed < 0) τ_allowed = 0
σ_shear = τ_allowed · σ_shear/|σ_shear|   ! 按比例缩减
```

### 5.4 Watertight接头 (IWJ模型)

**位置:** `Material.f90:56-66`, 17参数A数组

用于防水接头的低渗透性处理。

---

## 六、接触力→全局力的传递

### 6.1 ctfor_to_tofor 子程序

**位置:** `Solver.f90:4839-4910`

```fortran
subroutine ctfor_to_tofor(tofor0, toforx)
    toforx = tofor0   ! 复制基准力

    do igaps = 1, ngaps
        do ipairs = 1, gaps(igaps)%npairs
            ! 跳过不活跃的对
            if (gaps(igaps)%pair_process(ipairs) == 0) cycle

            ! 局部→全局旋转
            rot = gaps(igaps)%rot(:,:,ipairs)
            ctforl1 = gaps(igaps)%ctforce(:,ipairs)
            ctfori1 = transpose(rot) @ ctforl1    ! 表面1
            ctfori2 = -ctfori1                     ! 表面2 (等值反向)

            ! 分配到节点 (多点时平均)
            coef1 = 1.0/nnodej  (contactpe==2时)

            ! 表面1节点: 加正力
            toforx(nodfn(1:ijk, node1)) += ctfori1 * coef1
            ! 表面2节点: 加反力
            toforx(nodfn(1:ijk, node2)) += ctfori2 * coef1
        end do
    end do
end subroutine
```

**调用位置:** Fem.f90 中多处:
- 准静态: 2598, 3187
- 动力: 5210, 6104, 6331
- 频率: 7789
- 其它: 8832, 9085

---

## 七、块接触状态与标志

| 参数 | 值 | 含义 |
|------|-----|------|
| block_stab | 0 | 直接接触，无稳定 |
| | 1 | 完全刚体稳定 (含旋转) |
| | 2 | 点对稳定 (2019) |
| contactpe | 1 | 每对单点接触 |
| | 2 | 多点接触 (力分配到节点) |
| nonsbt | 0 | 对称组装 |
| | 1 | 非对称矩阵 |
| xlwsol | 0 | 标准Mohr-Coulomb |
| | 1 | 扩展滑移状态 |
| miter_state | >1 | 块内多次状态迭代 |
| istatec | 0 | 用初始状态(state0)计算刚度 |
| | 1 | 用当前状态计算刚度(切线) |

---

## 八、钢筋-混凝土粘结滑移

### 8.1 数据结构 (Global.f90:461-467)

```fortran
type line_steel_stick
    integer(ink) npairs_sc, nline_s
    integer(ink),pointer:: linenode_s(:,:)    ! 钢筋线节点
    integer(ink),pointer:: pairnode_sc(:)     ! 粘结对节点
    real(irk),pointer:: rot_sc(:,:)          ! 旋转矩阵
    real(irk),pointer:: aera_sc(:)           ! 接触面积
    real(irk),pointer:: tao_cs(:)            ! 粘结剪应力
    real(irk),pointer:: slip_sc(:)           ! 滑移位移
    real(irk),pointer:: cmatrix_c(:,:)       ! 混凝土柔度矩阵
    real(irk),pointer:: kmatrix_s(:,:)       ! 钢筋刚度矩阵
end type line_steel_stick
```

### 8.2 cmatrix_c_formation (Fem.f90:4392-4488)

构建混凝土侧柔度矩阵，用于钢筋-混凝土界面求解。

类似的 `Tcmatrix_c_formation` 用于冷却水管接触。
