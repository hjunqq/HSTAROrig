"""Validation rules for HSTAR YAML input.

Returns (errors, warnings) where errors are fatal and warnings are informational.
"""

from __future__ import annotations

from .schema import (
    HstarProject, MaterialType, AnalysisType, DofType, ElementType, Phase,
)
from .mesh_reader import MeshData


def validate_project(
    project: HstarProject,
    mesh: MeshData | None = None,
) -> tuple[list[str], list[str]]:
    """Validate the project configuration.

    Returns:
        (errors, warnings) — errors block generation, warnings are informational.
    """
    errors: list[str] = []
    warnings: list[str] = []

    mat_ids = {m.id for m in project.materials}

    # ── Material references ────────────────────────────────────────────
    for eg in project.element_groups:
        if eg.material not in mat_ids:
            errors.append(
                f"Element group {eg.group} references material {eg.material} "
                f"which does not exist. / 单元组 {eg.group} 引用了不存在的材料 {eg.material}"
            )

    # ── Node range checks ──────────────────────────────────────────────
    if mesh:
        node_ids = set(mesh.nodes.keys())
        for bc in project.boundary_conditions:
            for n in bc.nodes:
                if n not in node_ids:
                    errors.append(
                        f"Boundary condition '{bc.name}': node {n} not in mesh. "
                        f"/ 边界条件 '{bc.name}' 中节点 {n} 不在网格中"
                    )
        if project.loads:
            for pl in project.loads.point_loads:
                for n in pl.nodes:
                    if n not in node_ids:
                        errors.append(
                            f"Point load on node {n} not in mesh. "
                            f"/ 集中荷载节点 {n} 不在网格中"
                        )

    # ── DOF / analysis type compatibility ──────────────────────────────
    for eg in project.element_groups:
        if eg.dof_type == DofType.W and project.analysis.type == AnalysisType.STATIC:
            # Pure pore pressure DOF in static analysis — check if there's coupling
            pass  # Allow, as seepage analysis uses STATIC + W
        if eg.dof_type in (DofType.UW, DofType.UPC):
            has_fluid = any(
                m.phase == Phase.FLUID for m in project.materials
            )
            if not has_fluid:
                warnings.append(
                    f"Element group {eg.group} uses coupled DOF {eg.dof_type.value} "
                    f"but no FLUID material defined. / 单元组 {eg.group} 使用耦合自由度但无流体材料"
                )

    # ── Dimension / element type compatibility ─────────────────────────
    _3d_types = {ElementType.H4, ElementType.H10, ElementType.B8, ElementType.B20}
    _2d_types = {ElementType.T3, ElementType.T6, ElementType.Q4, ElementType.Q8}
    for eg in project.element_groups:
        if project.mesh.dimension == 2 and eg.element_type in _3d_types:
            errors.append(
                f"Element group {eg.group} uses 3D element {eg.element_type.value} "
                f"in 2D problem. / 单元组 {eg.group} 在2D问题中使用了3D单元"
            )
        if project.mesh.dimension == 3 and eg.element_type in _2d_types:
            warnings.append(
                f"Element group {eg.group} uses 2D element {eg.element_type.value} "
                f"in 3D problem. / 单元组 {eg.group} 在3D问题中使用了2D单元"
            )

    # ── Required material properties per type ──────────────────────────
    for m in project.materials:
        p = m.properties
        if m.type == MaterialType.ELASTIC_ISOTROPIC:
            if p.E <= 0:
                errors.append(
                    f"Material {m.id} ({m.name}): E must be > 0 for ELASTIC. "
                    f"/ 材料 {m.id}: 弹性模量 E 必须大于 0"
                )
        elif m.type == MaterialType.MOHR_COULOMB:
            if p.E <= 0:
                errors.append(
                    f"Material {m.id} ({m.name}): E must be > 0 for MOHR_COULOMB. "
                    f"/ 材料 {m.id}: 弹性模量 E 必须大于 0"
                )
            if p.friction_angle < 0 or p.friction_angle > 60:
                warnings.append(
                    f"Material {m.id}: friction_angle={p.friction_angle} "
                    f"outside typical range [0, 60]. / 内摩擦角超出常规范围"
                )
        elif m.type == MaterialType.GOODMAN:
            if p.kn <= 0:
                errors.append(
                    f"Material {m.id} ({m.name}): kn must be > 0 for GOODMAN joint. "
                    f"/ 材料 {m.id}: 法向刚度 kn 必须大于 0"
                )
        elif m.type == MaterialType.FLUID:
            if p.bulk_modulus <= 0:
                errors.append(
                    f"Material {m.id} ({m.name}): bulk_modulus must be > 0 for FLUID. "
                    f"/ 材料 {m.id}: 体积模量必须大于 0"
                )

    # ── Poisson's ratio warning ────────────────────────────────────────
    for m in project.materials:
        if m.type not in (MaterialType.FLUID, MaterialType.GOODMAN):
            if m.properties.nu > 0.49:
                warnings.append(
                    f"Material {m.id}: nu={m.properties.nu} > 0.49, "
                    f"near incompressible — may cause convergence issues. "
                    f"/ 泊松比接近不可压缩，可能导致收敛困难"
                )

    # ── Construction stage validation ──────────────────────────────────
    if project.construction_stages:
        total_steps = project.analysis.increments.total_steps
        prev_end = 0
        for cs in project.construction_stages:
            if len(cs.steps) != 2:
                errors.append(
                    f"Construction stage {cs.stage}: steps must be [start, end]. "
                    f"/ 施工阶段 {cs.stage}: steps 必须是 [起始, 终止]"
                )
                continue
            start, end = cs.steps
            if start > end:
                errors.append(
                    f"Construction stage {cs.stage}: start ({start}) > end ({end}). "
                    f"/ 施工阶段 {cs.stage}: 起始步 > 终止步"
                )
            if end > total_steps:
                errors.append(
                    f"Construction stage {cs.stage}: end ({end}) > total_steps ({total_steps}). "
                    f"/ 施工阶段 {cs.stage}: 终止步超出总步数"
                )
            if start <= prev_end and prev_end > 0:
                warnings.append(
                    f"Construction stage {cs.stage}: step range [{start}, {end}] "
                    f"overlaps with previous stage. / 施工阶段步数范围重叠"
                )
            prev_end = end

            for mc in cs.material_change:
                if mc.new_material not in mat_ids:
                    errors.append(
                        f"Construction stage {cs.stage}: material_change references "
                        f"material {mc.new_material} which does not exist. "
                        f"/ 施工阶段材料替换引用了不存在的材料"
                    )

    # ── Time curve references ──────────────────────────────────────────
    tc_ids = {tc.id for tc in project.time_curves}
    for bc in project.boundary_conditions:
        if bc.time_curve is not None and bc.time_curve not in tc_ids:
            errors.append(
                f"Boundary condition '{bc.name}': time_curve {bc.time_curve} not defined. "
                f"/ 边界条件 '{bc.name}' 引用了不存在的时间曲线"
            )
    if project.loads:
        for pl in project.loads.point_loads:
            if pl.time_curve is not None and pl.time_curve not in tc_ids:
                errors.append(
                    f"Point load: time_curve {pl.time_curve} not defined. "
                    f"/ 集中荷载引用了不存在的时间曲线"
                )
        for dl in project.loads.distributed_loads:
            if dl.time_curve is not None and dl.time_curve not in tc_ids:
                errors.append(
                    f"Distributed load: time_curve {dl.time_curve} not defined. "
                    f"/ 分布荷载引用了不存在的时间曲线"
                )

    # ── Unused materials warning ───────────────────────────────────────
    used_mats = {eg.material for eg in project.element_groups}
    for m in project.materials:
        if m.id not in used_mats:
            warnings.append(
                f"Material {m.id} ({m.name}) is defined but not used by any element group. "
                f"/ 材料 {m.id} 已定义但未被任何单元组使用"
            )

    # ── Contact stiffness warning ──────────────────────────────────────
    if project.contact:
        for gg in project.contact.gap_groups:
            if gg.normal_stiffness > 1.0e12:
                warnings.append(
                    f"Contact group {gg.id}: normal_stiffness={gg.normal_stiffness:.2e} "
                    f"is very large, may cause convergence issues. "
                    f"/ 接触法向刚度过大，可能导致收敛困难"
                )

    # ── Dynamic time step warning ──────────────────────────────────────
    if project.analysis.type == AnalysisType.DYNAMIC and project.analysis.dynamic:
        dyn = project.analysis.dynamic
        if dyn.dt > dyn.total_time / 10:
            warnings.append(
                f"Dynamic dt={dyn.dt} is > T/10={dyn.total_time/10}, "
                f"time step may be too large. / 时间步长可能过大"
            )

    return errors, warnings
