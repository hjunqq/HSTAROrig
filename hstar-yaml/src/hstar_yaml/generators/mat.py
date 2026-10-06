"""Generate the .mat (material) file.

Format (from Material.f90:material_set()):
    Line 1: text (title)
    Line 2: nscurve (number of stress-strain curves, usually 0)
    Line 3: nline (number of description lines)
    Lines 4..4+nline: description text
    Line: text
    Line: mmats (total material count)
    For each material:
        Line: material_serial  id
        Line: MECHANICAL/THERMAL  SOLID/FLUID/CONTACT  phase_count
        Line: number_of_phases
        Line: SOLID/FLUID
        Line: TYPE  density  ratio  thickness  E  nu  thermal_exp  ...
        Line: extra params (hardening, curves, etc.)
"""

from __future__ import annotations

import math

from .base import BaseGenerator
from ..schema import MaterialType, Phase


class MatGenerator(BaseGenerator):
    extension = ".mat"

    def build(self) -> str:
        p = self.project
        ndimn = p.mesh.dimension
        lines = []

        # Header
        lines.append(f" *({p.problem_name}.mat)输入材料信息material Property")

        # Stress-strain curves (0 for now)
        lines.append("  0")

        # Description lines
        desc = [
            "MATERIAL PROPERTIES-INPUT FOR 1 TO NMATS",
            "NAME FOR MATERIAL---FOR COMMENTS",
            "FIRST:NAME---MECHANICAL, TEMPERATURE, ETC",
            "FOR MECHNICAL, SECOND:PHASE---SOLID, FLUID, AIR OR OIL",
            "NUMBER OF PHASE",
            "FOR SOLID: MATERIAL(CLASSICALEP,CAMCLAY,ETC),density,ratio,thick,e,nu",
            "FOR MATERISL: CLASSICALEP:CRITERIA(TC,VM,DP,MC),sigm0,hardening",
            "FOR MC OR DP: FRICT_ANGLE, DILAN_ANGLE",
            "CAMCLAY    :Pc, lamda, Mg, Mf, D0, D1,gaama",
            "FOR FLUID: (1) density,ratio,bulkw (2)permiability(1:ndimn)",
        ]
        lines.append(f"   {len(desc)}")
        for d in desc:
            lines.append(d)

        # Total materials
        lines.append(" 总的材料参数输入次数")
        lines.append(f"  {len(p.materials)}")

        for mat in p.materials:
            props = mat.properties
            lines.append(f"     material_serial         {mat.id}")

            if mat.type == MaterialType.FLUID:
                # Fluid material
                lines.append(f"          MECHANICAL               WATER    {mat.id}")
                lines.append("  1")
                lines.append("               FLUID")
                # density, ratio, bulk_modulus, wave_speed
                wave_speed = 0
                if props.bulk_modulus > 0 and props.density > 0:
                    wave_speed = int(math.sqrt(props.bulk_modulus / props.density))
                lines.append(
                    f"   {self._fmtg(props.density)},{self._fmtg(1.0)},{self._fmtg(props.bulk_modulus)} {wave_speed}"
                )
                # permeability
                perm = props.permeability if props.permeability else [1.0e-3] * ndimn
                lines.append("   " + "      ".join(self._fmtg(v) for v in perm[:ndimn]))
            elif mat.type == MaterialType.GOODMAN:
                # Goodman joint — CONTACT phase
                lines.append(
                    f"          MECHANICAL               CONTACT   {mat.id}"
                )
                lines.append("  1")
                lines.append("               SOLID")

                # Property line with extra trailing zeros
                lines.append(
                    f"GOODMAN                   {self._fmtg(props.density)}      "
                    f"{self._fmtg(1.0)}      "
                    f"{self._fmtg(props.thickness)}      "
                    f"{self._fmtg(props.E)}      "
                    f"{self._fmtg(props.nu)}       "
                    f"{self._fmtg(props.thermal_expansion)}    0    0    0    0    0   "
                )

                # Hardening line with extra 0
                lines.append(
                    f"  0  0      {self._fmtg(props.hardening)} 0"
                )

                # Goodman params: 99, friction, kn, 1
                lines.append(
                    f"    99      {self._fmtg(props.friction)}      {self._fmtg(props.kn)}   1"
                )

                # JANBU line if applicable
                if props.janbu_modulus > 0:
                    lines.append(
                        f"      JANBU {props.janbu_exponent} {props.janbu_k0}   "
                        f"{self._fmtg(props.janbu_modulus)} {props.janbu_phi}"
                    )

                # Extended Goodman params
                lines.append(
                    f"{self._fmtg(props.kn)} {self._fmtg(props.ks)} "
                    f"{self._fmtg(props.E)} 0.99 {props.friction_angle} "
                    f"{self._fmtg(props.kn)} {self._fmtg(props.tensile_strength)} "
                    f"{self._fmtg(props.ks)} {self._fmtg(props.cohesion)}  {self._fmtg(props.kn)}"
                )

                # Trailing kn value (extra line)
                lines.append(f"{self._fmtg(props.kn)}")
            else:
                # Other solid materials
                phase_str = "SOLID"

                lines.append(
                    f"          MECHANICAL               {phase_str}    {mat.id}"
                )
                lines.append("  1")
                lines.append("               SOLID")

                # Material type name for HSTAR
                hstar_type = self._get_hstar_material_name(mat.type)

                # Common line: TYPE density ratio thickness E nu thermal_exp 0 0 0
                lines.append(
                    f"{hstar_type}         "
                    f"{self._fmtg(props.density)}      "
                    f"{self._fmtg(1.0)}      "
                    f"{self._fmtg(props.thickness)}      "
                    f"{self._fmtg(props.E)}     "
                    f"{self._fmtg(props.nu)}      "
                    f"{self._fmtg(props.thermal_expansion)}    0    0    0"
                )

                # Extra params depending on type
                if mat.type == MaterialType.ELASTIC_ISOTROPIC:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)}"
                    )
                elif mat.type == MaterialType.MOHR_COULOMB:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)} 0"
                    )
                    lines.append(
                        f"  MC  {self._fmtg(props.cohesion)}  {props.friction_angle}  {props.dilation_angle}"
                    )
                elif mat.type == MaterialType.DRUCKER_PRAGER:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)} 0"
                    )
                    lines.append(
                        f"  DP  {self._fmtg(props.cohesion)}  {props.friction_angle}  {props.dilation_angle}"
                    )
                elif mat.type == MaterialType.VON_MISES:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)} 0"
                    )
                    lines.append(
                        f"  VM  {self._fmtg(props.cohesion)}  0  0"
                    )
                elif mat.type == MaterialType.DUNCAN_CHANG:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)} 0"
                    )
                    lines.append(
                        f"  DUNCANCHANG {self._fmtg(props.cohesion)} {props.friction_angle} "
                        f"{self._fmtg(props.Rf)} {self._fmtg(props.K_dc)} "
                        f"{props.n_dc} {self._fmtg(props.pa)}"
                    )
                elif mat.type == MaterialType.CAMCLAY:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)} 0"
                    )
                    lines.append(
                        f"  CAMCLAY {self._fmtg(props.Pc)} {self._fmtg(props.lamda)} "
                        f"{self._fmtg(props.Mg)} {self._fmtg(props.Mf)} "
                        f"{self._fmtg(props.D0)} {self._fmtg(props.D1)} {self._fmtg(props.kappa)}"
                    )
                else:
                    lines.append(
                        f"  0  0      {self._fmtg(props.hardening)} 0"
                    )

        return "\n".join(lines) + "\n"

    @staticmethod
    def _get_hstar_material_name(mtype: MaterialType) -> str:
        """Map material type enum to HSTAR material name string."""
        mapping = {
            MaterialType.ELASTIC_ISOTROPIC: "ELASTIC_ISOTROPIC",
            MaterialType.MOHR_COULOMB: "CLASSICALEP",
            MaterialType.DRUCKER_PRAGER: "CLASSICALEP",
            MaterialType.VON_MISES: "CLASSICALEP",
            MaterialType.TRESCA: "CLASSICALEP",
            MaterialType.GOODMAN: "GOODMAN",
            MaterialType.DUNCAN_CHANG: "DUNCANCHANG",
            MaterialType.CAMCLAY: "CAMCLAY",
            MaterialType.CONCRETE: "CONCRETE",
            MaterialType.CREEP: "CREEP",
        }
        return mapping.get(mtype, "ELASTIC_ISOTROPIC")
