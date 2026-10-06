"""Generate the .glb (global) file.

This is the most complex HSTAR input file with ~50 parameters in strict order.
The format is read by Global.f90:global_data() starting at line 675.

Key sections (in order):
1. Header: npoin, nelem, ndimn, nmats, ngroup, etc.
2. Control: ninit, kinit, nblks, nlinks, etc.
3. Problem type: type_problem, type_solver, type_load, type_nl, etc.
4. Physical: nmass, nsmat, nhmat, etc.
5. Temperature params
6. MDOFN (master DOF)
7. Newmark params (beeta1, beeta2, theta1)
8. Construction stages: equvs_process, appear_level, appear_process, matno_process
9. Output flags (gid_u, gid_s, ...)
10. Water/IFS params
11. Contact params (from tension_joint to ngaps)
12. Group definitions
"""

from __future__ import annotations

from .base import BaseGenerator
from ..schema import (
    AnalysisType, DofType, ElementType, NonlinearMethod, SolverType,
    MaterialType, OutputFormat, OutputMode, Phase,
)


# Map element types to HSTAR element index
ELEMENT_INDEX = {
    ElementType.L2: 1,
    ElementType.L3: 1,
    ElementType.T3: 2,
    ElementType.T6: 3,
    ElementType.Q4: 4,
    ElementType.Q8: 5,
    ElementType.H4: 6,
    ElementType.H10: 7,
    ElementType.B8: 9,
    ElementType.B20: 10,
    ElementType.BEAM: 20,
    ElementType.CONTACT: 22,
}

# Map element types to HSTAR short names (always Q4 in .glb group line)
ELEMENT_NAMES = {
    ElementType.L2: "L2",
    ElementType.L3: "L3",
    ElementType.T3: "T3",
    ElementType.T6: "T6",
    ElementType.Q4: "Q4",
    ElementType.Q8: "Q4",
    ElementType.H4: "Q4",
    ElementType.H10: "Q4",
    ElementType.B8: "Q4",
    ElementType.B20: "Q4",
    ElementType.BEAM: "BEAM",
    ElementType.CONTACT: "CT",
}

# Map DOF types to HSTAR field identifiers and DOF lists
DOF_INFO = {
    # dof_type: (fieldid, nrfields, [(nfdof, listdof), ...])
    DofType.U: ("U", 1, [(None, None)]),     # ndimn DOFs: [1,2] or [1,2,3]
    DofType.UW: ("UW", 2, [(None, None), (1, [8])]),  # U + pore pressure
    DofType.UPC: ("UP", 2, [(None, None), (1, [8])]),
    DofType.T: ("T", 1, [(1, [10])]),
    DofType.W: ("W", 1, [(1, [8])]),
}


class GlbGenerator(BaseGenerator):
    extension = ".glb"

    def build(self) -> str:
        p = self.project
        a = p.analysis
        g = p.glb_control
        ngroup = len(p.element_groups)
        nmats = len(p.materials)
        ndimn = p.mesh.dimension
        nblks = max(len(p.construction_stages), 1)

        # Count nodes/elements from mesh
        npoin = self.mesh.npoin if self.mesh else 0
        nelem = self.mesh.nelem if self.mesh else 0
        # npoinb is typically same as nelem for solid problems
        # but for problems with boundary nodes it differs
        npoinb = nelem

        # outplot
        if p.output.format == OutputFormat.GID:
            outplot = "GIDA" if p.output.mode == OutputMode.APPEND else "GIDR"
        else:
            outplot = "COSMOSR"

        # type_problem
        if a.type == AnalysisType.STATIC:
            type_problem = "F"
        elif a.type == AnalysisType.DYNAMIC:
            type_problem = "E"
        elif a.type == AnalysisType.FREQUENCY:
            type_problem = "Q"
        elif a.type == AnalysisType.THERMAL:
            type_problem = "T"
        elif a.type == AnalysisType.SEEPAGE:
            type_problem = "S"
        else:
            type_problem = "F"

        # type_solver
        type_solver = p.solver.type.value

        # type_nl
        if a.nonlinear == NonlinearMethod.LINEAR:
            type_nl = 0
        elif a.nonlinear == NonlinearMethod.NEWTON_RAPHSON:
            type_nl = 5
        elif a.nonlinear == NonlinearMethod.MODIFIED_NR:
            type_nl = 1
        else:
            type_nl = 5

        # type_load
        type_load = "LOAD"

        # Has gravity?
        ngrav = g.ngrav

        # DOF configuration: determine mdofn and lmdofn
        # 10 possible DOFs: Ux Uy Uz Thxy Thyz Thzx P W A T
        lmdofn = [0] * 10
        # Enable displacement DOFs based on dimension
        for i in range(ndimn):
            lmdofn[i] = 1
        # Check if any group uses pore pressure or temperature
        for eg in p.element_groups:
            if eg.dof_type in (DofType.UW, DofType.UPC, DofType.W):
                lmdofn[7] = 1  # W (pore pressure)
            if eg.dof_type == DofType.T:
                lmdofn[9] = 1  # Temperature
        # mdofn = actual count of active DOFs
        mdofn = sum(lmdofn)

        # order_time for each DOF (0=none, 1=speed, 2=accel)
        order_time = [0] * 10
        if a.type == AnalysisType.DYNAMIC:
            for i in range(ndimn):
                order_time[i] = 2  # acceleration for displacement DOFs
            if lmdofn[7]:
                order_time[7] = 1  # speed for pore pressure

        # Newmark params — use glb_control overrides if set, else analysis.dynamic
        if g.beta1 is not None:
            beeta1 = g.beta1
        elif a.dynamic:
            beeta1 = a.dynamic.beta1
        else:
            beeta1 = 0.25
        if g.beta2 is not None:
            beeta2 = g.beta2
        elif a.dynamic:
            beeta2 = a.dynamic.beta2
        else:
            beeta2 = 0.50
        if g.theta1 is not None:
            theta1 = g.theta1
        elif a.dynamic:
            theta1 = a.dynamic.theta1
        else:
            theta1 = 1.0

        # Count element groups per group
        group_nelems = {}
        if self.mesh:
            for eg in p.element_groups:
                count = 0
                for elem in self.mesh.elements.values():
                    if elem.group == eg.group:
                        count += 1
                group_nelems[eg.group] = count
        else:
            for eg in p.element_groups:
                group_nelems[eg.group] = 0

        # Has contact?
        ngaps = 0
        ngapb = 0
        if p.contact:
            ngaps = len(p.contact.gap_groups)

        # uwcpl
        uwcpl = 0
        for eg in p.element_groups:
            if eg.dof_type in (DofType.UW, DofType.UPC):
                uwcpl = 1
                break

        lines = []

        # Line 1: header text
        lines.append("NPOIN  npoinb NELEM  NDIMN  NMATS  NGROUP NTLINK outplot  KSTAB   MAT_curve meshc rmesh  ,level_set_problem,ljdp,stab_matde")

        # Line 2: values
        lines.append(f"{npoin} \t{npoinb}\t  {nelem}    {ndimn}      {nmats}    {ngroup}     0   {outplot}         {g.kstab}    0    0    0    0    0 9999")

        # Line 3-4: valv1/ndivide
        lines.append("valv1")
        lines.append("ndivide")

        # Line 5: text
        lines.append("NINIT KINIT winit NBLKS NLINK NONSY OUTIP outir outiw neumn equvs,type_ABC,block_sta nbackf nbspring ebody outind nbackd")

        # Line 6: control values
        nonsym = 0
        equvs = 0
        nbackf = 0
        nbspring = 0
        type_abc = "FIX"
        lines.append(
            f"     {g.ninit}    {g.kinit}     {g.winit}    {nblks}     0     {nonsym}     0     0     0     0     {equvs}   {type_abc}            0      {nbackf}        {nbspring}     0      0      0      0"
        )

        # Line 7: text
        lines.append("TYPE_PROBLEM TYPE_SOLVER TYPE_LOAD TYPE_NL STABPW nlayer,kglb,state_change,Bparameter balgor upliftin")

        # Line 8: problem type values
        bparam = 0
        if p.back_analysis and p.back_analysis.enabled:
            bparam = len(p.back_analysis.parameters)
        lines.append(
            f"{type_problem}         {type_solver}   {type_load}          {type_nl}    0    0    0    0    {bparam}    0    0"
        )

        # Line 9: text (layer info)
        lines.append("type_layer1")

        # Line 10: text
        lines.append("NMASS NSMAT NHMAT NQMAT NLDFL KGMAT NSWKW UWCPL NGRAV nflow ECWPIPE")

        # Line 11: values - use glb_control values
        nldfl = g.nmass  # NLDFL typically same as nmass
        lines.append(
            f"   {g.nmass}   {g.nsmat}   {g.nhmat}   {g.nqmat}   {nldfl}   0\t   {g.nswkw}\t  {uwcpl}\t   {ngrav}   0      0      ! =2,drained(u-pw); "
        )

        # Line 12: text (freeflownode)
        lines.append("nfreeflownode")

        # Line 13: text (temperature)
        lines.append("NTSMAT NTHMAT  KSTAT ground_inf src,nextrf,submodel")

        # Line 14: temperature values
        lines.append(f"       {g.ntsmat}       {g.nthmat}         {g.kstat}         0         0         0         0")

        # Line 15-16: MDOFN
        lines.append("MDOFN(/Ux Uy Uz Thxy Thyz Thzx P W A T)")
        lines.append(f"        {mdofn}")
        lines.append("         " + "         ".join(str(v) for v in lmdofn))
        lines.append("         " + "         ".join(str(v) for v in order_time))

        # Newmark params
        lines.append("BEETA1 BEETA2 THETA1")
        lines.append(f"     {beeta1:.3f}     {beeta2:.3f}     {theta1:.3f}")

        # equvs_process - use shorthand
        lines.append("equvs_process(1:ngroup)")
        lines.append("   1000*0 ")

        # appear_level - use shorthand
        lines.append("appear_level(1:ngroup)")
        lines.append("   1000*0 ")

        # APPEAR_PROCESS
        lines.append("APPEAR_PROCESS(1-0nblks/1-ngroup)")
        if p.construction_stages:
            for cs in p.construction_stages:
                appear = []
                for eg in p.element_groups:
                    if eg.group in cs.active_groups:
                        appear.append(1)
                    elif eg.group in cs.deactivate_groups:
                        appear.append(0)
                    else:
                        appear.append(1)
                lines.append("   " + " ".join(str(v) for v in appear))
        else:
            lines.append("   " + " ".join(["1"] * ngroup))

        # MATNO_PROCESS
        lines.append("MATNO_PROCESS(1-nblks/1-ngroup)")
        if p.construction_stages:
            for cs in p.construction_stages:
                matno_row = []
                for eg in p.element_groups:
                    changed = False
                    for mc in cs.material_change:
                        if mc.group == eg.group:
                            matno_row.append(mc.new_material)
                            changed = True
                            break
                    if not changed:
                        matno_row.append(eg.material)
                lines.append("   " + " ".join(str(v) for v in matno_row))
        else:
            lines.append("   " + " ".join(str(eg.material) for eg in p.element_groups))

        # force_process - use shorthand
        lines.append("force_process(1:ngroup)")
        lines.append("   1000*0")

        # average_appear
        lines.append("average_appear(1:ngroup)")
        lines.append(" " + " ".join(["-2"] * ngroup) + " -2 ")

        # GID output flags
        lines.append("gid_u,gid_s,gid_ms,gid_f,gid_rot,gid_v,gid_a,gid_T,gid_P,gid_Pv,gid_ep,gid_Y,gid_FC,gid_Ns,gid_Ss,gid_mxy,gid_bem,gid_wh,gid_wv,gid_bcs")
        results = set(p.output.results) if p.output.results else set()
        gid_u = 1 if "displacement" in results else 0
        gid_s = 1 if "stress" in results else 0
        gid_ms = 1 if "mises_stress" in results else 0
        gid_f = 1 if "force" in results or "stress" in results else 0
        gid_v = 1 if "velocity" in results else 0
        gid_a = 1 if "acceleration" in results else 0
        gid_Y = 1 if "yield" in results else 0
        gid_P = 1 if "pore_pressure" in results else 0
        gid_ep = 1 if "strain" in results else 0
        lines.append(
            f"   {gid_u}    {gid_s}      {gid_ms}     {gid_f}      0      {gid_v}     {gid_a}     0     {gid_P}      0      {gid_ep}     {gid_Y}      0      0      0      0      0      0      0      0"
        )

        # res output flags
        lines.append("res_u,res_s,res_ms,res_f,res_rot,res_v,res_a,res_T,res_P,res_Pv,res_ep,res_Y,res_FC,res_Ns,res_Ss,res_Tv,res_Pa")
        res_u = 1 if "displacement" in results else 0
        res_s = 1 if "stress" in results else 0
        res_ms = 1 if "mises_stress" in results else 0
        res_f = 1 if "force" in results or "stress" in results else 0
        lines.append(
            f"   {res_u}    {res_s}      {res_ms}     {res_f}      0      0     0     0     0      0      0     0      0      0      0      0      0"
        )

        # IFS water parameters
        lines.append("Icaddmass,swlifs2006,toth,ifswater,ifsgravity,absorb,alfa_p4,stiff_p4")
        # Format toth: integer if whole number
        toth_str = str(int(g.toth)) if g.toth == int(g.toth) else f"{g.toth}"
        lines.append(
            f"     {g.icaddmass}   {g.swlifs2006}     {toth_str}     {g.ifswater}     {g.ifsgravity:.3f}     {g.absorb:.3f}   {g.alfa_p4:.3f} {self._fmtg(g.stiff_p4)}"
        )

        # Steel/crack params
        lines.append("ftcrack,coefMpa,ikindks,doubsig,ktan1,  ktan2,  nlocalbeam,ndimnrt,listglocbeam(1:nlocalbeam),lelenrt(1:ndimnrt)")
        lines.append(f" {self._fmtg(g.ftcrack)} {self._fmtg(g.coefMpa)}     {g.ikindks}     {g.doubsig} {self._fmtg(g.ktan1)} {self._fmtg(g.ktan2)}     0    0")

        # MIF params
        lines.append("ntrans,nlaymif,epsMIFb,gamaMIF,ifixvar0_inpb,camif,dxmif")
        lines.append("     0     0 0.100E+01 0.200E-01     2 0.198E+04 0.250E+02")

        # hdam — always output max(nblks, 2) values
        lines.append("hdam")
        nvals = max(nblks, 2)
        lines.append(" ".join(["0."] * nvals))

        # water_level
        lines.append("water_level")
        lines.append(" ".join(["0."] * nvals))

        # modf_dis_blocks
        lines.append("modf_dis_blocks")
        lines.append("         " + " ".join(["0"] * nvals) + " ")

        # uinitial
        lines.append("uinitial")
        lines.append("        " + " ".join(["0"] * nvals))

        # backf section
        lines.append("backf()%")

        # links
        lines.append("1-ILINKS(I0, Freedom, node1,node2)")
        lines.append("1-tLINKS(I0, node1,node2)")

        # Group definitions header
        lines.append("1-NGROUP--GROUP INFORMATION")
        lines.append("INCLUDE (1) NAME KNAME INDEX CLASS NRFIELDS FIELDID SPECIAL")
        lines.append(" SPTYPE NELGROUP MATNO TYPE_ALGO TYPE_STIFF TYPE_ECOINT ilayer,elcod_local group_inf uplift_ic liquj")
        lines.append("\t\t (2) TYPE_MASS(1:NRFIELDS)(0-lumped)")
        lines.append("\t\t (3) for each field: nfdof-number of freedom,listdof(1:nfdof)")

        # Group definitions
        for eg in p.element_groups:
            etype = eg.element_type
            ename = ELEMENT_NAMES.get(etype, "Q4")
            index = ELEMENT_INDEX.get(etype, 4)
            nelgroup = group_nelems.get(eg.group, 0)
            matno = eg.material

            # DOF info
            dof_info = DOF_INFO.get(eg.dof_type, DOF_INFO[DofType.U])
            fieldid, nrfields, dof_fields = dof_info

            eclass = eg.element_class
            algo = eg.algo_type
            stiff = eg.stiffness_type

            gname = f"GROUP{eg.group}"
            ecoint = eg.gauss_order
            elcod_local = eg.elcod_local

            lines.append(
                f"{ename}   {gname}     {index}  {eclass}     {nrfields} {fieldid}    {algo}   {stiff}          "
                f"{nelgroup}  {matno}  0  1  {ecoint}  1 {elcod_local:.3E}  0  0  0"
            )

            # Type mass line
            mass_vals = []
            for _ in range(nrfields):
                mass_vals.append(str(eg.mass_type))
            alfa_str = f"{eg.rayleigh_alpha:.15E}" if eg.rayleigh_alpha != 0 else "0.000000000000000E+000"
            beta_str = f"{eg.rayleigh_beta}" if eg.rayleigh_beta != 0 else "0"
            lines.append(f"      {'  '.join(mass_vals)}  {alfa_str}   {beta_str}")

            # Order time line (for each field)
            lines.append("         0         0")

            # DOF list for each field
            for ifield in range(nrfields):
                if ifield == 0 and fieldid[0] == "U":
                    nfdof = ndimn
                    dofs = list(range(1, ndimn + 1))
                elif ifield == 0 and fieldid == "W":
                    nfdof = 1
                    dofs = [8]
                elif ifield == 0 and fieldid == "T":
                    nfdof = 1
                    dofs = [10]
                elif ifield == 1:
                    nfdof = 1
                    dofs = [8]
                else:
                    nfdof = ndimn
                    dofs = list(range(1, ndimn + 1))

                lines.append(f"         {nfdof}")
                lines.append("         " + "         ".join(str(d) for d in dofs))

        # tension_joint
        lines.append("tension_joint")
        lines.append("         0")

        # contact_joint
        lines.append("contact_joint")
        lines.append("         0")

        # Contact parameters
        if p.contact and ngaps > 0:
            ct = p.contact
            ctt_solver = ct.solver_type.value
            lines.append("ngaps ngapb ctt_pe miter_bt   torbt iblkbt nonsbt xlwsol method_gapi miter_state type_solver_ctt restart_ctt damp_ctt istatec")
            lines.append(
                f"  {ngaps}     {ngapb}     1     {ct.max_iterations}     {ct.tolerance:.1E}   1      0      0        0           1          {ctt_solver}           0       0.0E+00    1"
            )
        else:
            lines.append("ngaps ngapb ctt_pe miter_bt   torbt iblkbt nonsbt xlwsol method_gapi miter_state type_solver_ctt restart_ctt damp_ctt istatec")
            lines.append("  0     0     1     500     0.1E-05   1      0      0        0           1          PROFILE           0       0.0E+00    1")

        # nrcsteel
        lines.append(" nrcsteel")
        lines.append("         0")

        # nwcpipe
        lines.append(" nwcpipe")
        lines.append("         0")

        return "\n".join(lines) + "\n"
