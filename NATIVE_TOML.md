# TOML-03 native Cook input

This branch starts at `a203e823e687f3c53e7a6d0aaaf2fe43a0811008`.
The same executable accepts the original no-argument invocation (reads `inp`),
or `hstar --input case.toml`. TOML mode reads TOML directly in Fortran; it does
not emit legacy control cards or use an in-memory formatted-card adapter.
Mesh `1.cor` and `1.ele` are resolved relative to the TOML file. Outputs are
written in the process working directory: always use a fresh run directory.

Only the fixed 16x16 Cook membrane **self-weight**, Q4 plane-strain profile from
HSTAR_NEXT `docs/fortran-interface/toml-01/examples/cooks_membrane` is supported.
This is not the literature Cook shear-load benchmark. The staged dam remains
available through legacy input; native staged input is deliberately rejected.

Editable fields: E (0,1e15] Pa, density (0,1e6] kg/m3, nu (-1,0.5), gravity
(0,1e4] m/s2, substeps and max_iterations integer [1,100]. Bounds define accepted
input, not engineering validation over every combination. All other values,
including reference names and metadata, are fixed by `HSTAR/NativeContract.inc`.
Unknown/missing/duplicate keys, unsupported tables, bad references and invalid
mesh rows are rejected before opening solver output, with exit 2. The TOML
subset uses one-line arrays and double-quoted unescaped strings, a 64 KiB file
limit and 1023-character lines; it is not a general TOML implementation.

`yl_authoring_toml.f90` reuses the standalone parser from HSTAR evolution
`c65bfe4`; this copy adds header tracking and strict token/length checks. It has
no ProblemState/RuntimeState dependency. `NativeInput` validates the bounded
contract; existing reader seams assign the validated material/load/control
values to the existing solver state. Shared allocations and numerical routines
remain in place. `native-consumed.txt` reports actual material, gravity and
step-control state at STATIC_U, after the corresponding readers ran.

Build with the existing `build_linux.sh release <new-output-directory>` under
Intel ifx + MKL. Set `GIDPOST_STUB` to HSTAR_NEXT's
`fem-chat/solver-repro/baseline/gidpost_stub.c`. The build script replaces its
output directory; never point it at a source tree or a previous run directory.
The Windows project includes the new sources, but Windows compilation has not
been tested in TOML-03.

Acceptance tools and records live in HSTAR_NEXT:
`tools/fortran-interface/{run_native,run_toml03_acceptance,check_toml03}.py`,
`test_toml_native.py`, and `docs/fortran-interface/toml-03/README.md`.
They use the existing authorized solver harness and a dedicated candidate binary.
No frontend, production binary, old deck or evolution source is modified.

Before extending the profile, retain the native no-card file-access gate,
full-field legacy/bridge/native comparisons, consumption probes and rejection
gates. Add only the already selected two-stage dam in TOML-04; do not turn this
bounded reader into a full-state migration or general legacy-card framework.
