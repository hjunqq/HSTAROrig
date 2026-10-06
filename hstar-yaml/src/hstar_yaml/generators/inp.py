"""Generate the `inp` control file.

Format (from Fem.f90 L92-100):
    Line 1: 'input'
    Line 2: restart relis sysrelis ADINA Uopt_R gamax  (all 0 for new analysis)
    Line 3: 'input'
    Line 4: problem_name
    Line 5: (blank or '1' for single run)
"""

from __future__ import annotations

from .base import BaseGenerator


class InpGenerator(BaseGenerator):
    extension = ""  # file is literally named "inp"

    def build(self) -> str:
        lines = [
            "input",
            "0 0 0 0 0 0  ! restart,relis,sysrelis,ADINA,Uopt_R,gamax",
            "input",
            self.project.problem_name,
            "1",
            "",
        ]
        return "\n".join(lines)

    def generate(self, output_dir):
        content = self.build()
        path = output_dir / "inp"
        path.write_text(content, encoding="utf-8")
        return path
