"""Command-line interface for hstar-yaml."""

from pathlib import Path

import click
import yaml

from .schema import HstarProject
from .validator import validate_project
from .mesh_reader import read_gid_mesh, read_cor_ele
from .generators import generate_all


@click.command()
@click.argument("yaml_file", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output-dir", type=click.Path(path_type=Path), default=None,
              help="Output directory (default: same as YAML file)")
@click.option("--validate", "validate_only", is_flag=True,
              help="Only validate, do not generate files")
@click.option("-v", "--verbose", is_flag=True, help="Verbose output")
def main(yaml_file: Path, output_dir: Path | None, validate_only: bool, verbose: bool):
    """Generate HSTAR input files from a YAML project file.

    \b
    Usage:
        hstar-yaml project.yaml
        hstar-yaml project.yaml -o ./run/
        hstar-yaml project.yaml --validate
    """
    if output_dir is None:
        output_dir = yaml_file.parent

    # Load YAML
    click.echo(f"Reading {yaml_file} ...")
    with open(yaml_file, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    # Parse into schema
    try:
        project = HstarProject.model_validate(raw)
    except Exception as e:
        click.echo(f"Schema error / 格式错误:\n{e}", err=True)
        raise SystemExit(1)

    if verbose:
        click.echo(f"  title: {project.title}")
        click.echo(f"  problem_name: {project.problem_name}")

    # Read mesh if specified
    mesh_data = None
    if project.mesh:
        if project.mesh.cor_file and project.mesh.ele_file:
            # Use .cor + .ele format
            cor_path = yaml_file.parent / project.mesh.cor_file
            ele_path = yaml_file.parent / project.mesh.ele_file
            if cor_path.exists() and ele_path.exists():
                click.echo(f"Reading mesh {cor_path} + {ele_path} ...")
                mesh_data = read_cor_ele(cor_path, ele_path, project.mesh.dimension)
                click.echo(f"  {mesh_data.npoin} nodes, {mesh_data.nelem} elements, "
                           f"{len(mesh_data.groups)} groups")
            else:
                if not cor_path.exists():
                    click.echo(f"Warning: cor file {cor_path} not found", err=True)
                if not ele_path.exists():
                    click.echo(f"Warning: ele file {ele_path} not found", err=True)
        elif project.mesh.file:
            mesh_path = yaml_file.parent / project.mesh.file
            if mesh_path.exists():
                click.echo(f"Reading mesh {mesh_path} ...")
                mesh_data = read_gid_mesh(mesh_path, project.mesh.dimension)
                click.echo(f"  {mesh_data.npoin} nodes, {mesh_data.nelem} elements, "
                           f"{len(mesh_data.groups)} groups")
            else:
                click.echo(f"Warning: mesh file {mesh_path} not found", err=True)

    # Validate
    click.echo("Validating / 验证中 ...")
    errors, warnings = validate_project(project, mesh_data)

    for w in warnings:
        click.echo(f"  WARNING: {w}")
    for e in errors:
        click.echo(f"  ERROR: {e}", err=True)

    if errors:
        click.echo(f"\n{len(errors)} error(s) found. Fix them before generating files.")
        click.echo(f"发现 {len(errors)} 个错误，请修正后重新运行。")
        raise SystemExit(1)

    if validate_only:
        click.echo("Validation passed. / 验证通过。")
        return

    # Generate
    click.echo(f"Generating files in {output_dir} ...")
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = generate_all(project, mesh_data, output_dir)

    for path in generated:
        click.echo(f"  -> {path.name}")

    click.echo(f"\nDone! Generated {len(generated)} files. / 完成！生成了 {len(generated)} 个文件。")


if __name__ == "__main__":
    main()
