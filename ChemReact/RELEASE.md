# Release Guide

## Minimal Release Prerequisites

- `ChemReact` folder is complete and validated.
- `README.md` and `USAGE.md` are up to date.
- Sample pipeline run passes.

## Pre-Release Check

```bash
python skills/ChemReact/scripts/run_closed_loop.py \
  --output-dir skills/ChemReact/sample_output_template \
  --emit-vis-plan-template-only

python skills/ChemReact/scripts/run_closed_loop.py \
  --target-smiles "CCCCC1=NC(Cl)=C(CO)N1CC2=CC=C(C3=CC=CC=C3C4=NNN=N4)C=C2" \
  --routes-file skills/ChemReact/samples/routes.json \
  --strategy-file skills/ChemReact/samples/strategy.json \
  --vis-plan-file skills/ChemReact/samples/vis_plan.json \
  --output-dir skills/ChemReact/sample_output
```

## Version Tagging (Local)

If repository is not initialized:

```bash
git init
git add .
git commit -m "feat: initial ChemReact release"
```

Tag release:

```bash
git tag -a v0.1.0 -m "ChemReact v0.1.0"
```

## Publish to Platform

1. Create remote repository.
2. Push branch and tag:

```bash
git remote add origin <your-repo-url>
git branch -M main
git push -u origin main
git push origin v0.1.0
```

3. Create platform Release named `v0.1.0`.
4. Attach:
- `README.md`
- `USAGE.md`
- sample outputs/screenshot

## Release Notes Template

```
## ChemReact v0.1.0
- Added closed-loop retrosynthesis pipeline
- Added schema-validated input contracts
- Added host adapters for Claude Code / OpenCode / Cursor
- Added visualization tree with orthogonal target-to-precursor connectors
```

