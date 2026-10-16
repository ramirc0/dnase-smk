# CLAUDE.md

ENCODE DNase-seq pipeline v3.0.0, ported to Snakemake with maintained tools, with
ENCODE's 2025 hotspot3 and footprint outputs. README.md holds the layout, steps,
differences from ENCODE and conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run test                                # pytest -n 16 on the chr22 fixture
pixi run snakemake -n -p --profile profiles/local --configfile .tests/integration/config.yaml
pixi run snakemake --profile profiles/slurm --configfile .tests/validation/config.yaml
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` holds the rule envs. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- New rule files MUST sit in `workflow/rules/`, since the env paths resolve from there.
- The `onstart` hook installs every rule env before any job. Snakemake keeps only one `onstart:`.

## Gotchas

- Shell blocks MUST call scripts as `{workflow.basedir:q}/scripts/x.py`.
- Snakemake does not track scripts called from `shell:`. After a script edit, force its rule with `-R`.
- A process substitution in a bash assignment closes when the assignment ends. Put `<(...)` on the command's own line.
- In the hotspot3 env, pandas MUST stay below 3 and scipy below 1.14. `test_envs.py` checks the scipy one.
- bioconda's pysam has no PyPI mapping in pixi. Envs with PyPI packages that need pysam take pysam from PyPI.
- `no-build-isolation` MUST sit in `[feature.footprints.pypi-options]`. The workspace-level table does not apply.
- `ftd learn_dm` fails on shallow samples. The fixture sets `footprints` to `false` for UHT and RWU.
- Adapters come from ENCODE's trimming reports, not the v3 README.
- Auto mode may block `pixi lock`, `pixi install`, or SLURM submission. Ask the user to run them.

## Verifying a change

`pixi run test` runs every rule on the chr22 fixture, then ENCODE's v3.0.0 tasks
on the same inputs. Build the fixture with `.tests/download.sh` and
`make_fixture.sh` first. `conftest.encode()` reruns a task when its command or
inputs change. A `temp()` output added to a rule does not rerun the rule.
