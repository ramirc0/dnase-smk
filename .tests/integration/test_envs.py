"""Check rule envs for library behavior the tools rely on."""

import subprocess

from conftest import ROOT


def env_python(env, code):
    """Run Python code in a rule env.

    Parameters
    ----------
    env : str
        pixi environment of ``workflow/envs/pixi.toml``.
    code : str
        Code to run.

    Returns
    -------
    str
        Standard output.
    """
    return subprocess.run(
        ["pixi", "run", "--manifest-path", str(ROOT / "workflow/envs/pixi.toml")]
        + ["-e", env, "python", "-c", code],
        check=True,
        capture_output=True,
        text=True,
    ).stdout


def test_hotspot3_betainc_underflow():
    """float32 betainc underflows to 0, which hotspot3 recomputes, not to NaN, which it rejects."""
    out = env_python(
        "hotspot3",
        "import numpy as np; from scipy.special import betainc; "
        "print(betainc(np.float32(980), np.float32(23), np.float32(0.2129), dtype=np.float32))",
    )
    assert float(out) == 0
