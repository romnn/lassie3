"""Repair the OpenMP runtime conflict between Qseek's extensions and PyTorch.

Qseek's C extensions, and pyrocko's once built with OpenMP, link Homebrew's
`libomp`, while PyTorch ships and loads its own copy. Two LLVM OpenMP runtimes in one process corrupt each
other's thread-pool state, and Qseek segfaults inside a barrier the moment
PhaseNet touches a tensor.

Repointing Qseek's extensions at the copy PyTorch already loads leaves exactly
one runtime in the process. This cannot be done at build time, because PyTorch
is not installed yet when Qseek's wheel is built, so it runs as a post-install
repair. `uv sync` reinstalls Qseek from its cached wheel and reverts the fix, so
this is idempotent and safe to run before every search.
"""

from __future__ import annotations

import logging
import subprocess
import sysconfig
from pathlib import Path

logger = logging.getLogger(__name__)


def _site_packages() -> Path:
    return Path(sysconfig.get_paths()["purelib"])


def _linked_libraries(binary: Path) -> list[str]:
    output = subprocess.run(
        ["otool", "-L", str(binary)], capture_output=True, text=True, check=True
    ).stdout
    return [line.split(" (")[0].strip() for line in output.splitlines()[1:]]


def unify_openmp() -> list[Path]:
    """Point every OpenMP-linked extension at PyTorch's `libomp`; return those changed."""
    site = _site_packages()
    torch_omp = site / "torch" / "lib" / "libomp.dylib"
    # Pyrocko's extensions sit directly in its package directory; Qseek keeps
    # its own under `ext/`. Both packages are imported by the Qseek process.
    extensions = sorted(
        list((site / "qseek" / "ext").glob("*.so")) + list((site / "pyrocko").glob("*.so"))
    )

    if not torch_omp.exists() or not extensions:
        logger.debug("nothing to repair: torch or qseek extensions not installed")
        return []

    patched = []
    for extension in extensions:
        stale = [
            library
            for library in _linked_libraries(extension)
            if library.endswith("libomp.dylib") and Path(library) != torch_omp
        ]
        if not stale:
            continue

        for library in stale:
            subprocess.run(
                ["install_name_tool", "-change", library, str(torch_omp), str(extension)],
                check=True,
            )
        # install_name_tool invalidates the signature, and macOS refuses to load
        # an arm64 binary whose signature no longer matches.
        subprocess.run(["codesign", "-f", "-s", "-", str(extension)], check=True,
                       capture_output=True)
        patched.append(extension)
        logger.info("repointed %s to %s", extension.name, torch_omp)

    return patched


def doctor() -> None:
    """Run every environment repair, reporting what was needed."""
    patched = unify_openmp()
    if patched:
        print(f"Repaired OpenMP linkage for {len(patched)} extension(s): "
              + ", ".join(path.name for path in patched))
    else:
        print("Environment OK: every OpenMP extension shares PyTorch's runtime.")
