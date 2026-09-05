"""The step registry, and the label boundary it enforces.

A step is label-free when it cannot see the phenotype: filtering by gene name,
joining an annotation, computing a protein-language-model score. Such a step may
run ONCE, before the mart is built, and its result is valid for every fold.

A step is label-using when it reads the phenotype: ranking by chi-square, mutual
information, or any selection that scores features against the outcome. Such a
step must run INSIDE a training fold. Run it before the boundary and the column
universe is chosen while looking at the held-out lineage, which is selection
leakage -- milder than choosing the model's features that way, because it only
decides which columns exist, but the same error.

`assert_placement` is what makes that structural. A pipeline that tries to run a
label-using step in the pre-mart phase fails at configuration rather than
producing a mart that looks fine and is not.
"""

from __future__ import annotations

import sys as _sys

from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Step:
    name: str
    layer: str          # universe | annotate | derive | rank | select
    label_free: bool
    summary: str
    run: Callable

    @property
    def phase(self) -> str:
        return "pre-mart" if self.label_free else "per-fold"


_REGISTRY: dict[str, Step] = {}


def register(step: Step) -> Step:
    if step.name in _REGISTRY:
        raise ValueError(f"duplicate FE step {step.name!r}")
    _REGISTRY[step.name] = step
    return step


def get(name: str) -> Step:
    try:
        return _REGISTRY[name]
    except KeyError:
        known = ", ".join(sorted(_REGISTRY)) or "(none registered)"
        raise KeyError(f"unknown FE step {name!r}; registered: {known}") from None


def all_steps() -> list[Step]:
    return [_REGISTRY[k] for k in sorted(_REGISTRY)]


def assert_placement(names: list[str], *, phase: str) -> None:
    """Refuse a step placed on the wrong side of the label boundary."""
    if phase not in ("pre-mart", "per-fold"):
        raise ValueError(f"unknown phase {phase!r}")
    misplaced = [n for n in names if get(n).phase != phase and get(n).phase == "per-fold"]
    if misplaced:
        raise ValueError(
            f"step(s) {misplaced} read the phenotype label and cannot run in the "
            f"{phase} phase: the mart is built once and shared by every fold, so the "
            "step would see the held-out lineage. Move them to the per-fold phase "
            "(--fe_rank / --fe_select)."
        )


_BUILTIN_STEPS = ("callability_filter", "class_filter", "rank_association",
                  "rank_mic", "region_filter", "select_cross_model")


def load_builtin_steps(strict: bool = False) -> list[str]:
    """Register the shipped steps, and report any that could not be imported.

    Imported one at a time rather than as one statement. A single import of all
    six means a step nobody asked for can end the run: rank_mic imports scipy,
    which the feature-engineering image did not carry, so a run with
    `fe_pre_steps` EMPTY still died at registry load. A dependency missing for an
    unrequested step is a reduced menu, not a failure.

    The unavailable steps are returned and named on stderr rather than swallowed,
    because a step silently missing from the registry would let a requested step
    be reported as unknown, which is a confusing way to say "not installed".
    """
    import importlib
    unavailable = []
    for name in _BUILTIN_STEPS:
        try:
            importlib.import_module(f"analysis.scripts.fe_steps.{name}")
        except ImportError as exc:
            if strict:
                raise
            unavailable.append(name)
            print(f"[fe_steps] step {name!r} unavailable: {exc}", file=_sys.stderr)
    return unavailable
