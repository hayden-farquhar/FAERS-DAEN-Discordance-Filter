"""The 2x2 contingency table underlying every disproportionality metric.

                 event = Y      event != Y
    drug = X         a               b
    drug != X        c               d

`a` is the co-reported count (the cell the signalling rules gate on via `a >= 3`).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Contingency:
    a: float
    b: float
    c: float
    d: float

    def __post_init__(self) -> None:
        for name, v in (("a", self.a), ("b", self.b), ("c", self.c), ("d", self.d)):
            if v < 0:
                raise ValueError(f"contingency cell {name!r} must be >= 0, got {v}")

    @classmethod
    def from_marginals(cls, a: float, n_drug: float, n_event: float, n_total: float) -> "Contingency":
        """Build from the co-report count and the marginals (the form FAERS yields)."""
        b = n_drug - a
        c = n_event - a
        d = n_total - n_drug - n_event + a
        return cls(a, b, c, d)

    @property
    def n_drug(self) -> float:
        return self.a + self.b

    @property
    def n_event(self) -> float:
        return self.a + self.c

    @property
    def n_total(self) -> float:
        return self.a + self.b + self.c + self.d

    @property
    def expected(self) -> float:
        """Expected co-report count under reporting independence: n_drug * n_event / N."""
        N = self.n_total
        return (self.n_drug * self.n_event / N) if N else float("nan")

    def cells(self) -> tuple[float, float, float, float]:
        return (self.a, self.b, self.c, self.d)
