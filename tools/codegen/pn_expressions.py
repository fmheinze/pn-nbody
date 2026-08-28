"""Symbolic source of truth for post-Newtonian code generation.

For ordinary changes to the existing ordered-triple 2PN Hamiltonian, edit only the section marked
``EDITABLE HAMILTONIAN DEFINITIONS``. The generator differentiates those expressions, applies the
registered geometric chain rules, performs common-subexpression elimination, and emits C.

If a future Hamiltonian introduces a new scalar invariant, define its symbol and add its gradient
with respect to the affected unit vectors to ``UNIT_VECTOR_GRADIENTS``. Existing distance,
momentum-dot-product, and registered unit-vector invariants need no generator changes.
"""

from __future__ import annotations

from typing import Iterable

import sympy as sp


Vector = tuple[sp.Symbol, sp.Symbol, sp.Symbol]
R = sp.Rational


def vec(prefix: str) -> Vector:
    return sp.symbols(f"{prefix}0 {prefix}1 {prefix}2")


def dot(a: Iterable[sp.Expr], b: Iterable[sp.Expr]) -> sp.Expr:
    return sum((x*y for x, y in zip(a, b)), sp.S.Zero)


# Physical variables and cached scalar invariants available to Hamiltonian definitions.
ma, mb, mc = sp.symbols("ma mb mc", positive=True)
ima, imb, imc = sp.symbols("inv_ma inv_mb inv_mc", positive=True)
rab, rac, rbc = sp.symbols("r_ab r_ac r_bc", positive=True)
irab, irac, irbc = sp.symbols("inv_r_ab inv_r_ac inv_r_bc", positive=True)
perimeter, inv_perimeter = sp.symbols("perimeter inv_perimeter", positive=True)

pa = sp.symbols("pa0 pa1 pa2_component")
pb = sp.symbols("pb0 pb1 pb2_component")
pc = sp.symbols("pc0 pc1 pc2_component")
nab = vec("nab")
nac = vec("nac")
ncb = vec("ncb")

pa_sq, pb_sq, pc_sq = sp.symbols("pa_sq pb_sq pc_sq")
papb, papc, pbpc = sp.symbols("papb papc pbpc")
Na, Nb, Gc, Ec, q = sp.symbols("Na Nb Gc Ec q")
upa, upc, vpa, vpb, vpc = sp.symbols("upa upc vpa vpb vpc")


# Geometric registry used by the generic chain-rule engine. Each entry gives the vector gradient
# of a scalar invariant with respect to one unit separation vector.
UNIT_VECTOR_GRADIENTS: dict[Vector, dict[sp.Symbol, Vector]] = {
    nab: {
        Na: pa,
        Nb: pb,
        Gc: pc,
        q: nac,
        upa: pa,
        upc: pc,
        vpa: pa,
        vpb: pb,
        vpc: pc,
    },
    nac: {
        Ec: pc,
        q: nab,
        upa: pa,
        upc: pc,
    },
    ncb: {
        vpa: pa,
        vpb: pb,
        vpc: pc,
    },
}


# ================================================================================================
# EDITABLE HAMILTONIAN DEFINITIONS
# ================================================================================================

def ordered_triple_2pn_factorizable_hamiltonian() -> sp.Expr:
    """Two ordered-triple terms that also have a reducible b = c case."""
    first = (
        18*pa_sq*ima**2
        + 14*pb_sq*imb**2
        - 2*Nb**2*imb**2
        - 50*papb*ima*imb
        + 17*pbpc*imb*imc
        - 14*Na*Nb*ima*imb
        + 14*Nb*Gc*imb*imc
        + q*Nb*Ec*imb*imc
    )
    second = (
        2*(Na + Nb)*Ec*ima*imc
        + 5*q*pc_sq*imc**2
        - q*Ec**2*imc**2
        - 14*Gc*Ec*imc**2
    )
    return R(1, 8)*ma*mb*mc*(first/(rab*rac) + second/rab**2)


def ordered_triple_2pn_genuine_hamiltonian() -> sp.Expr:
    """Three genuine-triangle ordered-triple terms, defined only for distinct bodies."""
    triangle_momentum = (
        8*upa*vpc*ima*imc
        - 16*upc*vpa*ima*imc
        + 3*upa*vpb*ima*imb
        + 4*upc*vpc*imc**2
        + upa*vpa*ima**2
    )
    projected_momentum = (
        8*(papc - Na*Gc)*ima*imc
        - 3*(papb - Na*Nb)*ima*imb
        - 4*(pc_sq - Gc**2)*imc**2
        - (pa_sq - Na**2)*ima**2
    )
    distance_polynomial = (
        18*rab**2*rac**2
        - 60*rab**2*rbc**2
        - 24*rab**2*rac*(rab + rbc)
        + 60*rab*rac*rbc**2
        + 56*rab**3*rbc
        - 72*rab*rbc**3
        + 35*rbc**4
        + 6*rab**4
    )
    return (
        R(1, 2)*ma*mb*mc*triangle_momentum/(rab + rac + rbc)**2
        + R(1, 2)*ma*mb*mc*projected_momentum/((rab + rac + rbc)*rab)
        - R(1, 64)*ma**2*mb*mc*distance_polynomial/(rab**3*rac**3*rbc)
    )


# ================================================================================================
# DERIVED TOPOLOGY-SPECIFIC HAMILTONIANS -- NORMALLY NO EDITS BELOW THIS LINE
# ================================================================================================

REDUCIBLE_B_EQUALS_C = {
    mc: mb,
    imc: imb,
    rac: rab,
    pc_sq: pb_sq,
    pbpc: pb_sq,
    Gc: Nb,
    Ec: Nb,
    q: 1,
}


def ordered_triple_2pn_distinct_hamiltonian() -> sp.Expr:
    return (
        ordered_triple_2pn_factorizable_hamiltonian()
        + ordered_triple_2pn_genuine_hamiltonian()
    )


def ordered_triple_2pn_reducible_hamiltonian() -> sp.Expr:
    return ordered_triple_2pn_factorizable_hamiltonian().subs(
        REDUCIBLE_B_EQUALS_C, simultaneous=True
    )
