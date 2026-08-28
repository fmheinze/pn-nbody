#!/usr/bin/env python3
"""Validate symbolic 2PN triple forces against high-precision complex-step derivatives."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import mpmath as mp
import sympy as sp

import pn_expressions as model
from generate_2pn_triple_momentum import distinct_forces, reducible_forces


Vector = Sequence[mp.mpf | mp.mpc]
StateBuilder = Callable[[list[list[mp.mpf | mp.mpc]]], dict[sp.Symbol, mp.mpf | mp.mpc]]


def dot(a: Vector, b: Vector) -> mp.mpf | mp.mpc:
    return sum((x*y for x, y in zip(a, b)), mp.mpf("0"))


def separation(a: Vector, b: Vector) -> tuple[mp.mpf | mp.mpc, list[mp.mpf | mp.mpc]]:
    delta = [a[i] - b[i] for i in range(3)]
    radius = mp.sqrt(dot(delta, delta))
    return radius, [component/radius for component in delta]


def assign_vector(
    values: dict[sp.Symbol, mp.mpf | mp.mpc], symbols: Sequence[sp.Symbol], vector: Vector
) -> None:
    values.update(zip(symbols, vector))


def distinct_values(
    positions: list[list[mp.mpf | mp.mpc]], momenta: list[list[mp.mpf]], masses: list[mp.mpf]
) -> dict[sp.Symbol, mp.mpf | mp.mpc]:
    xa, xb, xc = positions
    pa, pb, pc = momenta
    ma, mb, mc = masses
    rab, nab = separation(xa, xb)
    rac, nac = separation(xa, xc)
    rbc, ncb = separation(xc, xb)
    u = [nab[i] + nac[i] for i in range(3)]
    v = [nab[i] + ncb[i] for i in range(3)]

    values: dict[sp.Symbol, mp.mpf | mp.mpc] = {
        model.ma: ma,
        model.mb: mb,
        model.mc: mc,
        model.ima: 1/ma,
        model.imb: 1/mb,
        model.imc: 1/mc,
        model.rab: rab,
        model.rac: rac,
        model.rbc: rbc,
        model.irab: 1/rab,
        model.irac: 1/rac,
        model.irbc: 1/rbc,
        model.perimeter: rab + rac + rbc,
        model.inv_perimeter: 1/(rab + rac + rbc),
        model.pa_sq: dot(pa, pa),
        model.pb_sq: dot(pb, pb),
        model.pc_sq: dot(pc, pc),
        model.papb: dot(pa, pb),
        model.papc: dot(pa, pc),
        model.pbpc: dot(pb, pc),
        model.Na: dot(nab, pa),
        model.Nb: dot(nab, pb),
        model.Gc: dot(nab, pc),
        model.Ec: dot(nac, pc),
        model.q: dot(nab, nac),
        model.upa: dot(u, pa),
        model.upc: dot(u, pc),
        model.vpa: dot(v, pa),
        model.vpb: dot(v, pb),
        model.vpc: dot(v, pc),
    }
    assign_vector(values, model.pa, pa)
    assign_vector(values, model.pb, pb)
    assign_vector(values, model.pc, pc)
    assign_vector(values, model.nab, nab)
    assign_vector(values, model.nac, nac)
    assign_vector(values, model.ncb, ncb)
    return values


def reducible_values(
    positions: list[list[mp.mpf | mp.mpc]], momenta: list[list[mp.mpf]], masses: list[mp.mpf]
) -> dict[sp.Symbol, mp.mpf | mp.mpc]:
    xa, xb = positions
    pa, pb = momenta
    ma, mb = masses
    rab, nab = separation(xa, xb)
    values: dict[sp.Symbol, mp.mpf | mp.mpc] = {
        model.ma: ma,
        model.mb: mb,
        model.ima: 1/ma,
        model.imb: 1/mb,
        model.rab: rab,
        model.irab: 1/rab,
        model.pa_sq: dot(pa, pa),
        model.pb_sq: dot(pb, pb),
        model.papb: dot(pa, pb),
        model.Na: dot(nab, pa),
        model.Nb: dot(nab, pb),
    }
    assign_vector(values, model.pa, pa)
    assign_vector(values, model.pb, pb)
    assign_vector(values, model.nab, nab)
    return values


def compile_expressions(
    expressions: Sequence[sp.Expr],
) -> Callable[[dict[sp.Symbol, mp.mpf | mp.mpc]], tuple[mp.mpf | mp.mpc, ...]]:
    symbols = sorted(set().union(*(expression.free_symbols for expression in expressions)), key=str)
    function = sp.lambdify(symbols, tuple(expressions), modules="mpmath")

    def evaluate(values: dict[sp.Symbol, mp.mpf | mp.mpc]) -> tuple[mp.mpf | mp.mpc, ...]:
        return tuple(function(*(values[symbol] for symbol in symbols)))

    return evaluate


def validate_case(
    name: str,
    hamiltonian: sp.Expr,
    forces: Sequence[sp.Expr],
    positions: list[list[mp.mpf]],
    value_builder: StateBuilder,
) -> mp.mpf:
    evaluate_hamiltonian = compile_expressions((hamiltonian,))
    evaluate_forces = compile_expressions(forces)
    analytic = evaluate_forces(value_builder(positions))
    step = mp.mpf("1e-40")
    worst = mp.mpf("0")

    for body in range(len(positions)):
        for axis in range(3):
            perturbed = [list(position) for position in positions]
            perturbed[body][axis] = mp.mpc(perturbed[body][axis], step)
            value = evaluate_hamiltonian(value_builder(perturbed))[0]
            complex_step_force = -mp.im(value)/step
            expected = analytic[3*body + axis]
            scale = max(mp.mpf("1"), abs(expected), abs(complex_step_force))
            worst = max(worst, abs(expected - complex_step_force)/scale)

    print(f"{name}: worst scaled complex-step error {mp.nstr(worst, 6)}")
    return worst


def main() -> None:
    mp.mp.dps = 80

    distinct_ab = distinct_forces()
    distinct_all = distinct_ab + tuple(
        -distinct_ab[axis] - distinct_ab[3 + axis] for axis in range(3)
    )
    reducible_a = reducible_forces()
    reducible_all = reducible_a + tuple(-component for component in reducible_a)

    distinct_momenta = [
        [mp.mpf("0.07"), mp.mpf("-0.03"), mp.mpf("0.02")],
        [mp.mpf("-0.04"), mp.mpf("0.06"), mp.mpf("0.01")],
        [mp.mpf("0.03"), mp.mpf("0.02"), mp.mpf("-0.05")],
    ]
    distinct_masses = [mp.mpf("0.8"), mp.mpf("1.1"), mp.mpf("0.9")]
    distinct_builder = lambda positions: distinct_values(
        positions, distinct_momenta, distinct_masses
    )
    distinct_error = validate_case(
        "distinct triple",
        model.ordered_triple_2pn_distinct_hamiltonian(),
        distinct_all,
        [
            [mp.mpf("-1.2"), mp.mpf("0.4"), mp.mpf("0.1")],
            [mp.mpf("0.7"), mp.mpf("-0.8"), mp.mpf("0.3")],
            [mp.mpf("1.4"), mp.mpf("0.9"), mp.mpf("-0.6")],
        ],
        distinct_builder,
    )

    reducible_momenta = [distinct_momenta[0], distinct_momenta[1]]
    reducible_masses = [distinct_masses[0], distinct_masses[1]]
    reducible_builder = lambda positions: reducible_values(
        positions, reducible_momenta, reducible_masses
    )
    reducible_error = validate_case(
        "reducible b = c triple",
        model.ordered_triple_2pn_reducible_hamiltonian(),
        reducible_all,
        [
            [mp.mpf("-0.9"), mp.mpf("0.2"), mp.mpf("0.5")],
            [mp.mpf("1.3"), mp.mpf("-0.7"), mp.mpf("-0.2")],
        ],
        reducible_builder,
    )

    tolerance = mp.mpf("1e-35")
    if max(distinct_error, reducible_error) > tolerance:
        raise SystemExit("complex-step validation failed")


if __name__ == "__main__":
    main()
