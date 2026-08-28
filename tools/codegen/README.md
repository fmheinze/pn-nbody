# Post-Newtonian code generation

This directory contains development-time symbolic tools. They are deliberately outside `src/`:
the integrator does not need Python or SymPy at runtime, while regenerating analytic expressions
does.

## Editing the 2PN ordered-triple Hamiltonian

The symbolic source of truth is `pn_expressions.py`. For changes involving the expressions already
defined there, edit only the section labelled `EDITABLE HAMILTONIAN DEFINITIONS`.

The generator automatically constructs the reducible and distinct Hamiltonians, differentiates
them with respect to pair distances and unit-vector invariants, applies the Cartesian chain rule,
derives the final body's force from translation invariance, performs common-subexpression
elimination, and updates the marked region in `src/eom/eom_conservative.c`.

If a future Hamiltonian introduces a new scalar invariant, add its symbol and unit-vector gradient
to `UNIT_VECTOR_GRADIENTS` in `pn_expressions.py`. Changing coefficients or combinations of existing
invariants needs no change to the generator.

## Commands

SymPy is required for all commands; its normal installation also provides `mpmath` for the
high-precision validator.

```sh
python3 -m pip install sympy
make regenerate-eom
make check-generated
make validate-generated
```

- `regenerate-eom` updates the generated C region.
- `check-generated` makes no changes and fails if the checked-in C differs from the symbolic
  source. It is suitable for CI.
- `validate-generated` first checks freshness, then compares the symbolic analytic forces with
  independent high-precision complex-step derivatives of both the distinct and reducible
  Hamiltonians.
