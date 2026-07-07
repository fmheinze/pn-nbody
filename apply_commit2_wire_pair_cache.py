#!/usr/bin/env python3
"""
Commit 2 helper: wire the ActiveList / PairCache infrastructure into eom.c and hamiltonian.c.

Run this from the pn-nbody repository root, after applying the Commit 1 infrastructure patch:

    python3 apply_commit2_wire_pair_cache.py
    make clean
    make USE_CUBA=0

The script is intentionally search/replace based, because it is more robust than a hand-written
multi-file patch when the previous patch was already applied locally.
"""

from pathlib import Path
import sys


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{label}: expected exactly one match, found {count}")
    return text.replace(old, new, 1)


def insert_include_once(text: str, anchor: str, include_line: str, label: str) -> str:
    if include_line in text:
        return text
    if anchor not in text:
        raise RuntimeError(f"{label}: include anchor not found")
    return text.replace(anchor, anchor + "\n" + include_line, 1)


repo_root = Path.cwd()
eom_path = repo_root / "src" / "eom.c"
ham_path = repo_root / "src" / "hamiltonian.c"

if not eom_path.exists() or not ham_path.exists():
    print("Please run this script from the pn-nbody repository root.", file=sys.stderr)
    sys.exit(1)


# --------------------------------------------------------------------------------------
# eom.c
# --------------------------------------------------------------------------------------

eom = eom_path.read_text()

eom = insert_include_once(
    eom,
    '#include "hamiltonian.h"',
    '#include "pair_cache.h"',
    "src/eom.c",
)

old_eom_init = """    // Masses
    double m[num_bodies];
    double inv_m[num_bodies];
    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) continue;
        m[a] = ode_params->masses[a];
        inv_m[a] = 1.0 / m[a];
    }

    // Momenta
    double p[num_bodies][num_dim];
    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) continue;
        for (int i = 0; i < num_dim; i++)
            p[a][i] = w[array_half + a * num_dim + i];
    }

    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double inv_r[num_bodies][num_bodies];
    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) continue;
        for (int b = a; b < num_bodies; b++) {
            if (!ode_params->active[b]) continue;
            for (int i = 0; i < num_dim; i++){
                x_rel[a][b][i] = w[a * num_dim + i] - w[b * num_dim + i];
                x_rel[b][a][i] = -x_rel[a][b][i];
            } 
            if (a == b) {
                r[a][b] = 0.0;
                inv_r[a][b] = 0.0;
            } else {
                r[a][b] = norm(x_rel[a][b], num_dim);
                inv_r[a][b] = 1.0 / r[a][b];
            }
            r[b][a] = r[a][b];
            inv_r[b][a] = inv_r[a][b];

            for (int i = 0; i < num_dim; i++){
                if (a == b){
                        n[a][b][i] = n[b][a][i] = 0.0;
                } else {
                        n[a][b][i] = x_rel[a][b][i] * inv_r[a][b];
                        n[b][a][i] = -n[a][b][i];
                }
            }
        }
    }
"""

new_eom_init = """    // Cached active-body, momentum and pair-geometry data.
    //
    // This keeps the current RHS algebra unchanged, but builds the common quantities
    // through the reusable PairCache infrastructure.  Later analytic 2PN terms can use
    // the same cache directly instead of rebuilding r_ab, n_ab and scalar products.
    PairCache cache;
    pair_cache_build(&cache, w, ode_params);

    // Masses
    double m[num_bodies];
    double inv_m[num_bodies];
    for (int ia = 0; ia < cache.active.num_active; ia++) {
        int a = cache.active.ids[ia];
        m[a] = cache.m[a];
        inv_m[a] = cache.inv_m[a];
    }

    // Momenta
    double p[num_bodies][num_dim];
    for (int ia = 0; ia < cache.active.num_active; ia++) {
        int a = cache.active.ids[ia];
        for (int i = 0; i < num_dim; i++)
            p[a][i] = pair_cache_p(&cache, a, i);
    }

    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double inv_r[num_bodies][num_bodies];
    for (int ia = 0; ia < cache.active.num_active; ia++) {
        int a = cache.active.ids[ia];
        for (int ib = ia; ib < cache.active.num_active; ib++) {
            int b = cache.active.ids[ib];

            for (int i = 0; i < num_dim; i++) {
                x_rel[a][b][i] = pair_cache_xrel(&cache, a, b, i);
                x_rel[b][a][i] = pair_cache_xrel(&cache, b, a, i);
                n[a][b][i] = pair_cache_n(&cache, a, b, i);
                n[b][a][i] = pair_cache_n(&cache, b, a, i);
            }

            r[a][b] = pair_cache_r(&cache, a, b);
            r[b][a] = pair_cache_r(&cache, b, a);
            inv_r[a][b] = pair_cache_inv_r(&cache, a, b);
            inv_r[b][a] = pair_cache_inv_r(&cache, b, a);
        }
    }
"""

eom = replace_once(eom, old_eom_init, new_eom_init, "src/eom.c initialization block")

old_eom_end = """    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) {
            for (int i = 0; i < num_dim; i++) {
                dwdt[a * num_dim + i] = 0.0;
                dwdt[array_half + a * num_dim + i] = 0.0;
            }
        }
    }
}
"""

new_eom_end = """    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) {
            for (int i = 0; i < num_dim; i++) {
                dwdt[a * num_dim + i] = 0.0;
                dwdt[array_half + a * num_dim + i] = 0.0;
            }
        }
    }

    pair_cache_free(&cache);
}
"""

eom = replace_once(eom, old_eom_end, new_eom_end, "src/eom.c rhs cleanup block")

eom_path.write_text(eom)


# --------------------------------------------------------------------------------------
# hamiltonian.c
# --------------------------------------------------------------------------------------

ham = ham_path.read_text()

ham = insert_include_once(
    ham,
    '#include "utils.h"',
    '#include "pair_cache.h"',
    "src/hamiltonian.c",
)

old_h0_h1 = """// Computes the 0PN (Newtonian) part of the N-body Hamiltonian
double H0PN(double* w, struct ode_params* ode_params)
{
    int num_bodies = ode_params->num_bodies_initial;
    int num_dim = ode_params->num_dim;
    int array_half = num_dim * num_bodies;
    double rel_dist2, p2;
    double H = 0;

    for (int a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) continue;
        // Kinetic energy
        p2 = 0.0;
        for (int i = 0; i < num_dim; i++)
            p2 += w[array_half + num_dim * a + i] * w[array_half + num_dim * a + i];
        H += p2/(2 * ode_params->masses[a]);

        // Potential energy
        for (int b = a+1; b < num_bodies; b++) {
            if (!ode_params->active[b]) continue;
            rel_dist2 = 0.0;
            for (int i = 0; i < num_dim; i++)
                rel_dist2 += (w[num_dim * a + i] - w[num_dim * b + i]) * 
                    (w[num_dim * a + i] - w[num_dim * b + i]);
            H -= ode_params->masses[a] * ode_params->masses[b] / sqrt(rel_dist2);
        }
    }
    return H;
}


// Computes the 1PN part of the N-body Hamiltonian
double H1PN(double* w, struct ode_params* ode_params) 
{
    int num_bodies = ode_params->num_bodies_initial;
    int num_dim = ode_params->num_dim;
    int array_half = num_dim * num_bodies;
    int a, b, c, i;
    double m_a, m_b, m_c;
    double pa_dot_pa, pa_dot_pb;
    double dx, r_ab, r_ac, ni, na_dot_pa, na_dot_pb;
    double p_ai, p_bi, dx_ac;
    double H = 0.0;

    // Compute kinetic and potential energy
    for (a = 0; a < num_bodies; a++) {
        if (!ode_params->active[a]) continue;
        m_a = ode_params->masses[a];
        pa_dot_pa = 0.0;

        for (i = 0; i < num_dim; i++) {
            p_ai = w[array_half + a * num_dim + i];
            pa_dot_pa += p_ai * p_ai;
        }

        H += -0.125 * m_a * pa_dot_pa * pa_dot_pa / (m_a * m_a * m_a * m_a);

        for (b = 0; b < num_bodies; b++) {
            if (b == a || !ode_params->active[b]) continue;

            m_b = ode_params->masses[b];
            r_ab = 0.0;
            pa_dot_pb = 0.0;
            na_dot_pa = 0.0;
            na_dot_pb = 0.0;

            for (i = 0; i < num_dim; i++) {
                dx = w[a * num_dim + i] - w[b * num_dim + i];
                r_ab += dx * dx;
            }

            r_ab = sqrt(r_ab);

            for (i = 0; i < num_dim; i++) {
                p_ai = w[array_half + a * num_dim + i];
                p_bi = w[array_half + b * num_dim + i];

                dx = w[a * num_dim + i] - w[b * num_dim + i];
                ni = dx / r_ab;

                pa_dot_pb += p_ai * p_bi;
                na_dot_pa += ni * p_ai;
                na_dot_pb += ni * p_bi;
            }

            H += -0.25 * m_a * m_b / r_ab * (6 * pa_dot_pa / (m_a * m_a) 
                - 7 * pa_dot_pb / (m_a * m_b) - (na_dot_pa * na_dot_pb) / (m_a * m_b));

            for (c = 0; c < num_bodies; c++) {
                if (c == a || !ode_params->active[c]) continue;

                m_c = ode_params->masses[c];
                r_ac = 0.0;

                for (i = 0; i < num_dim; i++) {
                dx_ac = w[a * num_dim + i] - w[c * num_dim + i];
                r_ac += dx_ac * dx_ac;
                }
                r_ac = sqrt(r_ac);

                H += 0.5 * m_a * m_b * m_c / (r_ab * r_ac);
            }
        }
    }
    return H;
}
"""

new_h0_h1 = """// Computes the 0PN (Newtonian) part of the N-body Hamiltonian
double H0PN(double* w, struct ode_params* ode_params)
{
    PairCache cache;
    pair_cache_build(&cache, w, ode_params);

    double H = 0.0;

    for (int ia = 0; ia < cache.active.num_active; ia++) {
        int a = cache.active.ids[ia];

        // Kinetic energy
        H += cache.p2[a] / (2.0 * cache.m[a]);

        // Potential energy
        for (int ib = ia + 1; ib < cache.active.num_active; ib++) {
            int b = cache.active.ids[ib];
            H -= cache.m[a] * cache.m[b] * pair_cache_inv_r(&cache, a, b);
        }
    }

    pair_cache_free(&cache);
    return H;
}


// Computes the 1PN part of the N-body Hamiltonian
double H1PN(double* w, struct ode_params* ode_params) 
{
    PairCache cache;
    pair_cache_build(&cache, w, ode_params);

    double H = 0.0;

    // Compute kinetic and potential energy
    for (int ia = 0; ia < cache.active.num_active; ia++) {
        int a = cache.active.ids[ia];

        double m_a = cache.m[a];
        double pa_dot_pa = cache.p2[a];

        H += -0.125 * m_a * pa_dot_pa * pa_dot_pa / (m_a * m_a * m_a * m_a);

        for (int ib = 0; ib < cache.active.num_active; ib++) {
            int b = cache.active.ids[ib];
            if (b == a) continue;

            double m_b = cache.m[b];
            double r_ab = pair_cache_r(&cache, a, b);
            double pa_dot_pb = pair_cache_p_dot(&cache, a, b);
            double na_dot_pa = pair_cache_n_dot_p(&cache, a, b, a);
            double na_dot_pb = pair_cache_n_dot_p(&cache, a, b, b);

            H += -0.25 * m_a * m_b / r_ab * (6.0 * pa_dot_pa / (m_a * m_a) 
                - 7.0 * pa_dot_pb / (m_a * m_b) - (na_dot_pa * na_dot_pb) / (m_a * m_b));

            for (int ic = 0; ic < cache.active.num_active; ic++) {
                int c = cache.active.ids[ic];
                if (c == a) continue;

                double m_c = cache.m[c];
                double r_ac = pair_cache_r(&cache, a, c);

                H += 0.5 * m_a * m_b * m_c / (r_ab * r_ac);
            }
        }
    }

    pair_cache_free(&cache);
    return H;
}
"""

ham = replace_once(ham, old_h0_h1, new_h0_h1, "src/hamiltonian.c H0PN/H1PN block")

ham_path.write_text(ham)

print("Updated src/eom.c and src/hamiltonian.c")
print("Next suggested commands:")
print("  make clean")
print("  make USE_CUBA=0")
print("  make USE_CUBA=1   # optional, if Cuba is installed")
