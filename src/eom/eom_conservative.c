/**
 * @file eom_conservative.c
 * @brief Conservative 0PN, 1PN, and 2PN equations-of-motion kernels.
 */

#include <string.h>

#include "cache.h"
#include "eom_internal.h"
#include "utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void add_velocity_0pn(const PairCache *cache, double *velocity, double *lower_xdot)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];
        for (int i = 0; i < D; i++) {
            const int idx = a * D + i;
            const double value = pair_cache_p(cache, a, i) * cache->inv_m[a];
            if (velocity != NULL) velocity[idx] += value;
            if (lower_xdot != NULL) lower_xdot[idx] += value;
        }
    }
}


void add_momentum_0pn(const PairCache *cache, double *momentum, double *lower_pdot)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];
        const double ma = cache->m[a];

        for (int ib = ia + 1; ib < cache->active.num_active; ib++) {
            const int b = cache->active.ids[ib];
            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double pref = -ma * cache->m[b] * inv_r * inv_r;

            for (int i = 0; i < D; i++) {
                const double value = pref * pair_cache_n(cache, a, b, i);
                const int ai = a * D + i;
                const int bi = b * D + i;
                if (momentum != NULL) {
                    momentum[ai] += value;
                    momentum[bi] -= value;
                }
                if (lower_pdot != NULL) {
                    lower_pdot[ai] += value;
                    lower_pdot[bi] -= value;
                }
            }
        }
    }
}


void add_velocity_1pn(const PairCache *cache, double *velocity, double *lower_xdot)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];
        const double inv_ma = cache->inv_m[a];
        const double self_coeff = -0.5 * cache->p2[a] * inv_ma * inv_ma * inv_ma;

        for (int i = 0; i < D; i++) {
            const int idx = a * D + i;
            const double value = self_coeff * pair_cache_p(cache, a, i);
            if (velocity != NULL) velocity[idx] += value;
            if (lower_xdot != NULL) lower_xdot[idx] += value;
        }

        for (int ib = 0; ib < cache->active.num_active; ib++) {
            const int b = cache->active.ids[ib];
            if (b == a)
                continue;

            const double pref = -0.5 * pair_cache_inv_r(cache, a, b);
            const double nab_dot_pb = pair_cache_n_dot_p(cache, a, b, b);

            for (int i = 0; i < D; i++) {
                const int idx = a * D + i;
                const double value = pref * (
                    6.0 * cache->m[b] * inv_ma * pair_cache_p(cache, a, i)
                  - 7.0 * pair_cache_p(cache, b, i)
                  - nab_dot_pb * pair_cache_n(cache, a, b, i));
                if (velocity != NULL) velocity[idx] += value;
                if (lower_xdot != NULL) lower_xdot[idx] += value;
            }
        }
    }
}


void add_momentum_1pn(const PairCache *cache, double *momentum, double *lower_pdot)
{
    const ActiveList *active = &cache->active;
    const int D = cache->num_dim;
    const int N = cache->num_bodies;
    double Phi[N];

    for (int ia = 0; ia < active->num_active; ia++)
        Phi[active->ids[ia]] = 0.0;

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];
        for (int ib = ia + 1; ib < active->num_active; ib++) {
            const int b = active->ids[ib];
            const double inv_r = pair_cache_inv_r(cache, a, b);
            Phi[a] += cache->m[b] * inv_r;
            Phi[b] += cache->m[a] * inv_r;
        }
    }

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];
        const double ma = cache->m[a];
        const double inv_ma = cache->inv_m[a];

        for (int ib = 0; ib < active->num_active; ib++) {
            const int b = active->ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];
            const double inv_mb = cache->inv_m[b];
            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double pref = -0.5 * inv_r * inv_r;
            const double nab_dot_pa = pair_cache_n_dot_p(cache, a, b, a);
            const double nab_dot_pb = pair_cache_n_dot_p(cache, a, b, b);
            const double radial_scalar = pref * (
                3.0 * mb * inv_ma * cache->p2[a]
              + 3.0 * ma * inv_mb * cache->p2[b]
              - 7.0 * pair_cache_p_dot(cache, a, b)
              - 3.0 * nab_dot_pa * nab_dot_pb)
              + ma * mb * inv_r * inv_r * (Phi[a] + Phi[b]);

            for (int i = 0; i < D; i++) {
                const int idx = a * D + i;
                const double value = radial_scalar * pair_cache_n(cache, a, b, i)
                    + pref * (nab_dot_pb * pair_cache_p(cache, a, i)
                    + nab_dot_pa * pair_cache_p(cache, b, i));
                if (momentum != NULL) momentum[idx] += value;
                if (lower_pdot != NULL) lower_pdot[idx] += value;
            }
        }
    }
}


// ------------------------------------------------------------------------------------------------
// 2PN analytic building blocks
// ------------------------------------------------------------------------------------------------

/**
 * @brief Add the analytic one-body 2PN p^6 Hamiltonian contribution to dw/dt.
 *
 * The one-body 2PN Hamiltonian term is
 *
 *     H_2PN,1body = sum_a (p_a^2)^3 / (16 m_a^5).
 *
 * Therefore
 *
 *     dx_a^i/dt += dH/dp_{a i} = (3/8) (p_a^2)^2 p_a^i / m_a^5,
 *
 * and there is no dp/dt contribution because the term is independent of positions.
 */
void add_velocity_2pn_onebody_analytic(const PairCache *cache, double *velocity)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];

        const double inv_m = cache->inv_m[a];
        const double inv_m5 = inv_m * inv_m * inv_m * inv_m * inv_m;
        const double p2 = cache->p2[a];

        const double coeff = 3.0 / 8.0 * p2 * p2 * inv_m5;

        for (int i = 0; i < D; i++)
            velocity[a * D + i] += coeff * pair_cache_p(cache, a, i);
    }
}


// Add the velocity part of the ordered-pair 2PN contribution to the velocity.
void add_velocity_2pn_pair_analytic(const PairCache *cache, double *velocity)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];
        const double ma = cache->m[a];
        const double inv_ma = cache->inv_m[a];
        const double ma2 = ma * ma;
        const double pa2 = cache->p2[a];

        for (int ib = 0; ib < cache->active.num_active; ib++) {
            const int b = cache->active.ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];
            const double inv_mb = cache->inv_m[b];
            const double inv_mamb = inv_ma * inv_mb;
            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double inv_r2 = inv_r * inv_r;
            const double pb2 = cache->p2[b];
            const double papb = pair_cache_p_dot(cache, a, b);
            const double Na = pair_cache_n_dot_p(cache, a, b, a);
            const double Nb = pair_cache_n_dot_p(cache, a, b, b);

            const double dF1_dpa2 = 20.0 * mb * inv_ma * inv_ma * inv_ma * pa2
                - 11.0 * inv_mamb * pb2 + 10.0 * inv_mamb * Nb * Nb;
            const double dF1_dpb2 = -11.0 * inv_mamb * pa2;
            const double dF1_dpapb = -4.0 * inv_mamb * papb - 12.0 * inv_mamb * Na * Nb;
            const double dF1_dNa = -12.0 * inv_mamb * papb * Nb - 6.0 * inv_mamb * Na * Nb * Nb;
            const double dF1_dNb = 20.0 * inv_mamb * pa2 * Nb
                - 12.0 * inv_mamb * papb * Na - 6.0 * inv_mamb * Na * Na * Nb;

            const double dH_dpa2 = (1.0 / 16.0) * inv_r * dF1_dpa2 + 0.25 * mb * inv_r2;
            const double dH_dpb2 = (1.0 / 16.0) * inv_r * dF1_dpb2 + 0.25 * ma2 * inv_mb * inv_r2;
            const double dH_dpapb = (1.0 / 16.0) * inv_r * dF1_dpapb - 0.5 * ma * inv_r2;
            const double dH_dNa = (1.0 / 16.0) * inv_r * dF1_dNa;
            const double dH_dNb = (1.0 / 16.0) * inv_r * dF1_dNb;

            for (int i = 0; i < D; i++) {
                const double pai = pair_cache_p(cache, a, i);
                const double pbi = pair_cache_p(cache, b, i);
                const double ni = pair_cache_n(cache, a, b, i);

                velocity[a * D + i] += 2.0 * dH_dpa2 * pai + dH_dpapb * pbi + dH_dNa * ni;
                velocity[b * D + i] += 2.0 * dH_dpb2 * pbi + dH_dpapb * pai + dH_dNb * ni;
            }
        }
    }
}


/**
 * @brief Add the momentum part of the ordered-pair 2PN Hamiltonian contribution to dp/dt.
 *
 * This implements the ordered pair terms from H2PN/H2PN_base_complex:
 *
 *   H_ab = 1/(16 r_ab) * F1 + 1/4 * m_a^2 m_b / r_ab^2 * F2 - 1/4 * m_a^2 m_b^2 / r_ab^3,
 *
 * with
 *
 *   F1 = 10 m_a m_b (p_a^2/m_a^2)^2
 *     - 11 (p_a^2/(m_a m_b)) p_b^2
 *     -  2 (p_a.p_b)^2/(m_a m_b)
 *     + 10 (p_a^2/(m_a m_b)) (n_ab.p_b)^2
 *     - 12 (p_a.p_b/(m_a m_b)) (n_ab.p_a)(n_ab.p_b)
 *     -  3 (n_ab.p_a)^2 (n_ab.p_b)^2/(m_a m_b),
 *
 *   F2 = p_a^2/m_a^2 + p_b^2/m_b^2 - 2 (p_a.p_b)/(m_a m_b).
 */
void add_momentum_2pn_pair_analytic(const PairCache *cache, double *momentum)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];

        const double ma = cache->m[a];
        const double inv_ma = cache->inv_m[a];
        const double ma2 = ma * ma;
        const double inv_ma2 = inv_ma * inv_ma;

        const double pa2 = cache->p2[a];

        for (int ib = 0; ib < cache->active.num_active; ib++) {
            const int b = cache->active.ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];
            const double inv_mb = cache->inv_m[b];
            const double mb2 = mb * mb;
            const double inv_mb2 = inv_mb * inv_mb;

            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double inv_r2 = inv_r * inv_r;
            const double inv_r3 = inv_r2 * inv_r;
            const double inv_r4 = inv_r2 * inv_r2;

            const double pb2 = cache->p2[b];
            const double papb = pair_cache_p_dot(cache, a, b);
            const double Na = pair_cache_n_dot_p(cache, a, b, a);
            const double Nb = pair_cache_n_dot_p(cache, a, b, b);

            // Scalar pair Hamiltonian pieces
            const double inv_mamb = inv_ma * inv_mb;

            const double F1 =
                10.0 * mb * inv_ma * inv_ma * inv_ma * pa2 * pa2
              - 11.0 * inv_mamb * pa2 * pb2
              -  2.0 * inv_mamb * papb * papb
              + 10.0 * inv_mamb * pa2 * Nb * Nb
              - 12.0 * inv_mamb * papb * Na * Nb
              -  3.0 * inv_mamb * Na * Na * Nb * Nb;

            const double F2 =
                pa2 * inv_ma2
              + pb2 * inv_mb2
              - 2.0 * papb * inv_mamb;

            // Only the position-dependent invariant derivatives are required for dp/dt.
            const double dF1_dNa = -12.0 * inv_mamb * papb * Nb - 6.0 * inv_mamb * Na * Nb * Nb;

            const double dF1_dNb = 20.0 * inv_mamb * pa2 * Nb - 12.0 * inv_mamb * papb * Na
                - 6.0 * inv_mamb * Na * Na * Nb;

            const double dH_dNa = (1.0 / 16.0) * inv_r * dF1_dNa;

            const double dH_dNb = (1.0 / 16.0) * inv_r * dF1_dNb;

            const double dH_dr = - (1.0 / 16.0) * inv_r2 * F1 - 0.5 * ma2 * mb * inv_r3 * F2
                + 0.75 * ma2 * mb2 * inv_r4;

            /*
             * Vector chain rule.
             *
             * Position derivatives:
             *
             *   dH/dx_a^i = dH/dr n_i + dH/dNa * (p_a^i - Na n_i)/r + dH/dNb * (p_b^i - Nb n_i)/r
             *
             * and dH/dx_b^i = -dH/dx_a^i.
             */
            for (int i = 0; i < D; i++) {
                const double pai = pair_cache_p(cache, a, i);
                const double pbi = pair_cache_p(cache, b, i);
                const double ni = pair_cache_n(cache, a, b, i);

                const double dNa_dxa_i = inv_r * (pai - Na * ni);
                const double dNb_dxa_i = inv_r * (pbi - Nb * ni);
                const double dH_dxa_i = dH_dr * ni + dH_dNa * dNa_dxa_i + dH_dNb * dNb_dxa_i;

                // dp/dt = -dH/dx
                momentum[a * D + i] -= dH_dxa_i;
                momentum[b * D + i] += dH_dxa_i;
            }
        }
    }
}


// Add the velocity part of the ordered-triple 2PN contribution to the velocity.
void add_velocity_2pn_triple_analytic(const PairCache *cache, double *velocity)
{
    const ActiveList *active = &cache->active;
    const int D = cache->num_dim;
    const int na = active->num_active;
    const int len = cache->array_half;
    double *const buf = cache->triple_accum;

    memset(buf, 0, (size_t)na*(size_t)len*sizeof(*buf));

#ifdef _OPENMP
    const int threads = cache->triple_accum_threads;
    const int run_parallel = (threads > 1) && !omp_in_parallel();

    #pragma omp parallel for schedule(static) num_threads(threads) if (run_parallel)
#endif
    for (int ia = 0; ia < na; ia++) {
        double *const velocity_row = buf + (size_t)ia*(size_t)len;
        const int a = active->ids[ia];
        const double ma = cache->m[a];
        const double inv_ma = cache->inv_m[a];
        const double inv_ma2 = inv_ma * inv_ma;

        for (int ib = 0; ib < active->num_active; ib++) {
            const int b = active->ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];
            const double inv_mb = cache->inv_m[b];
            const double inv_mb2 = inv_mb * inv_mb;
            const double inv_rab = pair_cache_inv_r(cache, a, b);

            for (int ic = 0; ic < active->num_active; ic++) {
                const int c = active->ids[ic];
                if (c == a)
                    continue;

                const double mc = cache->m[c];
                const double inv_mc = cache->inv_m[c];
                const double inv_mc2 = inv_mc * inv_mc;
                const double inv_rac = pair_cache_inv_r(cache, a, c);
                const double inv_mamb = inv_ma * inv_mb;
                const double inv_mamc = inv_ma * inv_mc;
                const double inv_mbmc = inv_mb * inv_mc;

                const double Na = pair_cache_n_dot_p(cache, a, b, a);
                const double Nb = pair_cache_n_dot_p(cache, a, b, b);
                const double Gc = pair_cache_n_dot_p(cache, a, b, c);
                const double Ec = pair_cache_n_dot_p(cache, a, c, c);

                double q = 0.0;
                for (int i = 0; i < D; i++)
                    q += pair_cache_n(cache, a, b, i) * pair_cache_n(cache, a, c, i);

                const double pref1 = 0.125 * ma * mb * mc * inv_rab * inv_rac;
                const double pref2 = 0.125 * ma * mb * mc * inv_rab * inv_rab;

                for (int i = 0; i < D; i++) {
                    const double pai = pair_cache_p(cache, a, i);
                    const double pbi = pair_cache_p(cache, b, i);
                    const double pci = pair_cache_p(cache, c, i);
                    const double nab_i = pair_cache_n(cache, a, b, i);
                    const double nac_i = pair_cache_n(cache, a, c, i);

                    velocity_row[a * D + i] += pref1 * (
                        36.0 * pai * inv_ma2
                      - 50.0 * pbi * inv_mamb
                      - 14.0 * Nb * nab_i * inv_mamb);

                    velocity_row[b * D + i] += pref1 * (
                        28.0 * pbi * inv_mb2
                      - 4.0 * Nb * nab_i * inv_mb2
                      - 50.0 * pai * inv_mamb
                      + 17.0 * pci * inv_mbmc
                      - 14.0 * Na * nab_i * inv_mamb
                      + 14.0 * Gc * nab_i * inv_mbmc
                      + q * Ec * nab_i * inv_mbmc);

                    velocity_row[c * D + i] += pref1 * (
                        17.0 * pbi * inv_mbmc
                      + 14.0 * Nb * nab_i * inv_mbmc
                      + q * Nb * nac_i * inv_mbmc);

                    velocity_row[a * D + i] += pref2 * (2.0 * Ec * nab_i * inv_mamc);
                    velocity_row[b * D + i] += pref2 * (2.0 * Ec * nab_i * inv_mamc);
                    velocity_row[c * D + i] += pref2 * (
                        2.0 * (Na + Nb) * nac_i * inv_mamc
                      + 10.0 * q * pci * inv_mc2
                      - 2.0 * q * Ec * nac_i * inv_mc2
                      - 14.0 * (Ec * nab_i + Gc * nac_i) * inv_mc2);
                }

                if (c == b)
                    continue;

                const double s = pair_cache_r(cache, a, b) + pair_cache_r(cache, b, c) 
                    + pair_cache_r(cache, c, a);
                double upa = 0.0, upc = 0.0;
                double vpa = 0.0, vpb = 0.0, vpc = 0.0;

                for (int i = 0; i < D; i++) {
                    const double ui = pair_cache_n(cache, a, b, i) + pair_cache_n(cache, a, c, i);
                    const double vi = pair_cache_n(cache, a, b, i) + pair_cache_n(cache, c, b, i);
                    upa += ui * pair_cache_p(cache, a, i);
                    upc += ui * pair_cache_p(cache, c, i);
                    vpa += vi * pair_cache_p(cache, a, i);
                    vpb += vi * pair_cache_p(cache, b, i);
                    vpc += vi * pair_cache_p(cache, c, i);
                }

                const double pref3 = 0.5 * ma * mb * mc / (s * s);
                const double pref4 = 0.5 * ma * mb * mc / (s * pair_cache_r(cache, a, b));

                for (int i = 0; i < D; i++) {
                    const double pai = pair_cache_p(cache, a, i);
                    const double pbi = pair_cache_p(cache, b, i);
                    const double pci = pair_cache_p(cache, c, i);
                    const double nab_i = pair_cache_n(cache, a, b, i);
                    const double ui = nab_i + pair_cache_n(cache, a, c, i);
                    const double vi = nab_i + pair_cache_n(cache, c, b, i);

                    velocity_row[a * D + i] += pref3 * (
                        8.0 * ui * vpc * inv_mamc
                      - 16.0 * vi * upc * inv_mamc
                      + 3.0 * ui * vpb * inv_mamb
                      + (ui * vpa + vi * upa) * inv_ma2);
                    velocity_row[b * D + i] += pref3 * (3.0 * vi * upa * inv_mamb);
                    velocity_row[c * D + i] += pref3 * (
                        8.0 * vi * upa * inv_mamc
                      - 16.0 * ui * vpa * inv_mamc
                      + 4.0 * (ui * vpc + vi * upc) * inv_mc2);

                    velocity_row[a * D + i] += pref4 * (
                        8.0 * (pci - nab_i * Gc) * inv_mamc
                      - 3.0 * (pbi - nab_i * Nb) * inv_mamb
                      - 2.0 * (pai - nab_i * Na) * inv_ma2);
                    velocity_row[b * D + i] += pref4 * (
                        -3.0 * (pai - nab_i * Na) * inv_mamb);
                    velocity_row[c * D + i] += pref4 * (
                        8.0 * (pai - nab_i * Na) * inv_mamc
                      - 8.0 * (pci - nab_i * Gc) * inv_mc2);
                }
            }
        }
    }

    for (int ia = 0; ia < na; ia++) {
        const double *const row = buf + (size_t)ia*(size_t)len;
        for (int i = 0; i < len; i++)
            velocity[i] += row[i];
    }
}


// BEGIN GENERATED 2PN TRIPLE MOMENTUM
/*
 * Generated by tools/codegen/generate_2pn_triple_momentum.py from the symbolic
 * Hamiltonian in tools/codegen/pn_expressions.py and the ordered-triple terms in
 * H2PN_cached(). Do not edit this region by hand.
 */

// Add the reducible b = c ordered-triple 2PN contribution to dp/dt.
static void add_momentum_2pn_triple_reducible_pair_analytic(const PairCache *cache,
    double *momentum)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ++ia) {
        const int a = cache->active.ids[ia];
        const double ma = cache->m[a];
        const double inv_ma = cache->inv_m[a];
        const double pa_sq = cache->p2[a];
        const double pa0 = pair_cache_p(cache, a, 0);
        const double pa1 = pair_cache_p(cache, a, 1);
        const double pa2_component = pair_cache_p(cache, a, 2);

        for (int ib = 0; ib < cache->active.num_active; ++ib) {
            const int b = cache->active.ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];
            const double inv_mb = cache->inv_m[b];
            const double pb_sq = cache->p2[b];
            const double papb = pair_cache_p_dot(cache, a, b);
            const double Na = pair_cache_n_dot_p(cache, a, b, a);
            const double Nb = pair_cache_n_dot_p(cache, a, b, b);
            const double inv_r_ab = pair_cache_inv_r(cache, a, b);
            const double nab0 = pair_cache_n(cache, a, b, 0);
            const double nab1 = pair_cache_n(cache, a, b, 1);
            const double nab2 = pair_cache_n(cache, a, b, 2);
            const double pb0 = pair_cache_p(cache, b, 0);
            const double pb1 = pair_cache_p(cache, b, 1);
            const double pb2_component = pair_cache_p(cache, b, 2);

            const double z0 = inv_mb*inv_mb;
            const double z1 = pb_sq*z0;
            const double z2 = z0*(Nb*Nb);
            const double z3 = 2*Nb;
            const double z4 = inv_ma*inv_mb;
            const double z5 = z4*(2*Na + z3);
            const double z6 = 2*(inv_r_ab*inv_r_ab*inv_r_ab);
            const double z7 = 14*Na*z4;
            const double z8 = mb*mb;
            const double z9 = ma*z8;
            const double z10 = (1.0/8.0)*z9*(-z6*(Nb*z5 + 5*z1 - 15*z2) - z6*(-Nb*z7 + 18*pa_sq*(inv_ma*inv_ma) - 50*papb*z4 + 31*z1 + 13*z2));
            const double z11 = inv_r_ab*inv_r_ab;
            const double z12 = (3.0/2.0)*Nb*z11*z4*z9;
            const double z13 = z11*(26*Nb*z0 - z7) + z11*(-30*Nb*z0 + z3*z4 + z5);
            const double z14 = -1.0/8.0*ma*pb0*z13*z8 + pa0*z12;
            const double z15 = -1.0/8.0*ma*pb1*z13*z8 + pa1*z12;
            const double z16 = -1.0/8.0*ma*pb2_component*z13*z8 + pa2_component*z12;
            const double z17 = -nab0*z14 - nab1*z15 - nab2*z16;

            const double force0 = -inv_r_ab*(-nab0*z17 - z14) - nab0*z10;
            const double force1 = -inv_r_ab*(-nab1*z17 - z15) - nab1*z10;
            const double force2 = -inv_r_ab*(-nab2*z17 - z16) - nab2*z10;

            momentum[a * D + 0] += force0;
            momentum[a * D + 1] += force1;
            momentum[a * D + 2] += force2;
            momentum[b * D + 0] -= force0;
            momentum[b * D + 1] -= force1;
            momentum[b * D + 2] -= force2;
        }
    }
}


// Add the complete pairwise-distinct ordered-triple 2PN contribution to dp/dt.
static void add_momentum_2pn_triple_distinct_one(const PairCache *cache, int a, int b,
    int c, double *momentum)
{
    const int D = cache->num_dim;
    const double ma = cache->m[a];
    const double mb = cache->m[b];
    const double mc = cache->m[c];
    const double inv_ma = cache->inv_m[a];
    const double inv_mb = cache->inv_m[b];
    const double inv_mc = cache->inv_m[c];
    const double r_ab = pair_cache_r(cache, a, b);
    const double r_ac = pair_cache_r(cache, a, c);
    const double r_bc = pair_cache_r(cache, b, c);
    const double inv_r_ab = pair_cache_inv_r(cache, a, b);
    const double inv_r_ac = pair_cache_inv_r(cache, a, c);
    const double inv_r_bc = pair_cache_inv_r(cache, b, c);
    const double perimeter = r_ab + r_ac + r_bc;
    const double inv_perimeter = 1.0 / perimeter;
    const double pa_sq = cache->p2[a];
    const double pb_sq = cache->p2[b];
    const double pc_sq = cache->p2[c];
    const double papb = pair_cache_p_dot(cache, a, b);
    const double papc = pair_cache_p_dot(cache, a, c);
    const double pbpc = pair_cache_p_dot(cache, b, c);
    const double Na = pair_cache_n_dot_p(cache, a, b, a);
    const double Nb = pair_cache_n_dot_p(cache, a, b, b);
    const double Gc = pair_cache_n_dot_p(cache, a, b, c);
    const double Ec = pair_cache_n_dot_p(cache, a, c, c);
    const double upa = Na + pair_cache_n_dot_p(cache, a, c, a);
    const double upc = Gc + Ec;
    const double vpa = Na + pair_cache_n_dot_p(cache, c, b, a);
    const double vpb = Nb + pair_cache_n_dot_p(cache, c, b, b);
    const double vpc = Gc + pair_cache_n_dot_p(cache, c, b, c);

    double q = 0.0;
    for (int i = 0; i < D; ++i)
        q += pair_cache_n(cache, a, b, i) * pair_cache_n(cache, a, c, i);

    const double pa0 = pair_cache_p(cache, a, 0);
    const double pa1 = pair_cache_p(cache, a, 1);
    const double pa2_component = pair_cache_p(cache, a, 2);
    const double pb0 = pair_cache_p(cache, b, 0);
    const double pb1 = pair_cache_p(cache, b, 1);
    const double pb2_component = pair_cache_p(cache, b, 2);
    const double pc0 = pair_cache_p(cache, c, 0);
    const double pc1 = pair_cache_p(cache, c, 1);
    const double pc2_component = pair_cache_p(cache, c, 2);
    const double nab0 = pair_cache_n(cache, a, b, 0);
    const double nab1 = pair_cache_n(cache, a, b, 1);
    const double nab2 = pair_cache_n(cache, a, b, 2);
    const double nac0 = pair_cache_n(cache, a, c, 0);
    const double nac1 = pair_cache_n(cache, a, c, 1);
    const double nac2 = pair_cache_n(cache, a, c, 2);
    const double ncb0 = pair_cache_n(cache, c, b, 0);
    const double ncb1 = pair_cache_n(cache, c, b, 1);
    const double ncb2 = pair_cache_n(cache, c, b, 2);

    const double z0 = inv_ma*inv_mb;
    const double z1 = inv_mb*inv_mc;
    const double z2 = inv_ma*inv_ma;
    const double z3 = inv_mb*inv_mb;
    const double z4 = 14*Gc;
    const double z5 = z1*z4;
    const double z6 = 14*z0;
    const double z7 = Na*z6;
    const double z8 = Ec*q;
    const double z9 = z1*z8;
    const double z10 = Nb*z5 - Nb*z7 + Nb*z9 + 18*pa_sq*z2 - 50*papb*z0 + 14*pb_sq*z3 + 17*pbpc*z1 - 2*z3*Nb*Nb;
    const double z11 = mb*mc;
    const double z12 = ma*z11;
    const double z13 = (1.0/8.0)*z12;
    const double z14 = inv_r_ab*inv_r_ab*inv_r_ab;
    const double z15 = ma*ma;
    const double z16 = r_bc*r_bc*r_bc;
    const double z17 = 72*z16;
    const double z18 = r_ab*r_ab*r_ab;
    const double z19 = 56*z18;
    const double z20 = r_bc*r_bc;
    const double z21 = r_ab*z20;
    const double z22 = 60*z21;
    const double z23 = r_ab*r_ab;
    const double z24 = r_ac*r_ac;
    const double z25 = 60*z20;
    const double z26 = r_ab + r_bc;
    const double z27 = 24*z23;
    const double z28 = r_ac*z27;
    const double z29 = -r_ab*z17 + r_ac*z22 + r_bc*z19 + 18*z23*z24 - z23*z25 - z26*z28 + 6*(r_ab*r_ab*r_ab*r_ab) + 35*(r_bc*r_bc*r_bc*r_bc);
    const double z30 = inv_r_ac*inv_r_ac*inv_r_ac;
    const double z31 = (1.0/64.0)*inv_r_bc*z11*z14*z15*z30;
    const double z32 = upa*z2;
    const double z33 = 3*vpb*z0;
    const double z34 = inv_ma*inv_mc;
    const double z35 = 8*z34;
    const double z36 = vpc*z35;
    const double z37 = 16*z34;
    const double z38 = upc*z37;
    const double z39 = inv_mc*inv_mc;
    const double z40 = 4*upc*z39;
    const double z41 = inv_perimeter*inv_perimeter;
    const double z42 = 3*Nb;
    const double z43 = 8*Gc;
    const double z44 = (1.0/2.0)*z12;
    const double z45 = z44*(inv_ma*inv_mc*(-Na*z43 + 8*papc) - z0*(-Na*z42 + 3*papb) - z2*(pa_sq - (Na*Na)) - z39*(4*pc_sq - 4*Gc*Gc));
    const double z46 = inv_r_ab*z41*z45 + z12*(upa*z33 + upa*z36 + vpa*z32 - vpa*z38 + vpc*z40)*(inv_perimeter*inv_perimeter*inv_perimeter);
    const double z47 = -inv_r_ab*z10*z13*inv_r_ac*inv_r_ac + (3.0/64.0)*inv_r_bc*mb*mc*z14*z15*z29*(inv_r_ac*inv_r_ac*inv_r_ac*inv_r_ac) - z31*(36*r_ac*z23 + z22 - z26*z27) - z46;
    const double z48 = inv_r_ab*inv_r_ac;
    const double z49 = Nb*z48;
    const double z50 = z1*z49;
    const double z51 = inv_r_ab*inv_r_ab;
    const double z52 = z39*(Ec*Ec);
    const double z53 = z13*(Ec*z50 + z51*(5*pc_sq*z39 - z52));
    const double z54 = z39*z4;
    const double z55 = 2*Na;
    const double z56 = 2*Nb + z55;
    const double z57 = z13*(q*z50 + z51*(inv_ma*inv_mc*z56 - 2*z39*z8 - z54));
    const double z58 = vpa*z2 + z33 + z36;
    const double z59 = z41*z44;
    const double z60 = pa0*z59;
    const double z61 = -vpa*z37 + 4*vpc*z39;
    const double z62 = pc0*z59;
    const double z63 = z58*z60 + z61*z62;
    const double z64 = nab0*z53 + pc0*z57 + z63;
    const double z65 = pa1*z59;
    const double z66 = pc1*z59;
    const double z67 = z58*z65 + z61*z66;
    const double z68 = nab1*z53 + pc1*z57 + z67;
    const double z69 = pa2_component*z59;
    const double z70 = pc2_component*z59;
    const double z71 = z58*z69 + z61*z70;
    const double z72 = nab2*z53 + pc2_component*z57 + z71;
    const double z73 = nac0*z64 + nac1*z68 + nac2*z72;
    const double z74 = -z28;
    const double z75 = r_bc*z23;
    const double z76 = r_ab*r_ac;
    const double z77 = -inv_perimeter*z45*z51 + (3.0/64.0)*inv_r_bc*mb*mc*z15*z29*z30*(inv_r_ab*inv_r_ab*inv_r_ab*inv_r_ab) + (1.0/8.0)*ma*mb*mc*(-inv_r_ac*z10*z51 - 2*z14*(Ec*inv_ma*inv_mc*z56 - Ec*z54 + 5*pc_sq*q*z39 - q*z52)) - z31*(36*r_ab*z24 + r_ac*z25 - z17 + 24*z18 - 120*z21 - 48*z26*z76 + z74 + 168*z75) - z46;
    const double z78 = Ec*z51;
    const double z79 = 2*z34*z78;
    const double z80 = inv_perimeter*inv_r_ab;
    const double z81 = z44*z80;
    const double z82 = z13*(-z49*z6 + z79) + z81*(z0*z42 + z2*z55 - z34*z43);
    const double z83 = (3.0/2.0)*z0*z12;
    const double z84 = Na*z80*z83 + z13*(z48*(-4*Nb*z3 + z5 - z7 + z9) + z79);
    const double z85 = z13*(14*Nb*inv_mb*inv_mc*inv_r_ab*inv_r_ac - 14*z39*z78) + z81*(-Na*z35 + z39*z43);
    const double z86 = upa*z41*z83;
    const double z87 = z32 - z38;
    const double z88 = upa*z35 + z40;
    const double z89 = pb0*z86 + z60*z87 + z62*z88;
    const double z90 = nac0*z53 + pa0*z82 + pb0*z84 + pc0*z85 + z63 + z89;
    const double z91 = pb1*z86 + z65*z87 + z66*z88;
    const double z92 = nac1*z53 + pa1*z82 + pb1*z84 + pc1*z85 + z67 + z91;
    const double z93 = pb2_component*z86 + z69*z87 + z70*z88;
    const double z94 = nac2*z53 + pa2_component*z82 + pb2_component*z84 + pc2_component*z85 + z71 + z93;
    const double z95 = nab0*z90 + nab1*z92 + nab2*z94;
    const double z96 = inv_r_ab*(-nab0*z95 + z90) + nab0*z77;
    const double z97 = inv_r_ab*(-nab1*z95 + z92) + nab1*z77;
    const double z98 = inv_r_ab*(-nab2*z95 + z94) + nab2*z77;
    const double z99 = (1.0/64.0)*mb*mc*z14*z15*z29*z30*(inv_r_bc*inv_r_bc) - z31*(120*r_bc*z76 + 140*z16 + z19 - 216*z21 + z74 - 120*z75) - z46;
    const double z100 = ncb0*z89 + ncb1*z91 + ncb2*z93;

    const double force_a0 = -inv_r_ac*(-nac0*z73 + z64) - nac0*z47 - z96;
    const double force_a1 = -inv_r_ac*(-nac1*z73 + z68) - nac1*z47 - z97;
    const double force_a2 = -inv_r_ac*(-nac2*z73 + z72) - nac2*z47 - z98;
    const double force_b0 = inv_r_bc*(-ncb0*z100 + z89) + ncb0*z99 + z96;
    const double force_b1 = inv_r_bc*(-ncb1*z100 + z91) + ncb1*z99 + z97;
    const double force_b2 = inv_r_bc*(-ncb2*z100 + z93) + ncb2*z99 + z98;

    momentum[a * D + 0] += force_a0;
    momentum[a * D + 1] += force_a1;
    momentum[a * D + 2] += force_a2;
    momentum[b * D + 0] += force_b0;
    momentum[b * D + 1] += force_b1;
    momentum[b * D + 2] += force_b2;
    momentum[c * D + 0] -= force_a0 + force_b0;
    momentum[c * D + 1] -= force_a1 + force_b1;
    momentum[c * D + 2] -= force_a2 + force_b2;
}
// END GENERATED 2PN TRIPLE MOMENTUM


// Add the analytic combined ordered-triple 2PN Hamiltonian contribution to dp/dt.
void add_momentum_2pn_triple_analytic(const PairCache *cache, double *momentum)
{
    const ActiveList *active = &cache->active;

    if (cache->num_dim != 3)
        errorexit("add_momentum_2pn_triple_analytic currently supports only num_dim = 3");

    add_momentum_2pn_triple_reducible_pair_analytic(cache, momentum);

    const int na = active->num_active;
    const int len = cache->array_half;
    double *const buf = cache->triple_accum;

    memset(buf, 0, (size_t)na*(size_t)len*sizeof(*buf));

#ifdef _OPENMP
    const int threads = cache->triple_accum_threads;
    const int run_parallel = (threads > 1) && !omp_in_parallel();

    #pragma omp parallel for schedule(static) num_threads(threads) if (run_parallel)
#endif
    for (int ia = 0; ia < na; ia++) {
        double *const row = buf + (size_t)ia*(size_t)len;
        const int a = active->ids[ia];

        for (int ib = 0; ib < na; ib++) {
            const int b = active->ids[ib];
            if (b == a)
                continue;

            for (int ic = 0; ic < na; ic++) {
                const int c = active->ids[ic];
                if (c == a || c == b)
                    continue;

                add_momentum_2pn_triple_distinct_one(cache, a, b, c, row);
            }
        }
    }

    for (int ia = 0; ia < na; ia++) {
        const double *const row = buf + (size_t)ia*(size_t)len;
        for (int i = 0; i < len; ++i)
            momentum[i] += row[i];
    }
}


void add_momentum_2pn_fourbody_analytic(const PairCache *cache, double *momentum)
{
    const ActiveList *active = &cache->active;
    const int N = cache->num_bodies;
    const int D = cache->num_dim;

    double Phi[N];
    double Psi[N];

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];
        Phi[a] = 0.0;
        Psi[a] = 0.0;
    }

    // Phi_a = sum_{b != a} m_b / r_ab, compute it using unordered pairs and update both sides
    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        for (int ib = ia + 1; ib < active->num_active; ib++) {
            const int b = active->ids[ib];

            const double inv_r_ab = pair_cache_inv_r(cache, a, b);

            Phi[a] += cache->m[b] * inv_r_ab;
            Phi[b] += cache->m[a] * inv_r_ab;
        }
    }

    // Psi_a = sum_{b != a} m_b Phi_b / r_ab, again use unordered pairs and update both sides
    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        for (int ib = ia + 1; ib < active->num_active; ib++) {
            const int b = active->ids[ib];

            const double inv_r_ab = pair_cache_inv_r(cache, a, b);

            Psi[a] += cache->m[b] * Phi[b] * inv_r_ab;
            Psi[b] += cache->m[a] * Phi[a] * inv_r_ab;
        }
    }

    // Pairwise force accumulation from Q_ab = dH/d(1/r_ab)
    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];

        const double ma = cache->m[a];
        const double Phi_a = Phi[a];
        const double Psi_a = Psi[a];

        for (int ib = ia + 1; ib < active->num_active; ib++) {
            const int b = active->ids[ib];

            const double mb = cache->m[b];
            const double Phi_b = Phi[b];
            const double Psi_b = Psi[b];

            const double Q_star = -0.75 * ma * mb * (Phi_a * Phi_a + Phi_b * Phi_b);
            const double Q_chain = -0.75 * ma * mb * (Phi_a * Phi_b + Psi_a + Psi_b);
            const double Q = Q_star + Q_chain;

            const double inv_r_ab = pair_cache_inv_r(cache, a, b);
            const double coeff = Q * inv_r_ab * inv_r_ab;

            for (int i = 0; i < D; i++) {
                const double force_i = coeff * pair_cache_n(cache, a, b, i);

                momentum[a * D + i] += force_i;
                momentum[b * D + i] -= force_i;
            }
        }
    }
}
