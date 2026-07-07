/**
 * @file eom.c
 * @brief Functions for the right-hand side of the post-Newtonian equations of motion
 */

#include <complex.h>
#include <math.h>
#include "utils.h"
#include "eom.h"
#include "hamiltonian.h"
#include "pair_cache.h"


/**
 * @brief Add the analytic one-body 2PN p^6 Hamiltonian contribution.
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
static void add_eom_2pn_onebody_analytic(const PairCache *cache, double *dwdt)
{
    const int D = cache->num_dim;

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];

        const double inv_m = cache->inv_m[a];
        const double inv_m5 = inv_m * inv_m * inv_m * inv_m * inv_m;
        const double p2 = cache->p2[a];

        const double coeff = 3.0 / 8.0 * p2 * p2 * inv_m5;

        for (int i = 0; i < D; i++)
            dwdt[a * D + i] += coeff * pair_cache_p(cache, a, i);
    }
}


/**
 * @brief Add the analytic ordered-pair 2PN Hamiltonian contribution.
 *
 * This implements the ordered pair terms from H2PN/H2PN_base_complex:
 *
 *   H_ab =
 *     1/(16 r_ab) * F1
 *   + 1/4 * m_a^2 m_b / r_ab^2 * F2
 *   - 1/4 * m_a^2 m_b^2 / r_ab^3,
 *
 * with
 *
 *   F1 =
 *       10 m_a m_b (p_a^2/m_a^2)^2
 *     - 11 (p_a^2/(m_a m_b)) p_b^2
 *     -  2 (p_a.p_b)^2/(m_a m_b)
 *     + 10 (p_a^2/(m_a m_b)) (n_ab.p_b)^2
 *     - 12 (p_a.p_b/(m_a m_b)) (n_ab.p_a)(n_ab.p_b)
 *     -  3 (n_ab.p_a)^2 (n_ab.p_b)^2/(m_a m_b),
 *
 *   F2 =
 *       p_a^2/m_a^2
 *     + p_b^2/m_b^2
 *     - 2 (p_a.p_b)/(m_a m_b).
 *
 * The loop is over ordered pairs a != b, matching the existing Hamiltonian.
 */
static void add_eom_2pn_pair_analytic(const PairCache *cache, double *dwdt)
{
    const int D = cache->num_dim;
    const int P0 = cache->array_half;

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

            const double r = pair_cache_r(cache, a, b);
            const double inv_r = pair_cache_inv_r(cache, a, b);
            const double inv_r2 = cache->inv_r2[pair_cache_idx2(cache, a, b)];
            const double inv_r3 = inv_r2 * inv_r;
            const double inv_r4 = inv_r2 * inv_r2;

            const double pb2 = cache->p2[b];
            const double papb = pair_cache_p_dot(cache, a, b);
            const double Na = pair_cache_n_dot_p(cache, a, b, a);
            const double Nb = pair_cache_n_dot_p(cache, a, b, b);

            /*
             * Scalar pair Hamiltonian pieces.
             *
             * These formulas are intentionally expressed in terms of scalar invariants.
             * The vector chain rule is applied below.
             */
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

            /*
             * Scalar partial derivatives dH_ab/d(invariant).
             *
             * These are generated/checkable with generate_2pn_pair_rhs.py.
             */
            const double dF1_dpa2 =
                20.0 * mb * inv_ma * inv_ma * inv_ma * pa2
              - 11.0 * inv_mamb * pb2
              + 10.0 * inv_mamb * Nb * Nb;

            const double dF1_dpb2 =
              - 11.0 * inv_mamb * pa2;

            const double dF1_dpapb =
              -  4.0 * inv_mamb * papb
              - 12.0 * inv_mamb * Na * Nb;

            const double dF1_dNa =
              - 12.0 * inv_mamb * papb * Nb
              -  6.0 * inv_mamb * Na * Nb * Nb;

            const double dF1_dNb =
                20.0 * inv_mamb * pa2 * Nb
              - 12.0 * inv_mamb * papb * Na
              -  6.0 * inv_mamb * Na * Na * Nb;

            const double dH_dpa2 =
                (1.0 / 16.0) * inv_r * dF1_dpa2
              + 0.25 * mb * inv_r2;

            const double dH_dpb2 =
                (1.0 / 16.0) * inv_r * dF1_dpb2
              + 0.25 * ma2 * inv_mb * inv_r2;

            const double dH_dpapb =
                (1.0 / 16.0) * inv_r * dF1_dpapb
              - 0.5 * ma * inv_r2;

            const double dH_dNa =
                (1.0 / 16.0) * inv_r * dF1_dNa;

            const double dH_dNb =
                (1.0 / 16.0) * inv_r * dF1_dNb;

            const double dH_dr =
              - (1.0 / 16.0) * inv_r2 * F1
              - 0.5 * ma2 * mb * inv_r3 * F2
              + 0.75 * ma2 * mb2 * inv_r4;

            /*
             * Vector chain rule.
             *
             * Momentum derivatives:
             *
             *   dH/dp_a^i = 2 dH/dpa2 p_a^i + dH/dpapb p_b^i + dH/dNa n_ab^i
             *   dH/dp_b^i = 2 dH/dpb2 p_b^i + dH/dpapb p_a^i + dH/dNb n_ab^i
             *
             * Position derivatives:
             *
             *   dH/dx_a^i =
             *       dH/dr n_i
             *     + dH/dNa * (p_a^i - Na n_i)/r
             *     + dH/dNb * (p_b^i - Nb n_i)/r
             *
             * and dH/dx_b^i = -dH/dx_a^i.
             */
            for (int i = 0; i < D; i++) {
                const double pai = pair_cache_p(cache, a, i);
                const double pbi = pair_cache_p(cache, b, i);
                const double ni = pair_cache_n(cache, a, b, i);

                const double dH_dpa_i =
                    2.0 * dH_dpa2 * pai
                  + dH_dpapb * pbi
                  + dH_dNa * ni;

                const double dH_dpb_i =
                    2.0 * dH_dpb2 * pbi
                  + dH_dpapb * pai
                  + dH_dNb * ni;

                const double dNa_dxa_i = inv_r * (pai - Na * ni);
                const double dNb_dxa_i = inv_r * (pbi - Nb * ni);

                const double dH_dxa_i =
                    dH_dr * ni
                  + dH_dNa * dNa_dxa_i
                  + dH_dNb * dNb_dxa_i;

                // dx/dt = dH/dp
                dwdt[a * D + i] += dH_dpa_i;
                dwdt[b * D + i] += dH_dpb_i;

                // dp/dt = -dH/dx
                dwdt[P0 + a * D + i] -= dH_dxa_i;
                dwdt[P0 + b * D + i] += dH_dxa_i;
            }
        }
    }
}



/**
 * @brief Add the analytic ordered-triple 2PN Hamiltonian contribution.
 *
 * This function assumes num_dim == 3. It implements the ordered-triple terms from H_2PN
 */
static void add_eom_2pn_triple_analytic(const PairCache *cache, const double *w, double *dwdt)
{
    const int D = cache->num_dim;
    const int P0 = cache->array_half;

    if (D != 3)
        errorexit("add_eom_2pn_triple_analytic currently supports only num_dim = 3");

    for (int ia = 0; ia < cache->active.num_active; ia++) {
        const int a = cache->active.ids[ia];
        const double ma = cache->m[a];

        const double xa0 = w[a * D + 0];
        const double xa1 = w[a * D + 1];
        const double xa2 = w[a * D + 2];

        const double pa0 = pair_cache_p(cache, a, 0);
        const double pa1 = pair_cache_p(cache, a, 1);
        const double pa2 = pair_cache_p(cache, a, 2);

        for (int ib = 0; ib < cache->active.num_active; ib++) {
            const int b = cache->active.ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];

            const double xb0 = w[b * D + 0];
            const double xb1 = w[b * D + 1];
            const double xb2 = w[b * D + 2];

            const double pb0 = pair_cache_p(cache, b, 0);
            const double pb1 = pair_cache_p(cache, b, 1);
            const double pb2 = pair_cache_p(cache, b, 2);

            for (int ic = 0; ic < cache->active.num_active; ic++) {
                const int c = cache->active.ids[ic];
                if (c == a)
                    continue;

                const double mc = cache->m[c];

                const double xc0 = w[c * D + 0];
                const double xc1 = w[c * D + 1];
                const double xc2 = w[c * D + 2];

                const double pc0 = pair_cache_p(cache, c, 0);
                const double pc1 = pair_cache_p(cache, c, 1);
                const double pc2 = pair_cache_p(cache, c, 2);

                /*
                 * Block with condition b != a && c != a.
                 * This also includes the case c == b, matching the original Hamiltonian.
                 */
                {

            const double t0 = xa0 - xb0;
            const double t1 = pow(t0, 2);
            const double t2 = xa1 - xb1;
            const double t3 = pow(t2, 2);
            const double t4 = xa2 - xb2;
            const double t5 = pow(t4, 2);
            const double t6 = t1 + t3 + t5;
            const double t7 = pow(t6, -3.0/2.0);
            const double t8 = -t0*t7;
            const double t9 = xa0 - xc0;
            const double t10 = pow(t9, 2);
            const double t11 = xa1 - xc1;
            const double t12 = pow(t11, 2);
            const double t13 = xa2 - xc2;
            const double t14 = pow(t13, 2);
            const double t15 = t10 + t12 + t14;
            const double t16 = pow(t15, -1.0/2.0);
            const double t17 = pow(ma, -2);
            const double t18 = 1.0/mb;
            const double t19 = 1.0/ma;
            const double t20 = t18*t19;
            const double t21 = 50*t20;
            const double t22 = pow(mb, -2);
            const double t23 = 17*pc0;
            const double t24 = 17*pc1;
            const double t25 = 17*pc2;
            const double t26 = 1.0/mc;
            const double t27 = t18*t26;
            const double t28 = pow(t6, -1.0/2.0);
            const double t29 = pb0*t28;
            const double t30 = t0*t29;
            const double t31 = pb1*t28;
            const double t32 = t2*t31;
            const double t33 = pb2*t28;
            const double t34 = t33*t4;
            const double t35 = t30 + t32 + t34;
            const double t36 = 2*t22;
            const double t37 = 14*t30 + 14*t32 + 14*t34;
            const double t38 = pc0*t28;
            const double t39 = pc1*t28;
            const double t40 = t2*t39;
            const double t41 = pc2*t28;
            const double t42 = t4*t41;
            const double t43 = t27*(t0*t38 + t40 + t42);
            const double t44 = pa0*t28;
            const double t45 = 14*t44;
            const double t46 = pa1*t28;
            const double t47 = 14*t46;
            const double t48 = pa2*t28;
            const double t49 = 14*t48;
            const double t50 = t0*t45 + t2*t47 + t4*t49;
            const double t51 = t20*t35;
            const double t52 = t16*t28;
            const double t53 = t0*t52;
            const double t54 = t2*t52;
            const double t55 = t11*t54;
            const double t56 = t4*t52;
            const double t57 = t13*t56;
            const double t58 = t53*t9 + t55 + t57;
            const double t59 = pc0*t16;
            const double t60 = pc1*t16;
            const double t61 = pc2*t16;
            const double t62 = t11*t60 + t13*t61 + t59*t9;
            const double t63 = t27*t35;
            const double t64 = t62*t63;
            const double t65 = 18*t17*(pow(pa0, 2) + pow(pa1, 2) + pow(pa2, 2)) - t21*(pa0*pb0 + pa1*pb1 + pa2*pb2) + t22*(14*pow(pb0, 2) + 14*pow(pb1, 2) + 14*pow(pb2, 2)) + t27*(pb0*t23 + pb1*t24 + pb2*t25) - pow(t35, 2)*t36 + t37*t43 - t50*t51 + t58*t64;
            const double t66 = (1.0/8.0)*ma*mc;
            const double t67 = mb*t66;
            const double t68 = t65*t67;
            const double t69 = t16*t68;
            const double t70 = -t9;
            const double t71 = pow(t15, -3.0/2.0);
            const double t72 = t28*t71;
            const double t73 = t68*t72;
            const double t74 = 2*xa0 - 2*xb0;
            const double t75 = 5*t53;
            const double t76 = 5*t55 + 5*t57 + t75*t9;
            const double t77 = pow(mc, -2);
            const double t78 = t77*(pow(pc0, 2) + pow(pc1, 2) + pow(pc2, 2));
            const double t79 = 14*t38;
            const double t80 = t0*t79 + 14*t40 + 14*t42;
            const double t81 = t62*t77;
            const double t82 = 2*t44;
            const double t83 = 2*t46;
            const double t84 = 2*t48;
            const double t85 = t0*t82 + t2*t83 + t4*t84;
            const double t86 = t19*t26;
            const double t87 = t62*t86;
            const double t88 = 2*t29;
            const double t89 = t0*t88 + 2*t32 + 2*t34;
            const double t90 = pow(t62, 2)*t77;
            const double t91 = t67*(-t58*t90 + t76*t78 - t80*t81 + t85*t87 + t87*t89)/pow(t6, 2);
            const double t92 = t0*t8;
            const double t93 = pb0*t92;
            const double t94 = t2*t8;
            const double t95 = pb1*t94;
            const double t96 = t4*t8;
            const double t97 = pb2*t96;
            const double t98 = t88 + 2*t93 + 2*t95 + 2*t97;
            const double t99 = t35*t36;
            const double t100 = pc0*t92;
            const double t101 = pc1*t94;
            const double t102 = pc2*t96;
            const double t103 = t27*t37;
            const double t104 = t29 + t93 + t95 + t97;
            const double t105 = t20*t50;
            const double t106 = 14*t29;
            const double t107 = pa0*t92;
            const double t108 = pa1*t94;
            const double t109 = pa2*t96;
            const double t110 = pc0*t71;
            const double t111 = t110*t9;
            const double t112 = t111*t70;
            const double t113 = pc1*t71;
            const double t114 = t11*t113;
            const double t115 = t114*t70;
            const double t116 = pc2*t71;
            const double t117 = t116*t13;
            const double t118 = t117*t70;
            const double t119 = t112 + t115 + t118 + t59;
            const double t120 = t58*t63;
            const double t121 = t27*t58*t62;
            const double t122 = t52*t9;
            const double t123 = t16*t9;
            const double t124 = t123*t92;
            const double t125 = t11*t16;
            const double t126 = t125*t94;
            const double t127 = t13*t16;
            const double t128 = t127*t96;
            const double t129 = t0*t72;
            const double t130 = t129*t9;
            const double t131 = t130*t70;
            const double t132 = t2*t72;
            const double t133 = t11*t132;
            const double t134 = t133*t70;
            const double t135 = t4*t72;
            const double t136 = t13*t135;
            const double t137 = t136*t70;
            const double t138 = t122 + t124 + t126 + t128 + t131 + t134 + t137 + t53;
            const double t139 = t52*t67;
            const double t140 = t77*t80;
            const double t141 = t119*t86;
            const double t142 = 2*t59;
            const double t143 = t58*t81;
            const double t144 = 5*t122;
            const double t145 = t67/t6;
            const double t146 = -t2*t7;
            const double t147 = -t11;
            const double t148 = 2*xa1 - 2*xb1;
            const double t149 = 2*t31;
            const double t150 = t0*t146;
            const double t151 = pb0*t150;
            const double t152 = t146*t2;
            const double t153 = pb1*t152;
            const double t154 = t146*t4;
            const double t155 = pb2*t154;
            const double t156 = t149 + 2*t151 + 2*t153 + 2*t155;
            const double t157 = pc0*t150;
            const double t158 = pc1*t152;
            const double t159 = pc2*t154;
            const double t160 = t151 + t153 + t155 + t31;
            const double t161 = 14*t31;
            const double t162 = pa0*t150;
            const double t163 = pa1*t152;
            const double t164 = pa2*t154;
            const double t165 = t111*t147;
            const double t166 = t114*t147;
            const double t167 = t117*t147;
            const double t168 = t165 + t166 + t167 + t60;
            const double t169 = t11*t52;
            const double t170 = t123*t150;
            const double t171 = t130*t147;
            const double t172 = t125*t152;
            const double t173 = t127*t154;
            const double t174 = t133*t147;
            const double t175 = t136*t147;
            const double t176 = t169 + t170 + t171 + t172 + t173 + t174 + t175 + t54;
            const double t177 = 14*t39;
            const double t178 = t168*t86;
            const double t179 = 2*t60;
            const double t180 = 5*t54;
            const double t181 = 5*t169;
            const double t182 = -t4*t7;
            const double t183 = -t13;
            const double t184 = 2*xa2 - 2*xb2;
            const double t185 = 2*t33;
            const double t186 = t0*t182;
            const double t187 = pb0*t186;
            const double t188 = t182*t2;
            const double t189 = pb1*t188;
            const double t190 = t182*t4;
            const double t191 = pb2*t190;
            const double t192 = t185 + 2*t187 + 2*t189 + 2*t191;
            const double t193 = pc0*t186;
            const double t194 = pc1*t188;
            const double t195 = pc2*t190;
            const double t196 = t187 + t189 + t191 + t33;
            const double t197 = 14*t33;
            const double t198 = pa0*t186;
            const double t199 = pa1*t188;
            const double t200 = pa2*t190;
            const double t201 = t111*t183;
            const double t202 = t114*t183;
            const double t203 = t117*t183;
            const double t204 = t201 + t202 + t203 + t61;
            const double t205 = t13*t52;
            const double t206 = t123*t186;
            const double t207 = t130*t183;
            const double t208 = t125*t188;
            const double t209 = t133*t183;
            const double t210 = t127*t190;
            const double t211 = t136*t183;
            const double t212 = t205 + t206 + t207 + t208 + t209 + t210 + t211 + t56;
            const double t213 = 14*t41;
            const double t214 = t204*t86;
            const double t215 = 2*t61;
            const double t216 = 5*t56;
            const double t217 = 5*t205;
            const double t218 = t0*t7;
            const double t219 = mb*t218;
            const double t220 = t1*t7;
            const double t221 = pc0*t220;
            const double t222 = t2*t218;
            const double t223 = pc1*t222;
            const double t224 = t218*t4;
            const double t225 = pc2*t224;
            const double t226 = pa0*t220;
            const double t227 = pa1*t222;
            const double t228 = pa2*t224;
            const double t229 = pb0*t220;
            const double t230 = pb1*t222;
            const double t231 = pb2*t224;
            const double t232 = 2*t229 + 2*t230 + 2*t231 - t88;
            const double t233 = t123*t220;
            const double t234 = t125*t222;
            const double t235 = t127*t224;
            const double t236 = -t122 + t233 + t234 + t235;
            const double t237 = t229 + t230 + t231 - t29;
            const double t238 = t2*t7;
            const double t239 = t3*t7;
            const double t240 = pc1*t239;
            const double t241 = pc0*t222;
            const double t242 = t238*t4;
            const double t243 = pc2*t242;
            const double t244 = pa1*t239;
            const double t245 = pa0*t222;
            const double t246 = pa2*t242;
            const double t247 = pb1*t239;
            const double t248 = pb0*t222;
            const double t249 = pb2*t242;
            const double t250 = -t149 + 2*t247 + 2*t248 + 2*t249;
            const double t251 = t125*t239;
            const double t252 = t123*t222;
            const double t253 = t127*t242;
            const double t254 = -t169 + t251 + t252 + t253;
            const double t255 = t247 + t248 + t249 - t31;
            const double t256 = t4*t7;
            const double t257 = t5*t7;
            const double t258 = pc2*t257;
            const double t259 = pc0*t224;
            const double t260 = pc1*t242;
            const double t261 = pa2*t257;
            const double t262 = pa0*t224;
            const double t263 = pa1*t242;
            const double t264 = pb2*t257;
            const double t265 = pb0*t224;
            const double t266 = pb1*t242;
            const double t267 = -t185 + 2*t264 + 2*t265 + 2*t266;
            const double t268 = t127*t257;
            const double t269 = t123*t224;
            const double t270 = t125*t242;
            const double t271 = -t205 + t268 + t269 + t270;
            const double t272 = t264 + t265 + t266 - t33;
            const double t273 = t10*t110;
            const double t274 = t114*t9;
            const double t275 = t117*t9;
            const double t276 = t273 + t274 + t275 - t59;
            const double t277 = t10*t129;
            const double t278 = t133*t9;
            const double t279 = t136*t9;
            const double t280 = t277 + t278 + t279 - t53;
            const double t281 = t276*t86;
            const double t282 = t113*t12;
            const double t283 = t11*t111;
            const double t284 = t11*t117;
            const double t285 = t282 + t283 + t284 - t60;
            const double t286 = t12*t132;
            const double t287 = t11*t130;
            const double t288 = t11*t136;
            const double t289 = t286 + t287 + t288 - t54;
            const double t290 = t285*t86;
            const double t291 = t116*t14;
            const double t292 = t111*t13;
            const double t293 = t114*t13;
            const double t294 = t291 + t292 + t293 - t61;
            const double t295 = t135*t14;
            const double t296 = t13*t130;
            const double t297 = t13*t133;
            const double t298 = t295 + t296 + t297 - t56;
            const double t299 = t294*t86;
            const double t300 = (1.0/4.0)*t62;
            const double t301 = t219*t300;
            const double t302 = t0*t28;
            const double t303 = 14*t51;
            const double t304 = mb*t300;
            const double t305 = t238*t304;
            const double t306 = t2*t28;
            const double t307 = t256*t304;
            const double t308 = t28*t4;
            const double t309 = 28*t22;
            const double t310 = 4*t22*t35;
            const double t311 = 14*t43;
            const double t312 = 17*t27;
            const double t313 = 14*t81;
            const double t314 = t123*t86;
            const double t315 = 2*t76*t77;
            const double t316 = 2*t143;
            const double t317 = t125*t86;
            const double t318 = t127*t86;
            dwdt[P0 + a * D + 0] += -t139*(t103*(t100 + t101 + t102 + t38) - t104*t105 + t104*t121 + t119*t120 + t138*t64 + t43*(t106 + 14*t93 + 14*t95 + 14*t97) - t51*(14*t107 + 14*t108 + 14*t109 + t45) - t98*t99) - t145*(-t119*t140 - t138*t90 + t141*t85 + t141*t89 - t143*(2*t112 + 2*t115 + 2*t118 + t142) + t78*(5*t124 + 5*t126 + 5*t128 + 5*t131 + 5*t134 + 5*t137 + t144 + t75) - t81*(14*t100 + 14*t101 + 14*t102 + t79) + t87*t98 + t87*(2*t107 + 2*t108 + 2*t109 + t82)) - t69*t8 - t70*t73 + t74*t91;
            dwdt[P0 + a * D + 1] += -t139*(t103*(t157 + t158 + t159 + t39) - t105*t160 + t120*t168 + t121*t160 - t156*t99 + t176*t64 + t43*(14*t151 + 14*t153 + 14*t155 + t161) - t51*(14*t162 + 14*t163 + 14*t164 + t47)) - t145*(-t140*t168 - t143*(2*t165 + 2*t166 + 2*t167 + t179) + t156*t87 - t176*t90 + t178*t85 + t178*t89 + t78*(5*t170 + 5*t171 + 5*t172 + 5*t173 + 5*t174 + 5*t175 + t180 + t181) - t81*(14*t157 + 14*t158 + 14*t159 + t177) + t87*(2*t162 + 2*t163 + 2*t164 + t83)) - t146*t69 - t147*t73 + t148*t91;
            dwdt[P0 + a * D + 2] += -t139*(t103*(t193 + t194 + t195 + t41) - t105*t196 + t120*t204 + t121*t196 - t192*t99 + t212*t64 + t43*(14*t187 + 14*t189 + 14*t191 + t197) - t51*(14*t198 + 14*t199 + 14*t200 + t49)) - t145*(-t140*t204 - t143*(2*t201 + 2*t202 + 2*t203 + t215) + t192*t87 - t212*t90 + t214*t85 + t214*t89 + t78*(5*t206 + 5*t207 + 5*t208 + 5*t209 + 5*t210 + 5*t211 + t216 + t217) - t81*(14*t193 + 14*t194 + 14*t195 + t213) + t87*(2*t198 + 2*t199 + 2*t200 + t84)) - t182*t69 - t183*t73 + t184*t91;
            dwdt[P0 + b * D + 0] += -t139*(t103*(t221 + t223 + t225 - t38) - t105*t237 + t121*t237 - t232*t99 + t236*t64 + t43*(-t106 + 14*t229 + 14*t230 + 14*t231) - t51*(14*t226 + 14*t227 + 14*t228 - t45)) - t145*(t232*t87 - t236*t90 + t78*(-t144 + 5*t233 + 5*t234 + 5*t235) - t81*(14*t221 + 14*t223 + 14*t225 - t79) + t87*(2*t226 + 2*t227 + 2*t228 - t82)) - t16*t219*t65*t66 - t74*t91;
            dwdt[P0 + b * D + 1] += -t139*(t103*(t240 + t241 + t243 - t39) - t105*t255 + t121*t255 - t250*t99 + t254*t64 + t43*(-t161 + 14*t247 + 14*t248 + 14*t249) - t51*(14*t244 + 14*t245 + 14*t246 - t47)) - t145*(t250*t87 - t254*t90 + t78*(-t181 + 5*t251 + 5*t252 + 5*t253) - t81*(-t177 + 14*t240 + 14*t241 + 14*t243) + t87*(2*t244 + 2*t245 + 2*t246 - t83)) - t148*t91 - t238*t69;
            dwdt[P0 + b * D + 2] += -t139*(t103*(t258 + t259 + t260 - t41) - t105*t272 + t121*t272 - t267*t99 + t271*t64 + t43*(-t197 + 14*t264 + 14*t265 + 14*t266) - t51*(14*t261 + 14*t262 + 14*t263 - t49)) - t145*(t267*t87 - t271*t90 + t78*(-t217 + 5*t268 + 5*t269 + 5*t270) - t81*(-t213 + 14*t258 + 14*t259 + 14*t260) + t87*(2*t261 + 2*t262 + 2*t263 - t84)) - t184*t91 - t256*t69;
            dwdt[P0 + c * D + 0] += -t139*(t120*t276 + t280*t64) - t145*(-t140*t276 - t143*(-t142 + 2*t273 + 2*t274 + 2*t275) - t280*t90 + t281*t85 + t281*t89 + t78*(5*t277 + 5*t278 + 5*t279 - t75)) - t73*t9;
            dwdt[P0 + c * D + 1] += -t11*t73 - t139*(t120*t285 + t289*t64) - t145*(-t140*t285 - t143*(-t179 + 2*t282 + 2*t283 + 2*t284) - t289*t90 + t290*t85 + t290*t89 + t78*(-t180 + 5*t286 + 5*t287 + 5*t288));
            dwdt[P0 + c * D + 2] += -t13*t73 - t139*(t120*t294 + t298*t64) - t145*(-t140*t294 - t143*(-t215 + 2*t291 + 2*t292 + 2*t293) - t298*t90 + t299*t85 + t299*t89 + t78*(-t216 + 5*t295 + 5*t296 + 5*t297));
            dwdt[a * D + 0] += t139*(36*pa0*t17 - pb0*t21 - t302*t303) + t301;
            dwdt[a * D + 1] += t139*(36*pa1*t17 - pb1*t21 - t303*t306) + t305;
            dwdt[a * D + 2] += t139*(36*pa2*t17 - pb2*t21 - t303*t308) + t307;
            dwdt[b * D + 0] += t139*(-pa0*t21 + pb0*t309 - t105*t302 + t121*t302 + t23*t27 - t302*t310 + t302*t311) + t301;
            dwdt[b * D + 1] += t139*(-pa1*t21 + pb1*t309 - t105*t306 + t121*t306 + t24*t27 - t306*t310 + t306*t311) + t305;
            dwdt[b * D + 2] += t139*(-pa2*t21 + pb2*t309 - t105*t308 + t121*t308 + t25*t27 - t308*t310 + t308*t311) + t307;
            dwdt[c * D + 0] += t139*(pb0*t312 + t103*t302 + t120*t123) + t145*(pc0*t315 - t123*t140 - t123*t316 - t302*t313 + t314*t85 + t314*t89);
            dwdt[c * D + 1] += t139*(pb1*t312 + t103*t306 + t120*t125) + t145*(pc1*t315 - t125*t140 - t125*t316 - t306*t313 + t317*t85 + t317*t89);
            dwdt[c * D + 2] += t139*(pb2*t312 + t103*t308 + t120*t127) + t145*(pc2*t315 - t127*t140 - t127*t316 - t308*t313 + t318*t85 + t318*t89);

                }

                /*
                 * Block with condition b != a && c != a && c != b.
                 */
                if (c != b) {

            const double u0 = pow(ma, 2);
            const double u1 = 3*xa0;
            const double u2 = 3*xb0;
            const double u3 = u1 - u2;
            const double u4 = -u3;
            const double u5 = xa0 - xb0;
            const double u6 = pow(u5, 2);
            const double u7 = xa1 - xb1;
            const double u8 = pow(u7, 2);
            const double u9 = xa2 - xb2;
            const double u10 = pow(u9, 2);
            const double u11 = u10 + u6 + u8;
            const double u12 = pow(u11, -5.0/2.0);
            const double u13 = -xc0;
            const double u14 = u13 + xa0;
            const double u15 = pow(u14, 2);
            const double u16 = -xc1;
            const double u17 = u16 + xa1;
            const double u18 = pow(u17, 2);
            const double u19 = -xc2;
            const double u20 = u19 + xa2;
            const double u21 = pow(u20, 2);
            const double u22 = u15 + u18 + u21;
            const double u23 = pow(u22, -3.0/2.0);
            const double u24 = u13 + xb0;
            const double u25 = -u24;
            const double u26 = pow(u25, 2);
            const double u27 = u16 + xb1;
            const double u28 = -u27;
            const double u29 = pow(u28, 2);
            const double u30 = u19 + xb2;
            const double u31 = -u30;
            const double u32 = pow(u31, 2);
            const double u33 = u26 + u29 + u32;
            const double u34 = sqrt(u33);
            const double u35 = 1.0/u34;
            const double u36 = sqrt(u11);
            const double u37 = pow(u33, 3.0/2.0);
            const double u38 = 72*u37;
            const double u39 = pow(u11, 3.0/2.0);
            const double u40 = 56*u34;
            const double u41 = 18*u10 + 18*u6 + 18*u8;
            const double u42 = 60*u10 + 60*u6 + 60*u8;
            const double u43 = sqrt(u22);
            const double u44 = 60*u36*u43;
            const double u45 = 24*u10 + 24*u6 + 24*u8;
            const double u46 = u34 + u36;
            const double u47 = u43*u46;
            const double u48 = 6*pow(u11, 2) + u22*u41 + 35*pow(u33, 2) - u33*u42 + u33*u44 - u36*u38 + u39*u40 - u45*u47;
            const double u49 = -3*xc0;
            const double u50 = u1 + u49;
            const double u51 = 1.0/u39;
            const double u52 = pow(u22, -5.0/2.0);
            const double u53 = 36*xa0 - 36*xb0;
            const double u54 = 4*xb0;
            const double u55 = -u54 + 4*xa0;
            const double u56 = 6*u11;
            const double u57 = 120*xa0 - 120*xb0;
            const double u58 = -2*xc0;
            const double u59 = u58 + 2*xa0;
            const double u60 = 1.0/u36;
            const double u61 = u5*u60;
            const double u62 = u36*u40;
            const double u63 = 60*u33;
            const double u64 = u43*u63;
            const double u65 = 1.0/u43;
            const double u66 = u14*u65;
            const double u67 = u36*u63;
            const double u68 = 48*xa0 - 48*xb0;
            const double u69 = u43*u45;
            const double u70 = u45*u46;
            const double u71 = -u5;
            const double u72 = u51*u71;
            const double u73 = 1.0/ma;
            const double u74 = 1.0/mb;
            const double u75 = u73*u74;
            const double u76 = 3*u75;
            const double u77 = 1.0/u0;
            const double u78 = u60*u7;
            const double u79 = u60*u9;
            const double u80 = pa0*u61 + pa1*u78 + pa2*u79;
            const double u81 = pow(mc, -2);
            const double u82 = pc0*u61 + pc1*u78 + pc2*u79;
            const double u83 = pb0*u61 + pb1*u78 + pb2*u79;
            const double u84 = 1.0/mc;
            const double u85 = 8*pc0;
            const double u86 = 8*pc1;
            const double u87 = 8*pc2;
            const double u88 = 8*u80;
            const double u89 = 3*u73*u74*u80*u83 + u73*u84*(pa0*u85 + pa1*u86 + pa2*u87 - u82*u88) - u76*(pa0*pb0 + pa1*pb1 + pa2*pb2) - u77*(pow(pa0, 2) + pow(pa1, 2) + pow(pa2, 2) - pow(u80, 2)) - u81*(4*pow(pc0, 2) + 4*pow(pc1, 2) + 4*pow(pc2, 2) - 4*pow(u82, 2));
            const double u90 = u43 + u46;
            const double u91 = (1.0/2.0)*mc;
            const double u92 = ma*mb;
            const double u93 = u91*u92;
            const double u94 = u93/u90;
            const double u95 = u89*u94;
            const double u96 = u61 + u66;
            const double u97 = u91/pow(u90, 2);
            const double u98 = u92*u97;
            const double u99 = u60*u89*u98;
            const double u100 = pa0*u60;
            const double u101 = 2*u100;
            const double u102 = 2*pa0;
            const double u103 = u5*u72;
            const double u104 = 2*pa1;
            const double u105 = u7*u72;
            const double u106 = 2*pa2;
            const double u107 = u72*u9;
            const double u108 = u77*u80;
            const double u109 = pb0*u60;
            const double u110 = pb1*u105 + pb2*u107;
            const double u111 = u76*u80;
            const double u112 = pa1*u105 + pa2*u107;
            const double u113 = pa0*u103 + u100 + u112;
            const double u114 = u76*u83;
            const double u115 = pc0*u60;
            const double u116 = 2*u115;
            const double u117 = pc0*u103;
            const double u118 = pc1*u105;
            const double u119 = pc2*u107;
            const double u120 = 4*u81*u82;
            const double u121 = u118 + u119;
            const double u122 = 8*u82;
            const double u123 = u73*u84;
            const double u124 = u60*u94;
            const double u125 = 2*u61;
            const double u126 = u25*u35;
            const double u127 = u126 + u61;
            const double u128 = u28*u35;
            const double u129 = u128 + u78;
            const double u130 = u31*u35;
            const double u131 = u130 + u79;
            const double u132 = pa0*u127 + pa1*u129 + pa2*u131;
            const double u133 = pa0*u96;
            const double u134 = u17*u65;
            const double u135 = u134 + u78;
            const double u136 = pa1*u135;
            const double u137 = u20*u65;
            const double u138 = u137 + u79;
            const double u139 = pa2*u138;
            const double u140 = u77*(u133 + u136 + u139);
            const double u141 = pc0*u127 + pc1*u129 + pc2*u131;
            const double u142 = 4*pc0;
            const double u143 = 4*pc1;
            const double u144 = 4*pc2;
            const double u145 = u81*(u135*u143 + u138*u144 + u142*u96);
            const double u146 = pb0*u127 + pb1*u129 + pb2*u131;
            const double u147 = 3*u133 + 3*u136 + 3*u139;
            const double u148 = u147*u75;
            const double u149 = u123*(8*u133 + 8*u136 + 8*u139);
            const double u150 = 16*pc0;
            const double u151 = 16*pc1;
            const double u152 = 16*pc2;
            const double u153 = u123*(u135*u151 + u138*u152 + u150*u96);
            const double u154 = u93*(u132*u140 - u132*u153 + u141*u145 + u141*u149 + u146*u148)/pow(u90, 3);
            const double u155 = u103 + u60;
            const double u156 = pa0*u155 + u112;
            const double u157 = pc0*u155 + u121;
            const double u158 = -u14;
            const double u159 = u158*u23;
            const double u160 = u105 + u159*u17;
            const double u161 = pa1*u160;
            const double u162 = u107 + u159*u20;
            const double u163 = pa2*u162;
            const double u164 = u14*u159 + u155 + u65;
            const double u165 = pa0*u164;
            const double u166 = u132*u77;
            const double u167 = u141*u81;
            const double u168 = u146*u75;
            const double u169 = u123*u141;
            const double u170 = u123*u132;
            const double u171 = 3*xa1;
            const double u172 = 3*xb1;
            const double u173 = u171 - u172;
            const double u174 = -u173;
            const double u175 = -3*xc1;
            const double u176 = u171 + u175;
            const double u177 = 36*xa1 - 36*xb1;
            const double u178 = 4*xb1;
            const double u179 = -u178 + 4*xa1;
            const double u180 = 120*xa1 - 120*xb1;
            const double u181 = -2*xc1;
            const double u182 = u181 + 2*xa1;
            const double u183 = 48*xa1 - 48*xb1;
            const double u184 = -u7;
            const double u185 = u184*u51;
            const double u186 = pa1*u60;
            const double u187 = 2*u186;
            const double u188 = u185*u5;
            const double u189 = u185*u7;
            const double u190 = u185*u9;
            const double u191 = pb1*u60;
            const double u192 = pb0*u188 + pb2*u190;
            const double u193 = pa0*u188 + pa2*u190;
            const double u194 = pa1*u189 + u186 + u193;
            const double u195 = pc1*u60;
            const double u196 = 2*u195;
            const double u197 = pc0*u188;
            const double u198 = pc1*u189;
            const double u199 = pc2*u190;
            const double u200 = u197 + u199;
            const double u201 = 2*u78;
            const double u202 = u189 + u60;
            const double u203 = pa1*u202 + u193;
            const double u204 = pc1*u202 + u200;
            const double u205 = -u17;
            const double u206 = u205*u23;
            const double u207 = u14*u206 + u188;
            const double u208 = pa0*u207;
            const double u209 = u190 + u20*u206;
            const double u210 = pa2*u209;
            const double u211 = u17*u206 + u202 + u65;
            const double u212 = pa1*u211;
            const double u213 = 3*xa2;
            const double u214 = 3*xb2;
            const double u215 = u213 - u214;
            const double u216 = -u215;
            const double u217 = -3*xc2;
            const double u218 = u213 + u217;
            const double u219 = 36*xa2 - 36*xb2;
            const double u220 = 4*xb2;
            const double u221 = -u220 + 4*xa2;
            const double u222 = 120*xa2 - 120*xb2;
            const double u223 = -2*xc2;
            const double u224 = u223 + 2*xa2;
            const double u225 = 48*xa2 - 48*xb2;
            const double u226 = -u9;
            const double u227 = u226*u51;
            const double u228 = pa2*u60;
            const double u229 = 2*u228;
            const double u230 = u227*u5;
            const double u231 = u227*u7;
            const double u232 = u227*u9;
            const double u233 = pb2*u60;
            const double u234 = pb0*u230 + pb1*u231;
            const double u235 = pa0*u230 + pa1*u231;
            const double u236 = pa2*u232 + u228 + u235;
            const double u237 = pc2*u60;
            const double u238 = 2*u237;
            const double u239 = pc0*u230;
            const double u240 = pc1*u231;
            const double u241 = pc2*u232;
            const double u242 = u239 + u240;
            const double u243 = 2*u79;
            const double u244 = u232 + u60;
            const double u245 = pa2*u244 + u235;
            const double u246 = pc2*u244 + u242;
            const double u247 = -u20;
            const double u248 = u23*u247;
            const double u249 = u14*u248 + u230;
            const double u250 = pa0*u249;
            const double u251 = u17*u248 + u231;
            const double u252 = pa1*u251;
            const double u253 = u20*u248 + u244 + u65;
            const double u254 = pa2*u253;
            const double u255 = 1.0/u37;
            const double u256 = u5*u51;
            const double u257 = u54 - 4*xc0;
            const double u258 = 35*u33;
            const double u259 = u58 + 2*xb0;
            const double u260 = u60*u71;
            const double u261 = u24*u35;
            const double u262 = 56*u39;
            const double u263 = u2 + u49;
            const double u264 = 72*u34*u36;
            const double u265 = u260 + u261;
            const double u266 = u51*u6;
            const double u267 = u256*u7;
            const double u268 = u256*u9;
            const double u269 = pc0*u266;
            const double u270 = pc1*u267;
            const double u271 = pc2*u268;
            const double u272 = pa1*u267;
            const double u273 = pa2*u268;
            const double u274 = u272 + u273;
            const double u275 = pa0*u266 - u100 + u274;
            const double u276 = -u60;
            const double u277 = u266 + u276;
            const double u278 = pa0*u277;
            const double u279 = 8*pa1;
            const double u280 = 8*pa2;
            const double u281 = u25*u255;
            const double u282 = u267 + u28*u281;
            const double u283 = u268 + u281*u31;
            const double u284 = -u35;
            const double u285 = u255*u26 + u277 + u284;
            const double u286 = pa0*u285 + pa1*u282 + pa2*u283;
            const double u287 = pc0*u285 + pc1*u282 + pc2*u283;
            const double u288 = u51*u7;
            const double u289 = u178 - 4*xc1;
            const double u290 = u181 + 2*xb1;
            const double u291 = u184*u60;
            const double u292 = u27*u35;
            const double u293 = u172 + u175;
            const double u294 = u291 + u292;
            const double u295 = u51*u8;
            const double u296 = u288*u9;
            const double u297 = pc1*u295;
            const double u298 = pc0*u267;
            const double u299 = pc2*u296;
            const double u300 = pa0*u267;
            const double u301 = pa2*u296;
            const double u302 = u300 + u301;
            const double u303 = pa1*u295 - u186 + u302;
            const double u304 = u276 + u295;
            const double u305 = pa1*u304;
            const double u306 = 8*pa0;
            const double u307 = u255*u28;
            const double u308 = u296 + u307*u31;
            const double u309 = u255*u29 + u284 + u304;
            const double u310 = pa0*u282 + pa1*u309 + pa2*u308;
            const double u311 = pc0*u282 + pc1*u309 + pc2*u308;
            const double u312 = u220 - 4*xc2;
            const double u313 = u223 + 2*xb2;
            const double u314 = u226*u60;
            const double u315 = u30*u35;
            const double u316 = u214 + u217;
            const double u317 = u314 + u315;
            const double u318 = u10*u51;
            const double u319 = pc2*u318;
            const double u320 = pc0*u268;
            const double u321 = pc1*u296;
            const double u322 = pa0*u268;
            const double u323 = pa1*u296;
            const double u324 = u322 + u323;
            const double u325 = pa2*u318 - u228 + u324;
            const double u326 = u276 + u318;
            const double u327 = pa2*u326;
            const double u328 = u255*u32 + u284 + u326;
            const double u329 = pa0*u283 + pa1*u308 + pa2*u328;
            const double u330 = pc0*u283 + pc1*u308 + pc2*u328;
            const double u331 = (1.0/64.0)*mb*mc*u0*u48*u51;
            const double u332 = u23*u255*u331;
            const double u333 = u331*u35*u52;
            const double u334 = -u259;
            const double u335 = u158*u65;
            const double u336 = (1.0/64.0)*mb*mc*u0*u23*u35*u51;
            const double u337 = u14*u23;
            const double u338 = u17*u337;
            const double u339 = pa1*u338;
            const double u340 = pa2*u20;
            const double u341 = u337*u340;
            const double u342 = -u65;
            const double u343 = u15*u23 + u342;
            const double u344 = pa0*u343;
            const double u345 = u24*u307;
            const double u346 = u255*u31;
            const double u347 = u24*u346;
            const double u348 = u24*u281 + u35;
            const double u349 = pa0*u348 + pa1*u345 + pa2*u347;
            const double u350 = u20*u337;
            const double u351 = pc0*u348 + pc1*u345 + pc2*u347;
            const double u352 = -u290;
            const double u353 = u205*u65;
            const double u354 = pa0*u338;
            const double u355 = u17*u23;
            const double u356 = u340*u355;
            const double u357 = u18*u23 + u342;
            const double u358 = pa1*u357;
            const double u359 = u27*u281;
            const double u360 = u27*u346;
            const double u361 = u27*u307 + u35;
            const double u362 = pa0*u359 + pa1*u361 + pa2*u360;
            const double u363 = u20*u355;
            const double u364 = pc0*u359 + pc1*u361 + pc2*u360;
            const double u365 = -u313;
            const double u366 = u247*u65;
            const double u367 = pa0*u350;
            const double u368 = pa1*u363;
            const double u369 = u21*u23 + u342;
            const double u370 = pa2*u369;
            const double u371 = u281*u30;
            const double u372 = u30*u307;
            const double u373 = u30*u346 + u35;
            const double u374 = pa0*u371 + pa1*u372 + pa2*u373;
            const double u375 = pc0*u371 + pc1*u372 + pc2*u373;
            const double u376 = 3*u61;
            const double u377 = u75*u83;
            const double u378 = -u122*u61 + u85;
            const double u379 = 3*u78;
            const double u380 = -u122*u78 + u86;
            const double u381 = 3*u79;
            const double u382 = -u122*u79 + u87;
            const double u383 = u147*u97;
            dwdt[P0 + a * D + 0] += (1.0/64.0)*mb*mc*u0*u12*u23*u35*u4*u48 + (1.0/64.0)*mb*mc*u0*u23*u35*u51*(u22*u53 + u3*u62 - u33*u57 - u38*u61 + u41*u59 - u47*u68 + u55*u56 + u61*u64 - u61*u69 + u66*u67 - u66*u70) - 1.0/64.0*mb*mc*u0*u35*u48*u50*u51*u52 - u124*(u108*(u101 + u102*u103 + u104*u105 + u106*u107) + u111*(pb0*u103 + u109 + u110) + u113*u114 + u120*(u116 + 2*u117 + 2*u118 + 2*u119) + u123*(-u113*u122 - u88*(u115 + u117 + u121))) - u154*(-u125 - 2*u66) - u72*u95 + u96*u99 - u98*(u140*u156 + u145*u157 + u148*(pb0*u155 + u110) + u149*u157 - u153*u156 + u166*(u161 + u163 + u165) + u167*(u142*u164 + u143*u160 + u144*u162) + u168*(3*u161 + 3*u163 + 3*u165) + u169*(8*u161 + 8*u163 + 8*u165) - u170*(u150*u164 + u151*u160 + u152*u162));
            dwdt[P0 + a * D + 1] += (1.0/64.0)*mb*mc*u0*u12*u174*u23*u35*u48 - 1.0/64.0*mb*mc*u0*u176*u35*u48*u51*u52 + (1.0/64.0)*mb*mc*u0*u23*u35*u51*(u134*u67 - u134*u70 + u173*u62 + u177*u22 + u179*u56 - u180*u33 + u182*u41 - u183*u47 - u38*u78 + u64*u78 - u69*u78) - u124*(u108*(u102*u188 + u104*u189 + u106*u190 + u187) + u111*(pb1*u189 + u191 + u192) + u114*u194 + u120*(u196 + 2*u197 + 2*u198 + 2*u199) + u123*(-u122*u194 - u88*(u195 + u198 + u200))) + u135*u99 - u154*(-2*u134 - u201) - u185*u95 - u98*(u140*u203 + u145*u204 + u148*(pb1*u202 + u192) + u149*u204 - u153*u203 + u166*(u208 + u210 + u212) + u167*(u142*u207 + u143*u211 + u144*u209) + u168*(3*u208 + 3*u210 + 3*u212) + u169*(8*u208 + 8*u210 + 8*u212) - u170*(u150*u207 + u151*u211 + u152*u209));
            dwdt[P0 + a * D + 2] += (1.0/64.0)*mb*mc*u0*u12*u216*u23*u35*u48 - 1.0/64.0*mb*mc*u0*u218*u35*u48*u51*u52 + (1.0/64.0)*mb*mc*u0*u23*u35*u51*(u137*u67 - u137*u70 + u215*u62 + u219*u22 + u221*u56 - u222*u33 + u224*u41 - u225*u47 - u38*u79 + u64*u79 - u69*u79) - u124*(u108*(u102*u230 + u104*u231 + u106*u232 + u229) + u111*(pb2*u232 + u233 + u234) + u114*u236 + u120*(u238 + 2*u239 + 2*u240 + 2*u241) + u123*(-u122*u236 - u88*(u237 + u241 + u242))) + u138*u99 - u154*(-2*u137 - u243) - u227*u95 - u98*(u140*u245 + u145*u246 + u148*(pb2*u244 + u234) + u149*u246 - u153*u245 + u166*(u250 + u252 + u254) + u167*(u142*u249 + u143*u251 + u144*u253) + u168*(3*u250 + 3*u252 + 3*u254) + u169*(8*u250 + 8*u252 + 8*u254) - u170*(u150*u249 + u151*u251 + u152*u253));
            dwdt[P0 + b * D + 0] += (1.0/64.0)*mb*mc*u0*u12*u23*u3*u35*u48 + (1.0/64.0)*mb*mc*u0*u23*u25*u255*u48*u51 + (1.0/64.0)*mb*mc*u0*u23*u35*u51*(-u22*u53 + u257*u258 - u259*u42 + u259*u44 - u260*u38 + u260*u64 + u261*u262 - u263*u264 - u265*u69 + u33*u57 + u4*u62 + u47*u68 - u55*u56) - u124*(u108*(-u101 + u102*u266 + u104*u267 + u106*u268) + u111*(pb0*u266 + pb1*u267 + pb2*u268 - u109) + u114*u275 + u120*(-u116 + 2*u269 + 2*u270 + 2*u271) + u123*(-u122*u275 - u88*(-u115 + u269 + u270 + u271))) - u154*(-2*u260 - 2*u261) - u256*u95 + u265*u99 - u98*(u140*u286 + u145*u287 + u148*(pb0*u285 + pb1*u282 + pb2*u283) + u149*u287 - u153*u286 + u166*(u274 + u278) + u167*(u142*u277 + u143*u267 + u144*u268) + u168*(3*u272 + 3*u273 + 3*u278) + u169*(u267*u279 + u268*u280 + 8*u278) - u170*(u150*u277 + u151*u267 + u152*u268));
            dwdt[P0 + b * D + 1] += (1.0/64.0)*mb*mc*u0*u12*u173*u23*u35*u48 + (1.0/64.0)*mb*mc*u0*u23*u255*u28*u48*u51 + (1.0/64.0)*mb*mc*u0*u23*u35*u51*(u174*u62 - u177*u22 - u179*u56 + u180*u33 + u183*u47 + u258*u289 + u262*u292 - u264*u293 - u290*u42 + u290*u44 - u291*u38 + u291*u64 - u294*u69) - u124*(u108*(u102*u267 + u104*u295 + u106*u296 - u187) + u111*(pb0*u267 + pb1*u295 + pb2*u296 - u191) + u114*u303 + u120*(-u196 + 2*u297 + 2*u298 + 2*u299) + u123*(-u122*u303 - u88*(-u195 + u297 + u298 + u299))) - u154*(-2*u291 - 2*u292) - u288*u95 + u294*u99 - u98*(u140*u310 + u145*u311 + u148*(pb0*u282 + pb1*u309 + pb2*u308) + u149*u311 - u153*u310 + u166*(u302 + u305) + u167*(u142*u267 + u143*u304 + u144*u296) + u168*(3*u300 + 3*u301 + 3*u305) + u169*(u267*u306 + u280*u296 + 8*u305) - u170*(u150*u267 + u151*u304 + u152*u296));
            dwdt[P0 + b * D + 2] += (1.0/64.0)*mb*mc*u0*u12*u215*u23*u35*u48 + (1.0/64.0)*mb*mc*u0*u23*u255*u31*u48*u51 + (1.0/64.0)*mb*mc*u0*u23*u35*u51*(u216*u62 - u219*u22 - u221*u56 + u222*u33 + u225*u47 + u258*u312 + u262*u315 - u264*u316 - u313*u42 + u313*u44 - u314*u38 + u314*u64 - u317*u69) - u124*(u108*(u102*u268 + u104*u296 + u106*u318 - u229) + u111*(pb0*u268 + pb1*u296 + pb2*u318 - u233) + u114*u325 + u120*(-u238 + 2*u319 + 2*u320 + 2*u321) + u123*(-u122*u325 - u88*(-u237 + u319 + u320 + u321))) - u154*(-2*u314 - 2*u315) + u317*u99 - u51*u9*u95 - u98*(u140*u329 + u145*u330 + u148*(pb0*u283 + pb1*u308 + pb2*u328) + u149*u330 - u153*u329 + u166*(u324 + u327) + u167*(u142*u268 + u143*u296 + u144*u326) + u168*(3*u322 + 3*u323 + 3*u327) + u169*(u268*u306 + u279*u296 + 8*u327) - u170*(u150*u268 + u151*u296 + u152*u326));
            dwdt[P0 + c * D + 0] += -u154*(-2*u126 - 2*u335) + u24*u332 + u333*u50 + u336*(u126*u262 - u126*u69 - u257*u258 + u263*u264 - u334*u42 + u334*u44 + u335*u67 - u335*u70 - u41*u59) - u98*(u140*u349 + u145*u351 + u148*(pb0*u348 + pb1*u345 + pb2*u347) + u149*u351 - u153*u349 + u166*(u339 + u341 + u344) + u167*(u142*u343 + u143*u338 + u144*u350) + u168*(3*u339 + 3*u341 + 3*u344) + u169*(u279*u338 + u280*u350 + 8*u344) - u170*(u150*u343 + u151*u338 + u152*u350)) - u99*(-u126 - u335);
            dwdt[P0 + c * D + 1] += -u154*(-2*u128 - 2*u353) + u176*u333 + u27*u332 + u336*(u128*u262 - u128*u69 - u182*u41 - u258*u289 + u264*u293 - u352*u42 + u352*u44 + u353*u67 - u353*u70) - u98*(u140*u362 + u145*u364 + u148*(pb0*u359 + pb1*u361 + pb2*u360) + u149*u364 - u153*u362 + u166*(u354 + u356 + u358) + u167*(u142*u338 + u143*u357 + u144*u363) + u168*(3*u354 + 3*u356 + 3*u358) + u169*(u280*u363 + u306*u338 + 8*u358) - u170*(u150*u338 + u151*u357 + u152*u363)) - u99*(-u128 - u353);
            dwdt[P0 + c * D + 2] += -u154*(-2*u130 - 2*u366) + u218*u333 + u30*u332 + u336*(u130*u262 - u130*u69 - u224*u41 - u258*u312 + u264*u316 - u365*u42 + u365*u44 + u366*u67 - u366*u70) - u98*(u140*u374 + u145*u375 + u148*(pb0*u371 + pb1*u372 + pb2*u373) + u149*u375 - u153*u374 + u166*(u367 + u368 + u370) + u167*(u142*u350 + u143*u363 + u144*u369) + u168*(3*u367 + 3*u368 + 3*u370) + u169*(u279*u363 + u306*u350 + 8*u370) - u170*(u150*u350 + u151*u363 + u152*u369)) - u99*(-u130 - u366);
            dwdt[a * D + 0] += u124*(-pb0*u76 + u123*u378 + u376*u377 - u77*(u102 - u125*u80)) + u98*(u127*u140 - u127*u153 + u166*u96 + u168*(u376 + 3*u66) + u169*(8*u61 + 8*u66));
            dwdt[a * D + 1] += u124*(-pb1*u76 + u123*u380 + u377*u379 - u77*(u104 - u201*u80)) + u98*(u129*u140 - u129*u153 + u135*u166 + u168*(3*u134 + u379) + u169*(8*u134 + 8*u78));
            dwdt[a * D + 2] += u124*(-pb2*u76 + u123*u382 + u377*u381 - u77*(u106 - u243*u80)) + u98*(u131*u140 - u131*u153 + u138*u166 + u168*(3*u137 + u381) + u169*(8*u137 + 8*u79));
            dwdt[b * D + 0] += u124*(-pa0*u76 + 3*u5*u60*u73*u74*u80) + u127*u383;
            dwdt[b * D + 1] += u124*(-pa1*u76 + 3*u60*u7*u73*u74*u80) + u129*u383;
            dwdt[b * D + 2] += u124*(-pa2*u76 + 3*u60*u73*u74*u80*u9) + u131*u383;
            dwdt[c * D + 0] += u124*(-u378*u81 + u73*u84*(u306 - u61*u88)) + u98*(u127*u145 + u127*u149 + u167*(4*u61 + 4*u66) - u170*(16*u61 + 16*u66));
            dwdt[c * D + 1] += u124*(-u380*u81 + u73*u84*(u279 - u78*u88)) + u98*(u129*u145 + u129*u149 + u167*(4*u134 + 4*u78) - u170*(16*u134 + 16*u78));
            dwdt[c * D + 2] += u124*(-u382*u81 + u73*u84*(u280 - u79*u88)) + u98*(u131*u145 + u131*u149 + u167*(4*u137 + 4*u79) - u170*(16*u137 + 16*u79));

                }
            }
        }
    }
}


static inline void add_2pn_fourbody_inv_r_force(const PairCache *cache, double *dwdt, int u, int v,
    double prefactor_without_this_inv_r)
{
    const int D = cache->num_dim;
    const int P0 = cache->array_half;

    const double inv_r_uv = pair_cache_inv_r(cache, u, v);
    const double coeff = prefactor_without_this_inv_r * inv_r_uv * inv_r_uv;

    for (int i = 0; i < D; i++) {
        const double force_i = coeff * pair_cache_n(cache, u, v, i);

        dwdt[P0 + u * D + i] += force_i;
        dwdt[P0 + v * D + i] -= force_i;
    }
}


/**
 * @brief Add the analytic non-UTT4 four-body 2PN Hamiltonian contribution.
 *
 * This implements the position-only four-body terms from H2PN_base_complex:
 *
 *   H_4body =
 *     sum_{a,b,c,d} -3/8  m_a m_b m_c m_d / (r_ab r_bc r_cd)
 *       with b != a, c != b, d != c
 *
 *   + sum_{a,b,c,d} -1/4  m_a m_b m_c m_d / (r_ab r_ac r_ad)
 *       with b != a, c != a, d != a.
 *
 * These terms are independent of momenta, so they contribute only to dp/dt:
 *
 *   dx/dt += 0,
 *   dp/dt += -dH/dx.
 */
static void add_eom_2pn_fourbody_analytic(const PairCache *cache, double *dwdt)
{
    const ActiveList *active = &cache->active;

    for (int ia = 0; ia < active->num_active; ia++) {
        const int a = active->ids[ia];
        const double ma = cache->m[a];

        for (int ib = 0; ib < active->num_active; ib++) {
            const int b = active->ids[ib];
            if (b == a)
                continue;

            const double mb = cache->m[b];
            const double inv_r_ab = pair_cache_inv_r(cache, a, b);

            for (int ic = 0; ic < active->num_active; ic++) {
                const int c = active->ids[ic];
                const double mc = cache->m[c];

                /*
                 * Chain term:
                 *
                 *   -3/8 m_a m_b m_c m_d / (r_ab r_bc r_cd)
                 *
                 * Conditions exactly match H2PN_base_complex:
                 *   b != a, c != b, d != c.
                 */
                if (c != b) {
                    const double inv_r_bc = pair_cache_inv_r(cache, b, c);

                    for (int id = 0; id < active->num_active; id++) {
                        const int d = active->ids[id];
                        if (d == c)
                            continue;

                        const double md = cache->m[d];
                        const double inv_r_cd = pair_cache_inv_r(cache, c, d);

                        const double K = -0.375 * ma * mb * mc * md;

                        add_2pn_fourbody_inv_r_force(
                            cache, dwdt, a, b,
                            K * inv_r_bc * inv_r_cd);

                        add_2pn_fourbody_inv_r_force(
                            cache, dwdt, b, c,
                            K * inv_r_ab * inv_r_cd);

                        add_2pn_fourbody_inv_r_force(
                            cache, dwdt, c, d,
                            K * inv_r_ab * inv_r_bc);
                    }
                }

                /*
                 * Star term:
                 *
                 *   -1/4 m_a m_b m_c m_d / (r_ab r_ac r_ad)
                 *
                 * Conditions exactly match H2PN_base_complex:
                 *   b != a, c != a, d != a.
                 */
                if (c != a) {
                    const double inv_r_ac = pair_cache_inv_r(cache, a, c);

                    for (int id = 0; id < active->num_active; id++) {
                        const int d = active->ids[id];
                        if (d == a)
                            continue;

                        const double md = cache->m[d];
                        const double inv_r_ad = pair_cache_inv_r(cache, a, d);

                        const double K = -0.25 * ma * mb * mc * md;

                        add_2pn_fourbody_inv_r_force(
                            cache, dwdt, a, b,
                            K * inv_r_ac * inv_r_ad);

                        add_2pn_fourbody_inv_r_force(
                            cache, dwdt, a, c,
                            K * inv_r_ab * inv_r_ad);

                        add_2pn_fourbody_inv_r_force(
                            cache, dwdt, a, d,
                            K * inv_r_ab * inv_r_ac);
                    }
                }
            }
        }
    }
}


/**
 * @brief Right-hand side of the N-body equations of motion up to 2.5PN order
 * 
 * @param[in]       t           Time (currently not used, but kept for completeness)
 * @param[in]       w           State of the system, w = [positions, momenta]
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[out]      dwdt        Right-hand side of the equations of motion
 */
void rhs_pn_nbody(double t, double* w, struct ode_params* ode_params, double* dwdt)
{
    (void)t;    // Unused

    // --------------------------------------------------------------------------------------------
    // Initialize the arrays
    // --------------------------------------------------------------------------------------------
    int num_bodies = ode_params->num_bodies_initial;
    int num_dim = ode_params->num_dim;
    int array_half = num_bodies * num_dim; 
    double result, temp;

    // Cached active-body, momentum and pair-geometry data.
    //
    // PairCache owns an ActiveList in cache.active.  All expensive RHS loops below now iterate over
    // cache.active.ids instead of scanning num_bodies_initial and skipping inactive bodies.
    PairCache cache;
    pair_cache_build(&cache, w, ode_params);
    const ActiveList *active = &cache.active;

    // Masses
    double m[num_bodies];
    double inv_m[num_bodies];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        m[a] = cache.m[a];
        inv_m[a] = cache.inv_m[a];
    }

    // Momenta
    double p[num_bodies][num_dim];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        for (int i = 0; i < num_dim; i++)
            p[a][i] = pair_cache_p(&cache, a, i);
    }

    // Relative positions and distances
    double x_rel[num_bodies][num_bodies][num_dim]; 
    double n[num_bodies][num_bodies][num_dim];
    double r[num_bodies][num_bodies];
    double inv_r[num_bodies][num_bodies];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        for (int ib = ia; ib < active->num_active; ib++) {
            int b = active->ids[ib];

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

    // Time derivatives
    double p_dot[num_bodies][num_dim];
    double x_dot[num_bodies][num_dim];
    for (int ia = 0; ia < active->num_active; ia++) {
        int a = active->ids[ia];
        for (int i = 0; i < num_dim; i++) {
            p_dot[a][i] = 0.0;
            x_dot[a][i] = 0.0;
        }
    }

    // Set ODE right-hand side initially to zero.
    // This also guarantees that inactive bodies have zero RHS, because all contribution loops below
    // iterate only over active->ids.
    for (int i = 0; i < 2 * array_half; i++)   
        dwdt[i] = 0.0;

    // --------------------------------------------------------------------------------------------
    // Add 0PN (Newtonian) terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[0]) {
        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];

            // Velocities
            for (int i = 0; i < num_dim; i++) {
                result = w[array_half + a * num_dim + i] * inv_m[a];
                dwdt[a * num_dim + i] += result;
                x_dot[a][i] += result;
            }

            // Accelerations.  Use ib = ia + 1 to keep the original pair-once symmetric update.
            for (int ib = ia + 1; ib < active->num_active; ib++) {
                int b = active->ids[ib];

                for (int i = 0; i < num_dim; i++) {
                    result = -m[a] * m[b] * inv_r[a][b] * inv_r[a][b] * n[a][b][i];
                    dwdt[array_half + a * num_dim + i] += result;
                    dwdt[array_half + b * num_dim + i] += -result;
                    p_dot[a][i] += result;
                    p_dot[b][i] += -result;
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------------
    // Add 1PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[1]) {
        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];

            double pa_dot_pa = cache.p2[a];

            for (int i = 0; i < num_dim; i++) {

                // Velocities
                result = -0.5 * pa_dot_pa * inv_m[a] * inv_m[a] * inv_m[a] * p[a][i];
                x_dot[a][i] += result;
                dwdt[a * num_dim + i] += result;

                for (int ib = 0; ib < active->num_active; ib++) {
                    int b = active->ids[ib];

                    double pa_dot_pb = pair_cache_p_dot(&cache, a, b);
                    double pb_dot_pb = cache.p2[b];
                    double nab_dot_pa = pair_cache_n_dot_p(&cache, a, b, a);
                    double nab_dot_pb = pair_cache_n_dot_p(&cache, a, b, b);

                    if (b != a) {
                        temp = -0.5 * inv_r[a][b];
                        result = temp * (6 * m[b] * inv_m[a] * p[a][i] 
                            - 7 * p[b][i] - nab_dot_pb * n[a][b][i]);
                        x_dot[a][i] += result;
                        dwdt[a * num_dim + i] += result;

                        // Accelerations 
                        temp *= inv_r[a][b];
                        result = temp * (3 * m[b] * inv_m[a] * pa_dot_pa 
                            + 3 * m[a] * inv_m[b] * pb_dot_pb 
                            - 7 * pa_dot_pb - 3 * nab_dot_pa * nab_dot_pb) * n[a][b][i];
                        p_dot[a][i] += result;                       
                        dwdt[array_half + a * num_dim + i] += result;

                        result = temp * (nab_dot_pb * p[a][i] + nab_dot_pa * p[b][i]);
                        p_dot[a][i] += result;
                        dwdt[array_half + a * num_dim + i] += result;
                    }

                    for (int ic = 0; ic < active->num_active; ic++) {
                        int c = active->ids[ic];

                        temp = m[a] * m[b] * m[c] * inv_r[a][b] * inv_r[a][b] * n[a][b][i];

                        if (b != a && c != a) {
                            result = temp * inv_r[a][c];
                            p_dot[a][i] += result;
                            dwdt[array_half + a * num_dim + i] += result;
                        }

                        if (b != a && c != b) {
                            result = temp * inv_r[b][c];
                            p_dot[a][i] += result;
                            dwdt[array_half + a * num_dim + i] += result;
                        }
                    }
                }
            }
        }
    }

    // --------------------------------------------------------------------------------------------
    // Add 2PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[2]) {
        add_eom_2pn_onebody_analytic(&cache, dwdt);
        add_eom_2pn_pair_analytic(&cache, dwdt);
        add_eom_2pn_triple_analytic(&cache, w, dwdt);
        add_eom_2pn_fourbody_analytic(&cache, dwdt);

        // Add the remaining contributions from H2PN without UTT4.
        // update_eom_hamiltonian_cs(w, H2PN_base_complex, 1e-30, ode_params, dwdt);

        // If not using impulse splitting, add UTT4 contributions directly to dp/dt
        #if HAVE_CUBA
        if (ode_params->include_utt4 && !ode_params->use_impulse_method)
        {
            double dUdx[array_half];
            compute_dUTT4_dx(w, ode_params, dUdx);
            for (int ia = 0; ia < active->num_active; ia++) {
                int a = active->ids[ia];
                for (int i = 0; i < num_dim; i++) {
                    int idx = a * num_dim + i;
                    dwdt[array_half + idx] -= dUdx[idx];
                }
            }
        }
        #endif
    }

    // --------------------------------------------------------------------------------------------
    // Add 2.5PN terms
    // --------------------------------------------------------------------------------------------
    if (ode_params->pn_terms[3]) {

        // Initialize arrays
        double x_rel_dot[num_bodies][num_bodies][num_dim];

        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];
            for (int ib = 0; ib < active->num_active; ib++) {
                int b = active->ids[ib];

                for (int i = 0; i < num_dim; i++)
                    x_rel_dot[a][b][i] = x_dot[a][i] - x_dot[b][i];
            }
        }

        double chi_dot[num_dim][num_dim];
        double dp_chi[num_bodies][num_dim][num_dim][num_dim];
        double dx_chi[num_bodies][num_dim][num_dim][num_dim];

        for (int i = 0; i < num_dim; i++)
            for (int j = 0; j < num_dim; j++)
                chi_dot[i][j] = 0.0;

        for (int ia = 0; ia < active->num_active; ia++) { 
            int a = active->ids[ia];
            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[a][i][j][k] = 0.0;
                        dx_chi[a][i][j][k] = 0.0;
                    }
                }
            }
        }

        // Compute Chi values
        for (int i = 0; i < num_dim; i++) {
            for (int j = 0; j < num_dim; j++) {
                for (int ia = 0; ia < active->num_active; ia++) {
                    int a = active->ids[ia];

                    chi_dot[i][j] += 2.0 * inv_m[a] * (
                        2 * dot_product(p_dot[a], p[a], num_dim) * delta(i, j) 
                        - 3 * (p_dot[a][i] * p[a][j] + p[a][i] * p_dot[a][j]) );

                    for (int ib = 0; ib < active->num_active; ib++) {
                        int b = active->ids[ib];

                        if (b != a) {
                            chi_dot[i][j] += m[a] * m[b] * inv_r[a][b] * inv_r[a][b] * (3 * 
                                (x_rel_dot[a][b][i] * n[a][b][j] + n[a][b][i] * x_rel_dot[a][b][j])
                                + dot_product(n[a][b], x_rel_dot[a][b], num_dim) * (delta(i, j) 
                                - 9 * n[a][b][i] * n[a][b][j]) );
                        }
                    }
                }
            }
        }

        for (int ic = 0; ic < active->num_active; ic++) {
            int c = active->ids[ic];

            for (int i = 0; i < num_dim; i++) {
                for (int j = 0; j < num_dim; j++) {
                    for (int k = 0; k < num_dim; k++) {
                        dp_chi[c][i][j][k] += 2.0 / m[c] * (2 * p[c][k] * delta(i, j) 
                            - 3 * (p[c][j] * delta(i, k) + p[c][i] * delta(j, k)));

                        for (int ia = 0; ia < active->num_active; ia++) {
                            int a = active->ids[ia];

                            for (int ib = 0; ib < active->num_active; ib++) {
                                int b = active->ids[ib];

                                if (b != a) {
                                    dx_chi[c][i][j][k] += m[a] * m[b] * inv_r[a][b] * inv_r[a][b] *
                                        (delta(a, c) - delta(b, c)) * 
                                        (3 * (delta(i, k) * n[a][b][j] 
                                        + delta(j, k) * n[a][b][i]) 
                                        - 9 * n[a][b][k] * n[a][b][i] * n[a][b][j] 
                                        + delta(i, j) * n[a][b][k]);
                                }
                            }
                        }
                    }    
                }
            }
        }

        // Add contribution to the ODE right-hand side
        for (int ia = 0; ia < active->num_active; ia++) {
            int a = active->ids[ia];

            for (int k = 0; k < num_dim; k++) {
                for (int i = 0; i < num_dim; i++) {
                    for (int j = 0; j < num_dim; j++) {
                        dwdt[a * num_dim + k] += 1/45.0 * chi_dot[i][j] * dp_chi[a][i][j][k];
                        dwdt[array_half + a * num_dim + k] += -1/45.0 * chi_dot[i][j] 
                            * dx_chi[a][i][j][k];
                    }
                }
            }
        }
    }

    pair_cache_free(&cache);
}


/**
 * @brief Adds contribution from a Hamiltonian to the right-hand side of the equations of motion.
 * 
 * Adds contribution from a Hamiltonian to the right-hand side of the equations of motion,
 * according to dx/dt = dH/dp, dp/dt = -dH/dx. The derivatives of the Hamiltonian are computed
 * numerically using a complex-step derivative. The Hamiltonian must be of type complex double
 * with arguments (w, ode_params, p_flag), where p_flag just ignores all the terms that do not
 * have a momentum dependence for the computation of dH/dp.
 * 
 * @param[in]       w           State of the system, w = [positions, momenta]
 * @param[in]       H           Complex-valued Hamiltonian
 * @param[in]       h           Complex step size
 * @param[in]       ode_params  Parameter struct containing general information about the system
 * @param[out]      dwdt        Right-hand side of the ODE
 */
void update_eom_hamiltonian_cs(double *w, c_hamiltonian H, double h, struct ode_params* ode_params,
    double *dwdt)
{
    int num_dim = ode_params->num_dim;
    int num_bodies = ode_params->num_bodies_initial;
    int array_half = num_dim * num_bodies;
    int total_dim = 2 * array_half;
    complex double w_c[total_dim];
    complex double H_cs_val;
    double dHdw[total_dim];

    ActiveList active;
    active_list_init(&active, ode_params);

    // Copy original array to w_c and initialize dHdw
    for (int i = 0; i < total_dim; ++i) {
        w_c[i] = (complex double)w[i];
        dHdw[i] = 0.0;
    }

    // Position derivatives: dH/dx. These are needed for dp/dt = -dH/dx.
    for (int ia = 0; ia < active.num_active; ia++) {
        int a = active.ids[ia];

        for (int k = 0; k < num_dim; k++) {
            int idx = a * num_dim + k;

            // Add tiny imaginary step in coordinate idx
            w_c[idx] += I * h;

            // p_flag = 0: keep position-only Hamiltonian terms
            H_cs_val = H(w_c, ode_params, 0);
            dHdw[idx] = cimag(H_cs_val) / h;

            // Restore original value
            w_c[idx] = (complex double)w[idx];
        }
    }

    // Momentum derivatives: dH/dp. These are needed for dx/dt = dH/dp.
    for (int ia = 0; ia < active.num_active; ia++) {
        int a = active.ids[ia];

        for (int k = 0; k < num_dim; k++) {
            int idx = array_half + a * num_dim + k;

            // Add tiny imaginary step in momentum component idx
            w_c[idx] += I * h;

            // p_flag = 1: skip Hamiltonian terms without momentum dependence
            H_cs_val = H(w_c, ode_params, 1);
            dHdw[idx] = cimag(H_cs_val) / h;

            // Restore original value
            w_c[idx] = (complex double)w[idx];
        }
    }

    // Compute dwdt for active bodies only
    for (int ia = 0; ia < active.num_active; ia++) {
        int a = active.ids[ia];

        for (int k = 0; k < num_dim; k++) {
            int x_idx = a * num_dim + k;
            int p_idx = array_half + a * num_dim + k;

            dwdt[x_idx] += dHdw[p_idx];
            dwdt[p_idx] += -dHdw[x_idx];
        }
    }

    active_list_free(&active);
}
