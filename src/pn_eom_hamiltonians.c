#include <math.h>
#include "pn_eom.h"
#include "pn_eom_hamiltonians.h"
#include "utils.h"


double H0PN(double* w, struct ode_params* params) {
     int array_half = params->dim * params->num_bodies;
     double rel_dist2, v2;
     double H = 0;

     for (int a = 0; a < params->num_bodies; a++) {
          // Kinetic energy
          v2 = 0.0;
          for (int i = 0; i < params->dim; i++)
               v2 += pow(w[array_half + params->dim * a + i], 2);
          H += 0.5 * params->masses[a] * v2;

          // Potential energy
          for (int b = a+1; b < params->num_bodies; b++) {
               rel_dist2 = 0.0;
               for (int i = 0; i < params->dim; i++)
                    rel_dist2 += pow(w[params->dim * a + i] - w[params->dim * b + i], 2);
               H -= params->masses[a] * params->masses[b] / sqrt(rel_dist2);
          }
     }
     return H;
}


double H1PN(double* w, struct ode_params* params) {
    int a, b, c, i;
    double m_a, m_b, m_c;
    double pa_dot_pa, pb_dot_pb, pa_dot_pb;
    double dx, r_ab, r_ac, ni, na_dot_pa, na_dot_pb;
    double p_ai, p_bi, dx_ac;
    int array_half = params->dim * params->num_bodies;
    double H = 0.0;

    // Compute kinetic energy and potential energy
    for (a = 0; a < params->num_bodies; a++) {
        m_a = params->masses[a];
        pa_dot_pa = 0.0;

        for (i = 0; i < params->dim; i++) {
            p_ai = m_a * w[array_half + a * params->dim + i];
            pa_dot_pa += p_ai * p_ai;
        }

        H += -0.125 * m_a * pow(pa_dot_pa / (m_a * m_a), 2);

        for (b = 0; b < params->num_bodies; b++) {
            if (b == a) continue;

            m_b = params->masses[b];
            r_ab = 0.0;
            pb_dot_pb = 0.0;
            pa_dot_pb = 0.0;
            na_dot_pa = 0.0;
            na_dot_pb = 0.0;

            for (i = 0; i < params->dim; i++) {
                dx = w[a * params->dim + i] - w[b * params->dim + i];
                r_ab += dx * dx;
            }

            r_ab = sqrt(r_ab);

            for (i = 0; i < params->dim; i++) {
                p_ai = m_a * w[array_half + a * params->dim + i];
                p_bi = m_b * w[array_half + b * params->dim + i];

                dx = w[a * params->dim + i] - w[b * params->dim + i];
                ni = dx / r_ab;

                pb_dot_pb += p_bi * p_bi;
                pa_dot_pb += p_ai * p_bi;
                na_dot_pa += ni * p_ai;
                na_dot_pb += ni * p_bi;
            }

            H += -0.25 * m_a * m_b / r_ab * (6 * pa_dot_pa / (m_a * m_a) - 7 * pa_dot_pb / (m_a * m_b) - 
                  (na_dot_pa * na_dot_pb) / (m_a * m_b));

            if (params->num_bodies > 2) {
                for (c = 0; c < params->num_bodies; c++) {
                    if (c == a) continue;

                    m_c = params->masses[c];

                    r_ac = 0.0;

                    for (i = 0; i < params->dim; i++) {
                        dx_ac = w[a * params->dim + i] - w[c * params->dim + i];
                        r_ac += dx_ac * dx_ac;
                    }
                    r_ac = sqrt(r_ac);

                    H += 0.5 * m_a * m_b * m_c / (r_ab * r_ac);
                }
            }
        }
    }
    return H;
}




double H2PN_threebody(double* w, struct ode_params* params) {
     double m1 = params->masses[0];
     double m2 = params->masses[1];
     double m3 = params->masses[2];

     double x12 = w[0]-w[3];
     double x13 = w[0]-w[6];
     double x23 = w[3]-w[6];
     double y12 = w[1]-w[4];
     double y13 = w[1]-w[7];
     double y23 = w[4]-w[7];
     double z12 = w[2]-w[5];
     double z13 = w[2]-w[8];
     double z23 = w[5]-w[8];

     double p1x = w[9] * m1;
     double p1y = w[10] * m1;
     double p1z = w[11] * m1;
     double p2x = w[12] * m2;
     double p2y = w[13] * m2;
     double p2z = w[14] * m2;
     double p3x = w[15] * m3;
     double p3y = w[16] * m3;
     double p3z = w[17] * m3;

     double r12 = sqrt( pow(x12,2) + pow(y12,2) + pow(z12,2) );
     double r13 = sqrt( pow(x13,2) + pow(y13,2) + pow(z13,2) );
     double r23 = sqrt( pow(x23,2) + pow(y23,2) + pow(z23,2) );

     double n12x = x12/r12;
     double n12y = y12/r12;
     double n12z = z12/r12;
     double n13x = x13/r13;
     double n13y = y13/r13;
     double n13z = z13/r13;
     double n23x = x23/r23;
     double n23y = y23/r23;
     double n23z = z23/r23;

     return (pow(pow(p1x,2) + pow(p1y,2) + pow(p1z,2),3)/pow(m1,5) + pow(pow(p2x,2) + pow(p2y,2) + pow(p2z,2),3)/pow(m2,5) + 
          pow(pow(p3x,2) + pow(p3y,2) + pow(p3z,2),3)/pow(m3,5))/16. + 
     ((-2*pow(m1,2)*pow(m2,2))/pow(r12,3) - (2*pow(m1,2)*pow(m3,2))/pow(r13,3) - (2*pow(m2,2)*pow(m3,2))/pow(r23,3))/4. + 
     (3*((2*pow(m1,2)*pow(m2,2))/pow(r12,3) + (2*pow(m1,2)*pow(m3,2))/pow(r13,3) + (2*pow(m2,2)*pow(m3,2))/pow(r23,3)))/8. + 
     ((pow(m1,3)*m2)/pow(r12,3) + (m1*pow(m2,3))/pow(r12,3) + (pow(m1,3)*m3)/pow(r13,3) + (m1*pow(m3,3))/pow(r13,3) + 
          (pow(m2,3)*m3)/pow(r23,3) + (m2*pow(m3,3))/pow(r23,3))/2. + 
     ((pow(m1,2)*m2*((pow(p1x,2) + pow(p1y,2) + pow(p1z,2))/pow(m1,2) - (2*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + 
               (pow(p2x,2) + pow(p2y,2) + pow(p2z,2))/pow(m2,2)))/pow(r12,2) + 
          (m1*pow(m2,2)*((pow(p1x,2) + pow(p1y,2) + pow(p1z,2))/pow(m1,2) - (2*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + 
               (pow(p2x,2) + pow(p2y,2) + pow(p2z,2))/pow(m2,2)))/pow(r12,2) + 
          (pow(m1,2)*m3*((pow(p1x,2) + pow(p1y,2) + pow(p1z,2))/pow(m1,2) - (2*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) + 
               (pow(p3x,2) + pow(p3y,2) + pow(p3z,2))/pow(m3,2)))/pow(r13,2) + 
          (m1*pow(m3,2)*((pow(p1x,2) + pow(p1y,2) + pow(p1z,2))/pow(m1,2) - (2*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) + 
               (pow(p3x,2) + pow(p3y,2) + pow(p3z,2))/pow(m3,2)))/pow(r13,2) + 
          (pow(m2,2)*m3*((pow(p2x,2) + pow(p2y,2) + pow(p2z,2))/pow(m2,2) - (2*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + 
               (pow(p3x,2) + pow(p3y,2) + pow(p3z,2))/pow(m3,2)))/pow(r23,2) + 
          (m2*pow(m3,2)*((pow(p2x,2) + pow(p2y,2) + pow(p2z,2))/pow(m2,2) - (2*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + 
               (pow(p3x,2) + pow(p3y,2) + pow(p3z,2))/pow(m3,2)))/pow(r23,2))/4. + 
     ((pow(m1,2)*m2*((-14*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/pow(m1,2) + (2*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/(m1*m2) - 
               ((pow(n12x,2) + pow(n12y,2) + pow(n12z,2))*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/pow(m1,2) + 
               (5*(pow(n12x,2) + pow(n12y,2) + pow(n12z,2))*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) + 
               (2*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(-(n12x*p2x) - n12y*p2y - n12z*p2z))/(m1*m2)))/pow(r12,2) + 
          (m1*pow(m2,2)*((2*(n12x*p1x + n12y*p1y + n12z*p1z)*(n12x*p2x + n12y*p2y + n12z*p2z))/(m1*m2) - 
               (14*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/pow(m2,2) + (2*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/(m1*m2) - 
               ((pow(n12x,2) + pow(n12y,2) + pow(n12z,2))*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/pow(m2,2) + 
               (5*(pow(n12x,2) + pow(n12y,2) + pow(n12z,2))*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2)))/pow(r12,2) + 
          (m1*m2*m3*((2*(n12x*p1x + n12y*p1y + n12z*p1z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m1*m3) + 
               (2*(n12x*p2x + n12y*p2y + n12z*p2z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m1*m3) - 
               (14*(n12x*p3x + n12y*p3y + n12z*p3z)*(n13x*p3x + n13y*p3y + n13z*p3z))/pow(m3,2) - 
               ((n12x*n13x + n12y*n13y + n12z*n13z)*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/pow(m3,2) + 
               (5*(n12x*n13x + n12y*n13y + n12z*n13z)*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r12,2) + 
          (m1*m2*m3*((2*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m2*m3) + 
               (2*(-(n12x*p2x) - n12y*p2y - n12z*p2z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m2*m3) - 
               (14*(-(n12x*p3x) - n12y*p3y - n12z*p3z)*(n23x*p3x + n23y*p3y + n23z*p3z))/pow(m3,2) - 
               ((-(n12x*n23x) - n12y*n23y - n12z*n23z)*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/pow(m3,2) + 
               (5*(-(n12x*n23x) - n12y*n23y - n12z*n23z)*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r12,2) + 
          (pow(m1,2)*m3*((-14*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/pow(m1,2) + (2*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/(m1*m3) - 
               ((pow(n13x,2) + pow(n13y,2) + pow(n13z,2))*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/pow(m1,2) + 
               (5*(pow(n13x,2) + pow(n13y,2) + pow(n13z,2))*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) + 
               (2*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z))/(m1*m3)))/pow(r13,2) + 
          (m1*m2*m3*((2*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n23x*p2x) - n23y*p2y - n23z*p2z))/(m2*m3) - 
               (14*(-(n13x*p2x) - n13y*p2y - n13z*p2z)*(-(n23x*p2x) - n23y*p2y - n23z*p2z))/pow(m2,2) - 
               ((n13x*n23x + n13y*n23y + n13z*n23z)*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/pow(m2,2) + 
               (5*(n13x*n23x + n13y*n23y + n13z*n23z)*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) + 
               (2*(-(n23x*p2x) - n23y*p2y - n23z*p2z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z))/(m2*m3)))/pow(r13,2) + 
          (m1*m2*m3*((2*(n13x*p1x + n13y*p1y + n13z*p1z)*(n12x*p2x + n12y*p2y + n12z*p2z))/(m1*m2) - 
               ((n12x*n13x + n12y*n13y + n12z*n13z)*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/pow(m2,2) - 
               (14*(n12x*p2x + n12y*p2y + n12z*p2z)*(n13x*p2x + n13y*p2y + n13z*p2z))/pow(m2,2) + 
               (5*(n12x*n13x + n12y*n13y + n12z*n13z)*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) + 
               (2*(n12x*p2x + n12y*p2y + n12z*p2z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m1*m2)))/pow(r13,2) + 
          (m1*pow(m3,2)*((2*(n13x*p1x + n13y*p1y + n13z*p1z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m1*m3) - 
               (14*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/pow(m3,2) + (2*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/(m1*m3) - 
               ((pow(n13x,2) + pow(n13y,2) + pow(n13z,2))*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/pow(m3,2) + 
               (5*(pow(n13x,2) + pow(n13y,2) + pow(n13z,2))*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r13,2) + 
          (m1*m2*m3*(-(((n13x*n23x + n13y*n23y + n13z*n23z)*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/pow(m1,2)) - 
               (14*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n23x*p1x) - n23y*p1y - n23z*p1z))/pow(m1,2) + 
               (5*(n13x*n23x + n13y*n23y + n13z*n23z)*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) + 
               (2*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n23x*p2x) - n23y*p2y - n23z*p2z))/(m1*m3) + 
               (2*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z))/(m1*m3)))/pow(r23,2) + 
          (pow(m2,2)*m3*((-14*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/pow(m2,2) + (2*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/(m2*m3) - 
               ((pow(n23x,2) + pow(n23y,2) + pow(n23z,2))*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/pow(m2,2) + 
               (5*(pow(n23x,2) + pow(n23y,2) + pow(n23z,2))*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) + 
               (2*(-(n23x*p2x) - n23y*p2y - n23z*p2z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z))/(m2*m3)))/pow(r23,2) + 
          (m1*m2*m3*(-(((-(n12x*n23x) - n12y*n23y - n12z*n23z)*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/pow(m1,2)) - 
               (14*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(n23x*p1x + n23y*p1y + n23z*p1z))/pow(m1,2) + 
               (5*(-(n12x*n23x) - n12y*n23y - n12z*n23z)*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) + 
               (2*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(n23x*p2x + n23y*p2y + n23z*p2z))/(m1*m2) + 
               (2*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m1*m2)))/pow(r23,2) + 
          (m2*pow(m3,2)*((2*(n23x*p2x + n23y*p2y + n23z*p2z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m2*m3) - 
               (14*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/pow(m3,2) + (2*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/(m2*m3) - 
               ((pow(n23x,2) + pow(n23y,2) + pow(n23z,2))*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/pow(m3,2) + 
               (5*(pow(n23x,2) + pow(n23y,2) + pow(n23z,2))*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r23,2))/8. + 
     ((m1*m2*((10*pow(pow(p1x,2) + pow(p1y,2) + pow(p1z,2),2))/pow(m1,4) - 
               (2*pow(pow(p1x,2) + pow(p1y,2) + pow(p1z,2),2))/(pow(m1,2)*pow(m2,2)) - 
               (3*pow(n12x*p1x + n12y*p1y + n12z*p1z,2)*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/(pow(m1,2)*pow(m2,2)) + 
               (10*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2))*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/(pow(m1,2)*pow(m2,2)) - 
               (12*(n12x*p1x + n12y*p1y + n12z*p1z)*(n12x*p2x + n12y*p2y + n12z*p2z)*(p1x*p2x + p1y*p2y + p1z*p2z))/(pow(m1,2)*pow(m2,2)) - 
               (11*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2))*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/(pow(m1,2)*pow(m2,2))))/r12 + 
          (m1*m2*((-3*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2)*pow(-(n12x*p2x) - n12y*p2y - n12z*p2z,2))/(pow(m1,2)*pow(m2,2)) - 
               (12*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(-(n12x*p2x) - n12y*p2y - n12z*p2z)*(p1x*p2x + p1y*p2y + p1z*p2z))/(pow(m1,2)*pow(m2,2)) + 
               (10*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2)*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/(pow(m1,2)*pow(m2,2)) - 
               (11*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2))*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/(pow(m1,2)*pow(m2,2)) + 
               (10*pow(pow(p2x,2) + pow(p2y,2) + pow(p2z,2),2))/pow(m2,4) - 
               (2*pow(pow(p2x,2) + pow(p2y,2) + pow(p2z,2),2))/(pow(m1,2)*pow(m2,2))))/r12 + 
          (m1*m3*((10*pow(pow(p1x,2) + pow(p1y,2) + pow(p1z,2),2))/pow(m1,4) - 
               (2*pow(pow(p1x,2) + pow(p1y,2) + pow(p1z,2),2))/(pow(m1,2)*pow(m3,2)) - 
               (3*pow(n13x*p1x + n13y*p1y + n13z*p1z,2)*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/(pow(m1,2)*pow(m3,2)) + 
               (10*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2))*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/(pow(m1,2)*pow(m3,2)) - 
               (12*(n13x*p1x + n13y*p1y + n13z*p1z)*(n13x*p3x + n13y*p3y + n13z*p3z)*(p1x*p3x + p1y*p3y + p1z*p3z))/(pow(m1,2)*pow(m3,2)) - 
               (11*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2))*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/(pow(m1,2)*pow(m3,2))))/r13 + 
          (m1*m3*((-3*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2)*pow(-(n13x*p3x) - n13y*p3y - n13z*p3z,2))/(pow(m1,2)*pow(m3,2)) - 
               (12*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z)*(p1x*p3x + p1y*p3y + p1z*p3z))/(pow(m1,2)*pow(m3,2)) + 
               (10*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2)*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/(pow(m1,2)*pow(m3,2)) - 
               (11*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2))*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/(pow(m1,2)*pow(m3,2)) + 
               (10*pow(pow(p3x,2) + pow(p3y,2) + pow(p3z,2),2))/pow(m3,4) - 
               (2*pow(pow(p3x,2) + pow(p3y,2) + pow(p3z,2),2))/(pow(m1,2)*pow(m3,2))))/r13 + 
          (m2*m3*((10*pow(pow(p2x,2) + pow(p2y,2) + pow(p2z,2),2))/pow(m2,4) - 
               (2*pow(pow(p2x,2) + pow(p2y,2) + pow(p2z,2),2))/(pow(m2,2)*pow(m3,2)) - 
               (3*pow(n23x*p2x + n23y*p2y + n23z*p2z,2)*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/(pow(m2,2)*pow(m3,2)) + 
               (10*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2))*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/(pow(m2,2)*pow(m3,2)) - 
               (12*(n23x*p2x + n23y*p2y + n23z*p2z)*(n23x*p3x + n23y*p3y + n23z*p3z)*(p2x*p3x + p2y*p3y + p2z*p3z))/(pow(m2,2)*pow(m3,2)) - 
               (11*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2))*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/(pow(m2,2)*pow(m3,2))))/r23 + 
          (m2*m3*((-3*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2)*pow(-(n23x*p3x) - n23y*p3y - n23z*p3z,2))/(pow(m2,2)*pow(m3,2)) - 
               (12*(-(n23x*p2x) - n23y*p2y - n23z*p2z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z)*(p2x*p3x + p2y*p3y + p2z*p3z))/(pow(m2,2)*pow(m3,2)) + 
               (10*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2)*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/(pow(m2,2)*pow(m3,2)) - 
               (11*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2))*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/(pow(m2,2)*pow(m3,2)) + 
               (10*pow(pow(p3x,2) + pow(p3y,2) + pow(p3z,2),2))/pow(m3,4) - 
               (2*pow(pow(p3x,2) + pow(p3y,2) + pow(p3z,2),2))/(pow(m2,2)*pow(m3,2))))/r23)/16. + 
     (3*(-((pow(m1,3)*m2)/pow(r12,3)) - (m1*pow(m2,3))/pow(r12,3) - (pow(m1,3)*m3)/pow(r13,3) - (m1*pow(m3,3))/pow(r13,3) - 
          (m1*m2*pow(m3,2))/(r12*pow(r13,2)) - (m1*pow(m2,2)*m3)/(pow(r12,2)*r13) - (pow(m2,3)*m3)/pow(r23,3) - 
          (m2*pow(m3,3))/pow(r23,3) - (m1*m2*pow(m3,2))/(r12*pow(r23,2)) - (m1*pow(m2,2)*m3)/(r13*pow(r23,2)) - 
          (pow(m1,2)*m2*m3)/(pow(r12,2)*r23) - (pow(m1,2)*m2*m3)/(pow(r13,2)*r23)))/4. - 
     (3*((2*pow(m1,2)*pow(m2,2))/pow(r12,3) + (2*pow(m1,2)*pow(m3,2))/pow(r13,3) + (pow(m1,2)*m2*m3)/(r12*pow(r13,2)) + 
          (pow(m1,2)*m2*m3)/(pow(r12,2)*r13) + (2*pow(m2,2)*pow(m3,2))/pow(r23,3) + (m1*pow(m2,2)*m3)/(r12*pow(r23,2)) + 
          (m1*m2*pow(m3,2))/(r13*pow(r23,2)) + (m1*pow(m2,2)*m3)/(pow(r12,2)*r23) + (m1*m2*pow(m3,2))/(pow(r13,2)*r23)))/4. + 
     ((pow(m1,2)*m2*((12*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/pow(m1,2) + 
               ((pow(n12x,2) + pow(n12y,2) + pow(n12z,2))*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/pow(m1,2) + 
               (31*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(-(n12x*p2x) - n12y*p2y - n12z*p2z))/(m1*m2) - (50*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + 
               (18*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2)))/pow(r12,2) + 
          (m1*pow(m2,2)*((18*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(n12x*p1x + n12y*p1y + n12z*p1z)*(n12x*p2x + n12y*p2y + n12z*p2z))/(m1*m2) + 
               (12*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/pow(m2,2) + 
               ((pow(n12x,2) + pow(n12y,2) + pow(n12z,2))*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/pow(m2,2) - 
               (50*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + (31*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2)))/pow(r12,2) + 
          (pow(m1,2)*m3*((12*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/pow(m1,2) + 
               ((pow(n13x,2) + pow(n13y,2) + pow(n13z,2))*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/pow(m1,2) + 
               (31*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z))/(m1*m3) - (50*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) + 
               (18*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r13,2) + 
          (m1*pow(m3,2)*((18*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(n13x*p1x + n13y*p1y + n13z*p1z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m1*m3) + 
               (12*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/pow(m3,2) + 
               ((pow(n13x,2) + pow(n13y,2) + pow(n13z,2))*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/pow(m3,2) - 
               (50*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) + (31*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r13,2) + 
          (m1*m2*m3*((18*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(n12x*p1x + n12y*p1y + n12z*p1z)*(n12x*p2x + n12y*p2y + n12z*p2z))/(m1*m2) - 
               (2*pow(n12x*p2x + n12y*p2y + n12z*p2z,2))/pow(m2,2) - (50*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + 
               (14*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) + 
               (14*(n12x*p2x + n12y*p2y + n12z*p2z)*(n12x*p3x + n12y*p3y + n12z*p3z))/(m2*m3) + 
               ((n12x*n13x + n12y*n13y + n12z*n13z)*(n12x*p2x + n12y*p2y + n12z*p2z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m2*m3) + 
               (17*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3)))/(r12*r13) + 
          (m1*m2*m3*((18*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(n13x*p1x + n13y*p1y + n13z*p1z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m1*m3) + 
               ((n12x*n13x + n12y*n13y + n12z*n13z)*(n12x*p2x + n12y*p2y + n12z*p2z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m2*m3) + 
               (14*(n13x*p2x + n13y*p2y + n13z*p2z)*(n13x*p3x + n13y*p3y + n13z*p3z))/(m2*m3) - 
               (2*pow(n13x*p3x + n13y*p3y + n13z*p3z,2))/pow(m3,2) - (50*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) + 
               (17*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + (14*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/(r12*r13) + 
          (pow(m2,2)*m3*((12*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/pow(m2,2) + 
               ((pow(n23x,2) + pow(n23y,2) + pow(n23z,2))*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/pow(m2,2) + 
               (31*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) - 
               (14*(-(n23x*p2x) - n23y*p2y - n23z*p2z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z))/(m2*m3) - (50*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + 
               (18*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r23,2) + 
          (m2*pow(m3,2)*((18*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) - 
               (14*(n23x*p2x + n23y*p2y + n23z*p2z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m2*m3) + 
               (12*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/pow(m3,2) + 
               ((pow(n23x,2) + pow(n23y,2) + pow(n23z,2))*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/pow(m3,2) - 
               (50*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + (31*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/pow(r23,2) + 
          (m1*m2*m3*((-2*pow(-(n12x*p1x) - n12y*p1y - n12z*p1z,2))/pow(m1,2) + (14*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) - 
               (14*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(-(n12x*p2x) - n12y*p2y - n12z*p2z))/(m1*m2) - (50*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + 
               (18*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) + 
               (14*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(-(n12x*p3x) - n12y*p3y - n12z*p3z))/(m1*m3) + 
               ((-(n12x*n23x) - n12y*n23y - n12z*n23z)*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m1*m3) + 
               (17*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3)))/(r12*r23) + 
          (m1*m2*m3*((18*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) + 
               ((-(n12x*n23x) - n12y*n23y - n12z*n23z)*(-(n12x*p1x) - n12y*p1y - n12z*p1z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m1*m3) + 
               (14*(n23x*p1x + n23y*p1y + n23z*p1z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m1*m3) - 
               (14*(n23x*p2x + n23y*p2y + n23z*p2z)*(n23x*p3x + n23y*p3y + n23z*p3z))/(m2*m3) - 
               (2*pow(n23x*p3x + n23y*p3y + n23z*p3z,2))/pow(m3,2) + (17*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) - 
               (50*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + (14*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/(r12*r23) + 
          (m1*m2*m3*((-2*pow(-(n13x*p1x) - n13y*p1y - n13z*p1z,2))/pow(m1,2) + (14*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2)))/pow(m1,2) + 
               (14*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n13x*p2x) - n13y*p2y - n13z*p2z))/(m1*m2) + 
               ((n13x*n23x + n13y*n23y + n13z*n23z)*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n23x*p2x) - n23y*p2y - n23z*p2z))/(m1*m2) + 
               (17*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) - (14*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z))/(m1*m3) - 
               (50*(p1x*p3x + p1y*p3y + p1z*p3z))/(m1*m3) + (18*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/(r13*r23) + 
          (m1*m2*m3*(((n13x*n23x + n13y*n23y + n13z*n23z)*(-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n23x*p2x) - n23y*p2y - n23z*p2z))/(m1*m2) + 
               (14*(-(n23x*p1x) - n23y*p1y - n23z*p1z)*(-(n23x*p2x) - n23y*p2y - n23z*p2z))/(m1*m2) - 
               (2*pow(-(n23x*p2x) - n23y*p2y - n23z*p2z,2))/pow(m2,2) + (17*(p1x*p2x + p1y*p2y + p1z*p2z))/(m1*m2) + 
               (14*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2)))/pow(m2,2) - 
               (14*(-(n23x*p2x) - n23y*p2y - n23z*p2z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z))/(m2*m3) - (50*(p2x*p3x + p2y*p3y + p2z*p3z))/(m2*m3) + 
               (18*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2)))/pow(m3,2)))/(r13*r23))/8. - 
     (3*((2*pow(m1,2)*m2*m3)/(r12*r13*r23) + (2*m1*pow(m2,2)*m3)/(r12*r13*r23) + (2*m1*m2*pow(m3,2))/(r12*r13*r23)))/8. + 
     ((m1*m2*m3*((((n12x + n13x)*p1x + (n12y + n13y)*p1y + (n12z + n13z)*p1z)*((n12x - n23x)*p1x + (n12y - n23y)*p1y + (n12z - n23z)*p1z))/
               pow(m1,2) + (3*((n12x + n13x)*p1x + (n12y + n13y)*p1y + (n12z + n13z)*p1z)*
               ((n12x - n23x)*p2x + (n12y - n23y)*p2y + (n12z - n23z)*p2z))/(m1*m2) - 
               (16*((n12x - n23x)*p1x + (n12y - n23y)*p1y + (n12z - n23z)*p1z)*((n12x + n13x)*p3x + (n12y + n13y)*p3y + (n12z + n13z)*p3z))/(m1*m3) + 
               (8*((n12x + n13x)*p1x + (n12y + n13y)*p1y + (n12z + n13z)*p1z)*((n12x - n23x)*p3x + (n12y - n23y)*p3y + (n12z - n23z)*p3z))/(m1*m3) + 
               (4*((n12x + n13x)*p3x + (n12y + n13y)*p3y + (n12z + n13z)*p3z)*((n12x - n23x)*p3x + (n12y - n23y)*p3y + (n12z - n23z)*p3z))/pow(m3,2)))
          /pow(r12 + r13 + r23,2) + (m1*m2*m3*((4*((-n12x - n13x)*p2x + (-n12y - n13y)*p2y + (-n12z - n13z)*p2z)*
               ((-n13x - n23x)*p2x + (-n13y - n23y)*p2y + (-n13z - n23z)*p2z))/pow(m2,2) - 
               (16*((-n13x - n23x)*p2x + (-n13y - n23y)*p2y + (-n13z - n23z)*p2z)*((-n12x - n13x)*p3x + (-n12y - n13y)*p3y + (-n12z - n13z)*p3z))/
               (m2*m3) + (3*((-n12x - n13x)*p1x + (-n12y - n13y)*p1y + (-n12z - n13z)*p1z)*
               ((-n13x - n23x)*p3x + (-n13y - n23y)*p3y + (-n13z - n23z)*p3z))/(m1*m3) + 
               (8*((-n12x - n13x)*p2x + (-n12y - n13y)*p2y + (-n12z - n13z)*p2z)*((-n13x - n23x)*p3x + (-n13y - n23y)*p3y + (-n13z - n23z)*p3z))/
               (m2*m3) + (((-n12x - n13x)*p3x + (-n12y - n13y)*p3y + (-n12z - n13z)*p3z)*
               ((-n13x - n23x)*p3x + (-n13y - n23y)*p3y + (-n13z - n23z)*p3z))/pow(m3,2)))/pow(r12 + r13 + r23,2) + 
          (m1*m2*m3*((4*((n12x - n23x)*p1x + (n12y - n23y)*p1y + (n12z - n23z)*p1z)*((-n13x - n23x)*p1x + (-n13y - n23y)*p1y + (-n13z - n23z)*p1z))/
               pow(m1,2) - (16*((-n13x - n23x)*p1x + (-n13y - n23y)*p1y + (-n13z - n23z)*p1z)*
               ((n12x - n23x)*p3x + (n12y - n23y)*p3y + (n12z - n23z)*p3z))/(m1*m3) + 
               (8*((n12x - n23x)*p1x + (n12y - n23y)*p1y + (n12z - n23z)*p1z)*((-n13x - n23x)*p3x + (-n13y - n23y)*p3y + (-n13z - n23z)*p3z))/(m1*m3) + 
               (3*((n12x - n23x)*p2x + (n12y - n23y)*p2y + (n12z - n23z)*p2z)*((-n13x - n23x)*p3x + (-n13y - n23y)*p3y + (-n13z - n23z)*p3z))/(m2*m3) + 
               (((n12x - n23x)*p3x + (n12y - n23y)*p3y + (n12z - n23z)*p3z)*((-n13x - n23x)*p3x + (-n13y - n23y)*p3y + (-n13z - n23z)*p3z))/pow(m3,2))
          )/pow(r12 + r13 + r23,2) + (m1*m2*m3*((3*((-n12x - n13x)*p1x + (-n12y - n13y)*p1y + (-n12z - n13z)*p1z)*
               ((-n12x + n23x)*p2x + (-n12y + n23y)*p2y + (-n12z + n23z)*p2z))/(m1*m2) + 
               (((-n12x - n13x)*p2x + (-n12y - n13y)*p2y + (-n12z - n13z)*p2z)*((-n12x + n23x)*p2x + (-n12y + n23y)*p2y + (-n12z + n23z)*p2z))/
               pow(m2,2) + (8*((-n12x + n23x)*p2x + (-n12y + n23y)*p2y + (-n12z + n23z)*p2z)*
               ((-n12x - n13x)*p3x + (-n12y - n13y)*p3y + (-n12z - n13z)*p3z))/(m2*m3) - 
               (16*((-n12x - n13x)*p2x + (-n12y - n13y)*p2y + (-n12z - n13z)*p2z)*((-n12x + n23x)*p3x + (-n12y + n23y)*p3y + (-n12z + n23z)*p3z))/
               (m2*m3) + (4*((-n12x - n13x)*p3x + (-n12y - n13y)*p3y + (-n12z - n13z)*p3z)*
               ((-n12x + n23x)*p3x + (-n12y + n23y)*p3y + (-n12z + n23z)*p3z))/pow(m3,2)))/pow(r12 + r13 + r23,2) + 
          (m1*m2*m3*((((n12x + n13x)*p1x + (n12y + n13y)*p1y + (n12z + n13z)*p1z)*((n13x + n23x)*p1x + (n13y + n23y)*p1y + (n13z + n23z)*p1z))/
               pow(m1,2) - (16*((n13x + n23x)*p1x + (n13y + n23y)*p1y + (n13z + n23z)*p1z)*
               ((n12x + n13x)*p2x + (n12y + n13y)*p2y + (n12z + n13z)*p2z))/(m1*m2) + 
               (8*((n12x + n13x)*p1x + (n12y + n13y)*p1y + (n12z + n13z)*p1z)*((n13x + n23x)*p2x + (n13y + n23y)*p2y + (n13z + n23z)*p2z))/(m1*m2) + 
               (4*((n12x + n13x)*p2x + (n12y + n13y)*p2y + (n12z + n13z)*p2z)*((n13x + n23x)*p2x + (n13y + n23y)*p2y + (n13z + n23z)*p2z))/
               pow(m2,2) + (3*((n12x + n13x)*p1x + (n12y + n13y)*p1y + (n12z + n13z)*p1z)*
               ((n13x + n23x)*p3x + (n13y + n23y)*p3y + (n13z + n23z)*p3z))/(m1*m3)))/pow(r12 + r13 + r23,2) + 
          (m1*m2*m3*((4*((-n12x + n23x)*p1x + (-n12y + n23y)*p1y + (-n12z + n23z)*p1z)*((n13x + n23x)*p1x + (n13y + n23y)*p1y + (n13z + n23z)*p1z))/
               pow(m1,2) + (8*((n13x + n23x)*p1x + (n13y + n23y)*p1y + (n13z + n23z)*p1z)*
               ((-n12x + n23x)*p2x + (-n12y + n23y)*p2y + (-n12z + n23z)*p2z))/(m1*m2) - 
               (16*((-n12x + n23x)*p1x + (-n12y + n23y)*p1y + (-n12z + n23z)*p1z)*((n13x + n23x)*p2x + (n13y + n23y)*p2y + (n13z + n23z)*p2z))/
               (m1*m2) + (((-n12x + n23x)*p2x + (-n12y + n23y)*p2y + (-n12z + n23z)*p2z)*((n13x + n23x)*p2x + (n13y + n23y)*p2y + (n13z + n23z)*p2z))/
               pow(m2,2) + (3*((-n12x + n23x)*p2x + (-n12y + n23y)*p2y + (-n12z + n23z)*p2z)*
               ((n13x + n23x)*p3x + (n13y + n23y)*p3y + (n13z + n23z)*p3z))/(m2*m3)))/pow(r12 + r13 + r23,2))/2. + 
     ((m1*m2*m3*((-3*(p1x*p2x + p1y*p2y + p1z*p2z - (-(n12x*p1x) - n12y*p1y - n12z*p1z)*(-(n12x*p2x) - n12y*p2y - n12z*p2z)))/(m1*m2) - 
               (pow(p2x,2) + pow(p2y,2) + pow(p2z,2) - pow(-(n12x*p2x) - n12y*p2y - n12z*p2z,2))/pow(m2,2) + 
               (8*(p2x*p3x + p2y*p3y + p2z*p3z - (-(n12x*p2x) - n12y*p2y - n12z*p2z)*(-(n12x*p3x) - n12y*p3y - n12z*p3z)))/(m2*m3) - 
               (4*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2) - pow(-(n12x*p3x) - n12y*p3y - n12z*p3z,2)))/pow(m3,2)))/(r12*(r12 + r13 + r23)) + 
          (m1*m2*m3*(-((pow(p1x,2) + pow(p1y,2) + pow(p1z,2) - pow(n12x*p1x + n12y*p1y + n12z*p1z,2))/pow(m1,2)) - 
               (3*(p1x*p2x + p1y*p2y + p1z*p2z - (n12x*p1x + n12y*p1y + n12z*p1z)*(n12x*p2x + n12y*p2y + n12z*p2z)))/(m1*m2) + 
               (8*(p1x*p3x + p1y*p3y + p1z*p3z - (n12x*p1x + n12y*p1y + n12z*p1z)*(n12x*p3x + n12y*p3y + n12z*p3z)))/(m1*m3) - 
               (4*(pow(p3x,2) + pow(p3y,2) + pow(p3z,2) - pow(n12x*p3x + n12y*p3y + n12z*p3z,2)))/pow(m3,2)))/(r12*(r12 + r13 + r23)) + 
          (m1*m2*m3*((-4*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2) - pow(-(n13x*p2x) - n13y*p2y - n13z*p2z,2)))/pow(m2,2) - 
               (3*(p1x*p3x + p1y*p3y + p1z*p3z - (-(n13x*p1x) - n13y*p1y - n13z*p1z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z)))/(m1*m3) + 
               (8*(p2x*p3x + p2y*p3y + p2z*p3z - (-(n13x*p2x) - n13y*p2y - n13z*p2z)*(-(n13x*p3x) - n13y*p3y - n13z*p3z)))/(m2*m3) - 
               (pow(p3x,2) + pow(p3y,2) + pow(p3z,2) - pow(-(n13x*p3x) - n13y*p3y - n13z*p3z,2))/pow(m3,2)))/(r13*(r12 + r13 + r23)) + 
          (m1*m2*m3*(-((pow(p1x,2) + pow(p1y,2) + pow(p1z,2) - pow(n13x*p1x + n13y*p1y + n13z*p1z,2))/pow(m1,2)) + 
               (8*(p1x*p2x + p1y*p2y + p1z*p2z - (n13x*p1x + n13y*p1y + n13z*p1z)*(n13x*p2x + n13y*p2y + n13z*p2z)))/(m1*m2) - 
               (4*(pow(p2x,2) + pow(p2y,2) + pow(p2z,2) - pow(n13x*p2x + n13y*p2y + n13z*p2z,2)))/pow(m2,2) - 
               (3*(p1x*p3x + p1y*p3y + p1z*p3z - (n13x*p1x + n13y*p1y + n13z*p1z)*(n13x*p3x + n13y*p3y + n13z*p3z)))/(m1*m3)))/(r13*(r12 + r13 + r23))\
          + (m1*m2*m3*((-4*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2) - pow(-(n23x*p1x) - n23y*p1y - n23z*p1z,2)))/pow(m1,2) + 
               (8*(p1x*p3x + p1y*p3y + p1z*p3z - (-(n23x*p1x) - n23y*p1y - n23z*p1z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z)))/(m1*m3) - 
               (3*(p2x*p3x + p2y*p3y + p2z*p3z - (-(n23x*p2x) - n23y*p2y - n23z*p2z)*(-(n23x*p3x) - n23y*p3y - n23z*p3z)))/(m2*m3) - 
               (pow(p3x,2) + pow(p3y,2) + pow(p3z,2) - pow(-(n23x*p3x) - n23y*p3y - n23z*p3z,2))/pow(m3,2)))/(r23*(r12 + r13 + r23)) + 
          (m1*m2*m3*((-4*(pow(p1x,2) + pow(p1y,2) + pow(p1z,2) - pow(n23x*p1x + n23y*p1y + n23z*p1z,2)))/pow(m1,2) + 
               (8*(p1x*p2x + p1y*p2y + p1z*p2z - (n23x*p1x + n23y*p1y + n23z*p1z)*(n23x*p2x + n23y*p2y + n23z*p2z)))/(m1*m2) - 
               (pow(p2x,2) + pow(p2y,2) + pow(p2z,2) - pow(n23x*p2x + n23y*p2y + n23z*p2z,2))/pow(m2,2) - 
               (3*(p2x*p3x + p2y*p3y + p2z*p3z - (n23x*p2x + n23y*p2y + n23z*p2z)*(n23x*p3x + n23y*p3y + n23z*p3z)))/(m2*m3)))/(r23*(r12 + r13 + r23)))/2. + 
     (-((m1*pow(m2,2)*m3*(6*pow(r12,4) + 56*pow(r12,3)*r13 - 60*pow(r12,2)*pow(r13,2) - 72*r12*pow(r13,3) + 35*pow(r13,4) + 
               60*r12*pow(r13,2)*r23 - 24*pow(r12,2)*(r12 + r13)*r23 + 18*pow(r12,2)*pow(r23,2)))/(pow(r12,3)*r13*pow(r23,3))) - 
          (m1*m2*pow(m3,2)*(35*pow(r12,4) - 72*pow(r12,3)*r13 - 60*pow(r12,2)*pow(r13,2) + 56*r12*pow(r13,3) + 6*pow(r13,4) + 
               60*pow(r12,2)*r13*r23 - 24*pow(r13,2)*(r12 + r13)*r23 + 18*pow(r13,2)*pow(r23,2)))/(r12*pow(r13,3)*pow(r23,3)) - 
          (pow(m1,2)*m2*m3*(6*pow(r12,4) + 18*pow(r12,2)*pow(r13,2) + 56*pow(r12,3)*r23 - 60*pow(r12,2)*pow(r23,2) + 
               60*r12*r13*pow(r23,2) - 72*r12*pow(r23,3) + 35*pow(r23,4) - 24*pow(r12,2)*r13*(r12 + r23)))/(pow(r12,3)*pow(r13,3)*r23) - 
          (m1*m2*pow(m3,2)*(35*pow(r12,4) - 72*pow(r12,3)*r23 + 60*pow(r12,2)*r13*r23 - 60*pow(r12,2)*pow(r23,2) + 
               18*pow(r13,2)*pow(r23,2) + 56*r12*pow(r23,3) + 6*pow(r23,4) - 24*r13*pow(r23,2)*(r12 + r23)))/(r12*pow(r13,3)*pow(r23,3))\
          - (pow(m1,2)*m2*m3*(18*pow(r12,2)*pow(r13,2) + 6*pow(r13,4) + 56*pow(r13,3)*r23 + 60*r12*r13*pow(r23,2) - 
               60*pow(r13,2)*pow(r23,2) - 72*r13*pow(r23,3) + 35*pow(r23,4) - 24*r12*pow(r13,2)*(r13 + r23)))/(pow(r12,3)*pow(r13,3)*r23)\
          - (m1*pow(m2,2)*m3*(35*pow(r13,4) + 60*r12*pow(r13,2)*r23 - 72*pow(r13,3)*r23 + 18*pow(r12,2)*pow(r23,2) - 
               60*pow(r13,2)*pow(r23,2) + 56*r13*pow(r23,3) + 6*pow(r23,4) - 24*r12*pow(r23,2)*(r13 + r23)))/(pow(r12,3)*r13*pow(r23,3)))/64.;
}


void update_eom_hamiltonian(double *w, double *dwdt, double (*hamiltonian)(double*, struct ode_params*), double h, struct ode_params* params) {
     int array_half = params->dim * params->num_bodies;
     double w_copy[18];
     double dHdw[18];
     double forward_value, backward_value;

     // Copy original array to w_copy
     for (int i = 0; i < 2 * array_half; i++) {
          w_copy[i] = w[i];
     }

     for (int i = 0; i < 2 * array_half; i++) {
          // Perturb the i-th variable forward
          w_copy[i] = w[i] + h;
          forward_value = hamiltonian(w_copy, params);

          // Perturb the i-th variable backward
          w_copy[i] = w[i] - h;
          backward_value = hamiltonian(w_copy, params);

          // Restore original value
          w_copy[i] = w[i];

          // Compute the partial derivative
          dHdw[i] = (forward_value - backward_value) / (2 * h);
     }

     // Compute dwdt
     for (int i = 0; i < 2 * array_half; i++) {
          if (i < array_half)
               dwdt[i] += dHdw[i + array_half];
          else
               dwdt[i] += -dHdw[i - array_half];
     }
}