/*! Copyright UCL, 2019
 * author : Timothy Spain < t.spain@ucl.ac.uk >
 *
 * This software is a computer program whose purpose is to reconstruct mass maps
 * from weak gravitational lensing.
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 * spg_cpu.cpp
 *
 *  Created on: 2019-10-04
 */

#include <cmath>
#include <limits>
#include <cstring>
#include "spg_cpu.h"

void reduction(int nz, int nlos, float result[], float left_curr[], float left_old[], float right_curr[], float right_old[]);


spg_cpu::spg_cpu( int npix, int nz, int nframes, const double *P, const float *l1_weights ) :
npix ( npix ), nz ( nz ), nframes ( nframes ),
nlos ( npix * npix ), ncoeff ( nlos * nz ),
nwavelets (nlos * nframes), nwavcoeff ( nwavelets * nz )
{

    // Allocate arrays for the calculation variables, and zero the arrays
    u = new float[nwavcoeff]();
    u_pos = new float[ncoeff]();
    w = new float[nwavcoeff]();

    // Allocate arrays for and store the preconditioning matrix
    pp = new float[nz*nz]();
    p = new float[nz*nz];

    // copy the preconditioning data
    for (int zz = 0; zz < nz * nz; zz++) {
        p[zz] = P[zz];
    }

    // Calculate the squared matrix
    for (int z1 = 0; z1 < nz; z1++) {
        for (int z2 = 0; z2 < nz; z2++) {
            for (int z3 = 0; z3 < nz; z3++) {
                pp[z1 * nz + z2] += p[z1 * nz + z3] * p[z3 * nz + z2];
            }
        }
    }

    // Initialize the timer
    timer = NULL;
    sdkCreateTimer( &timer );
    sdkResetTimer( &timer );
}

spg_cpu::~spg_cpu()
{
    // Free allocated arrays
    delete[] u;
    delete[] u_pos;
    delete[] w;

    delete[] pp;
    delete[] p;

    StopWatchInterface *timer;
}

void spg_cpu::prox_pos ( float *delta, int niter )
{
    float conditioned[ncoeff];
    float condsq[ncoeff];

    float gg0[nlos];

    // zero fill the arrays
    for ( int i = 0; i < ncoeff; i++ ) {
        conditioned[i] = 0.;
    }
    for ( int x = 0; x < nlos; x++ ) {
        gg0[x] = 0;
    }

    // Multiply by the preconditioning matrix (spg.cu LL282—287)
    for ( int z1 = 0; z1 < nz; z1++ ) {
        for ( int z2 = 0; z2 < nz; z2++ ) {
            int pindex = z2*nz + z1;
            for ( int x = 0; x < nlos; x++ ) {
                int dindex = z2*nz + x;
                int cindex = z1*nz + x;
                conditioned[cindex] += delta[dindex] * p[pindex];
            }
        }
    }

    // Square the matrix element-by-element (spg.cu L290)
    for ( int z = 0; z < nz; z++ ){
        for ( int x = 0; x < nlos; x++ ) {
            int cindex = z*nz + x;
            condsq[cindex] = conditioned[cindex] * conditioned[cindex];
        }
    }
    // Sum along the redshift dimension (spg.cu LL292—297)
    for ( int z = 0; z < nz; z++ ){
        for ( int x = 0; x < nlos; x++ ) {
            int cindex = z*nz + x;
            gg0[x] += condsq[cindex];
        }
    }

    // Calculate gg0 (spg.cu L299)
    for ( int x = 0; x < nlos; x++ ) {
        gg0[x] = std::sqrt(gg0[x]);
    }

    iterate_prox_pos(niter, conditioned, u_pos, gg0, epsilon_0, epsilon);

    // Compute transformation of the resulting array (spg.cu LL417—430)
    for ( int z1 = 0; z1 < nz; z1++ ) {
        for ( int x = 0; x < nlos; x++ ) {
            int idx = z1 * nlos + x;
            delta[idx] = 0.;
            for ( int z2 = 0; z2 < nz; z2++ ) {
                int idx2 = z2 * nlos + x;
                if (u[idx2] < 0.) {
                    delta[idx] += u[idx2] * p[z2 * nz + z1];
                }
            }
        }
    }
}

void spg_cpu::iterate_prox_pos(int niter, float px[], float u[], float gg0[], float epsilon_0, float epsilon) {
    float g[ncoeff];
    float gold[ncoeff];

    float uold[ncoeff];

    float sy[nlos];
    float bb[nlos];

    float optim[ncoeff];

    for ( int iter = 0; iter < niter; iter++){
        // Zero fill the arrays
        for ( int i = 0; i < ncoeff; i++ ) {
            // Not g
            // Not u
            // Not gold
            // Not uold
            optim[i] = 0.;
        }
        for ( int x = 0; x < nlos; x++ ) {
            sy[x] = 0.;
            bb[x] = 0.;
        }
        // Initialize the gradient with A^t x (spg.cu L304)
        for ( int i = 0; i < ncoeff; i++ ) {
            g[i] = -px[i];
        }

        // Compute the matrix product which gives g[z] (spg.cu L314)
        for ( int z = 0; z < nz; z++ ) {
            for ( int x = 0; x < nlos; x++ ) {
                for ( int z1 = 0; z1 < nz; z1++) {
                    g[z * nz + x] += pp[z1 * nz + z] * u[z1 * nlos + x];
                }
            }
        }

        // Compute BB steps (spg.cu LL318—359)
        if ( iter > 0 ) {
            reduction(nz, nlos, sy, g, gold, u, uold);
            if (iter % 2) {
                reduction(nz, nlos, bb, g, gold, g, gold);
                for ( int x = 0; x < nlos; x++ ) {
                    bb[x] = (bb[x] == 0.) ? 0. : sy[x]/bb[x];
                }
            } else {
                reduction(nz, nlos, bb, u, uold, u, uold);
                for ( int x = 0; x < nlos; x++ ) {
                    bb[x] = (sy[x] == 0.) ? 0. : bb[x]/sy[x];
                }
            }
        }

        // Compute optimality checks (spg.cu LL362—368) and reduce over redshift (LL373—376)
        for ( int z  = 0; z < nz; z++ ) {
            for ( int x = 0; x < nlos; x++ ) {
                int idx = z * nlos + x;
                if ( std::sqrt(gg0[x]) >= epsilon_0 ) {
                    float vv = (u[x] == 0.) ? std::fmax(g[idx], 0.) : g[idx];
                    optim[x] += vv * vv / gg0[x];
                }
            }
        }
        // spg.cu LL381—384
        for ( int x = 0; x < nlos; x++ ) {
            optim[x] = std::sqrt(optim[x]);
        }

        float eps_max = -std::numeric_limits<float>::infinity();
        // Update variables only if optimality has not been achieved (spg.cu LL388-408)
        // Loop over lines of sight
        for ( int x = 0; x < nlos; x++ ) {
            if (optim[x] > epsilon) {
                // Prevent bb from increasing in magnitude without bound
                bb[x] = ( std::fabs(bb[x]) > 1.0 ) ? 0.001 : bb[x];
                // Update uold, gold
                for ( int z = 0; z < nz; z++ ) {
                    int idx = z * nlos + x;
                    uold[idx] = u[idx];
                    gold[idx] = g[idx];

                    u[idx] = u[idx] - bb[x] * g[idx];
                    // Apply projection
                    u[idx] = u[idx] - std::fmax(u[idx], 0.);
                }
            }
            // Reduce epsilon across lines of sight with the maximum function (spg.cu LL403—407)
            eps_max = std::fmax(eps_max, optim[x]);
        }

        // Exit the iteration if optimality has been reached for every l.o.s. (spg.cu LL411—412)
        if ( eps_max <= epsilon )
            break;
    }

}

void reduction(int nz, int nlos, float result[], float left_curr[], float left_old[], float right_curr[], float right_old[]) {
    for ( int z1 = 0; z1 < nz; z1++ ) {
        for ( int x = 0; x < nlos; x++ ) {
            int idx = z1 * nlos + x;
            result[x] += (left_curr[idx] - left_old[idx]) * (right_curr[idx] - right_old[idx]);
        }
    }
}

void spg_cpu::prox_l1 ( float *alpha, int niter )
{
    float conditioned[nwavcoeff];
    float condsq[nwavcoeff];

    float gg0[nwavelets];

    // zero fill the arrays (spg.cu LL67—70)
    for ( int i = 0; i < nwavcoeff; i++ ) {
        conditioned[i] = 0.;
    }
    // (spg.cu L62)
    for ( int l = 0; l < nwavelets; l++ ) {
        gg0[l] = 0;
    }

    // Multiply by the preconditioning matrix (spg.cu LL84—87)
    for ( int z1 = 0; z1 < nz; z1++ ) {
        for ( int z2 = 0; z2 < nz; z2++ ) {
            int pindex = z2*nz + z1;
            for ( int l = 0; l < nwavelets; l++ ) {
                int aindex = z2*nz + l;
                int cindex = z1*nz + l;
                conditioned[cindex] += alpha[aindex] * p[pindex];
            }
        }
    }

    // Square the matrix element-by-element (spg.cu L92)
    for ( int z = 0; z < nz; z++ ){
        for ( int x = 0; x < nlos; x++ ) {
            int cindex = z*nz + x;
            condsq[cindex] = conditioned[cindex] * conditioned[cindex];
        }
    }
    // Sum along the redshift dimension (spg.cu LL94—100)
    for ( int z = 0; z < nz; z++ ){
        for ( int l = 0; l < nwavelets; w++ ) {
            int cindex = z*nz + l;
            gg0[l] += condsq[cindex];
        }
    }

    // Calculate gg0 (spg.cu L101)
    for ( int l = 0; l < nwavelets; l++ ) {
        gg0[l] = std::sqrt(gg0[l]);
    }

    // Iterate (spg.cu LL104—222
    iterate_prox_l1(niter, conditioned, u, gg0, epsilon_0, epsilon_l1);

    // Compute transformation of the resulting array (spg.cu LL226—229)
    for (int i = 0; i < nwavcoeff; i++) {
        conditioned[i] = u[i] - std::copysign( std::fdim( std::fabs( u[i] ), w[i]), u[i] );
    }

    // Compute matrix product (spg.cu LL231—233
    for ( int z1 = 0; z1 < nz; z1++ ) {
        for ( int l = 0; l < nwavelets; l++ ) {
            int idx = z1 * nlos + l;
            alpha[idx] = 0.;
            for ( int z2 = 0; z2 < nz; z2++ ) {
                int idx2 = z2 * nlos + w;
                if (u[idx2] < 0.) {
                    alpha[idx] += conditioned[idx2] * p[z2 * nz + z1];
                }
            }
        }
    }
}

void spg_cpu::iterate_prox_l1(int niter, float px[], float u[], float gg0[], float epsilon_0, float epsilon) {
    float g[nwavcoeff];
    float gold[nwavcoeff];

    float uold[nwavcoeff];

    float sy[nwavelets];
    float bb[nwavelets];

    float optim[nwavcoeff];

    for ( int iter = 0; iter < niter; iter++){
        // Zero fill the arrays
        for ( int i = 0; i < nwavcoeff; i++ ) {
            // Not g
            // Not u
            // Not gold
            // Not uold
            optim[i] = 0.;
        }
        for ( int l = 0; l < nwavelets; l++ ) {
            sy[l] = 0.;
            bb[l] = 0.;
        }
        // Initialize the gradient with A^t x (spg.cu L106)
        for ( int i = 0; i < nwavcoeff; i++ ) {
            g[i] = -px[i];
        }

        // Compute the matrix product which gives g[z] (spg.cu L116)
        for ( int z = 0; z < nz; z++ ) {
            for ( int l = 0; l < nwavelets; l++ ) {
                for ( int z1 = 0; z1 < nz; z1++) {
                    g[z * nz + l] += pp[z1 * nz + z] * u[z1 * nwavelets + l];
                }
            }
        }

        // Compute BB steps (spg.cu LL120—161)
        if ( iter > 0 ) {
            reduction(nz, nwavelets, sy, g, gold, u, uold);
            if (iter % 2) {
                reduction(nz, nwavelets, bb, g, gold, g, gold);
                for ( int l = 0; l < nwavelets; l++ ) {
                    bb[l] = (bb[l] == 0.) ? 0. : sy[l]/bb[l];
                }
            } else {
                reduction(nz, nwavelets, bb, u, uold, u, uold);
                for ( int l = 0; l < nwavelets; l++ ) {
                    bb[l] = (sy[l] == 0.) ? 0. : bb[l]/sy[l];
                }
            }
        }

        // Compute optimality checks (spg.cu LL164—177) and reduce over redshift (LL179—184)
        for ( int z  = 0; z < nz; z++ ) {
            for ( int l = 0; l < nwavelets; l++ ) {
                int idx = z * nwavelets + l;
                if ( std::sqrt(gg0[l]) >= epsilon_0 ) {
                    float vv;
                    const float c = 1 - 1e-5;
                    if ( u[l] <= -c * w )
                        vv = std::fmin(0.f, g[idx]);
                    if ( u[l] >= c * w )
                        vv = std::fmax(0.f, g[idx]);
                    else
                        vv = g[idx];
                    optim[l] += vv * vv / gg0[l];
                }
            }
        }
        // spg.cu LL189—192
        for ( int l = 0; l < nwavelets; l++ ) {
            optim[l] = std::sqrt(optim[l]);
        }

        float eps_max = -std::numeric_limits<float>::infinity();
        // Update variables only if optimality has not been achieved (spg.cu LL196—216)
        // Loop over lines of sight
        for ( int l = 0; l < nwavelets; l++ ) {
            if (optim[l] > epsilon) {
                // Prevent bb from increasing in magnitude without bound
                bb[l] = ( std::fabs(bb[l]) > 1.0 ) ? 0.001 : bb[l];
                // Update uold, gold
                for ( int z = 0; z < nz; z++ ) {
                    int idx = z * nwavelets + l;
                    uold[idx] = u[idx];
                    gold[idx] = g[idx];

                    u[idx] = u[idx] - bb[l] * g[idx];
                    // Apply projection
                    u[idx] = u[idx] - std::copysign(std::fdim(std::fabs(u[idx]), w[idx]), u[idx]);
                }
            }
            // Reduce epsilon across lines of sight with the maximum function (spg.cu LL211—215)
            eps_max = std::fmax(eps_max, optim[l]);
        }

        // Exit the iteration if optimality has been reached for every l.o.s. (spg.cu LL219—221)
        if ( eps_max <= epsilon )
            break;
    }

}

void spg_cpu::update_weights ( float *l1_weights )
{
    std::memcpy(w, l1_weights, nwavcoeff * sizeof(float));
}

void spg_cpu::inject_u_pos ( float *h_u_pos )
{
    std::memcpy(u_pos, h_u_pos, ncoeff * sizeof(float));
}

void spg_cpu::extract_u_pos( float *h_u_pos )
{
    std::memcpy(h_u_pos, u_pos, ncoeff * sizeof(float));
}

void spg_cpu::inject_u( float *h_u )
{
    std::memcpy(u, h_u, nwavcoeff * sizeof(float));
}

void spg_cpu::extract_u( float *h_u )
{
    std::memcpy(h_u, u, nwavcoeff * sizeof(float));
}
