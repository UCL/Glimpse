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

#include "spg_cpu.h"

spg_cpu::spg_cpu( int npix, int nz, int nframes, const double *P, const float *l1_weights ) :
npix ( npix ), nz ( nz ), nframes ( nframes ),
nlos ( npix * npix ), ncoeff ( nlos * nz ), nwavcoeff ( ncoeff * nframes )
{

    // Coefficient plurality
    int nwavcoeff = npix * npix * nframes * nz;

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
    int nlos = npix * npix;
    int ncoeff = nlos * nz;

    float conditioned[] = new float[ncoeff]();
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
    float condsq[] = new float[ncoeff]();
    for ( int z = 0; z < nz; z++ ){
        for ( int x = 0; x < nlos; x++ ) {
            int cindex = z*nz + x;
            condsq[cindex] = conditioned[cindex] * conditioned[cindex];
        }
    }
    // Sum along the redshift dimension (spg.cu LL292—297)
    float gg0[] = new float[nlos]();
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

    iterate_prox_pos(niter, conditioned, u_pos, gg0);
}

void spg_cpu::iterate_prox_pos(int niter, float px[], float u[], float gg0[]) {
    float g[] = new float[ncoeff];
    float gold[] = new float[ncoeff];

    float uold[] = new float[ncoeff];

    for ( int iter = 0; iter < niter; iter++){
        // Initialize the gradient with A^t x (spg.cu L304)
        for ( int i = 0; i < ncoeff; i++ ) {
            g[i] = -px[i];
        }

        // Compute the matrix product which gives g[z] (spg.cu L314)
        for ( int z = 0; z < nz; z++ ) {

        }
    }

}


void spg_cpu::prox_l1 ( float *alpha, int niter )
{
// Do nothing
}

void spg_cpu::update_weights ( float *l1_weights )
{

}

void spg_cpu::inject_u_pos ( float *h_u_pos )
{

}

void spg_cpu::extract_u_pos( float *h_u_pos )
{

}

void spg_cpu::inject_u( float *h_u )
{

}

void spg_cpu::extract_u( float *h_u )
{

}
