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
 * spg_cpu.h
 *
 *  Created on: 4 October 2019
 *      Author: Timothy Spain, t.spain@ucl.ac.uk
 */

#ifndef SPG_CPU_H
#define SPG_CPU_H

#include "helper_timer.h"

class spg_cpu
{
    // Array dimensions
    static const int npix;
    static const int nz;
    static const int nframes;
    static const int nlos;
    static const int ncoeff;
    static const int nwavcoeff;

    // Pointers to the data arrays
    //float * x; // Not needed due to not copying data to the GPU
    float *u;
    float *u_pos;
    float *w;

    // Pre-conditioning matrix
    float *p;
    float *pp;

    static const float epsilon = 1e-4;
    static const float epsilon_0 = 1e-5; // Hard-coded in the original (spg.cu L64)

    StopWatchInterface *timer;

public:
    /* Constructor
     *
     */
    spg_cpu( int npix, int nz, int nframes, const double *P, const float *l1_weights );

    /* Destructor
     *
     */
    ~spg_cpu();

    /* Compute the proximity operator of the sparsity constraint
     *
     */
    void prox_l1(float *alpha, int niter=10000);

    /* Compute the proximity operator of the positivity constraint
     *
     */
    void prox_pos(float *delta, int niter=10000);

    /* Updates the l1 thresholds
     *
     */
    void update_weights(float *l1_weights);

    // Apply/retrieve the u arrays
    void inject_u_pos(float *h_u_pos);
    void extract_u_pos(float *h_u_pos);

    void inject_u(float *h_u);
    void extract_u(float *h_u);

private:
    // Iteration member function
    void iterate_prox_pos(int niter, float px[], float u[], float gg0[], float epsilon_0, float epsilon);

};
#endif /* SPG_CPU_H */
