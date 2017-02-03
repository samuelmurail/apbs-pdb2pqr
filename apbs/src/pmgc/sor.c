/**
 *  @ingroup PMGC
 *  @author  Juan Brandi
 *  @brief
 *  @version $Id:
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2014 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "sor.h"

VPUBLIC void Vsor(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint){

	int numdia; /*number of non-zero diags in matrix depending on
				 * method (finite or volume element).
				 */

	MAT2(ac, *nx * *ny * *nz, 1);

	//Do in one step
	numdia = VAT(ipc, 11);
	if (numdia == 7) {
		Vsor7x(nx, ny, nz,
                ipc, rpc,
                RAT2(ac, 1,1), cc, fc,
                RAT2(ac, 1,2), RAT2(ac, 1,3), RAT2(ac, 1,4),
                x, w1, w2, r,
                itmax, iters, errtol, omega, iresid, iadjoint);
	} else if (numdia == 27) {
		Vsor27x(nx, ny, nz,
				 ipc, rpc,
				 RAT2(ac, 1, 1), cc, fc,
				 RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
				 RAT2(ac, 1, 5), RAT2(ac, 1, 6),
				 RAT2(ac, 1, 7), RAT2(ac, 1, 8), RAT2(ac, 1, 9), RAT2(ac, 1,10),
				 RAT2(ac, 1,11), RAT2(ac, 1,12), RAT2(ac, 1,13), RAT2(ac, 1,14),
				 x, w1, w2, r,
				 itmax, iters, errtol, omega, iresid, iadjoint);
	} else {
		Vnm_print(2, "SOR: invalid stencil type given...\n");
	}
}

VPUBLIC void Vsor7x(int *nx,int *ny,int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *fc,
        double *oE, double *oN, double *uC,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint){

    int i, j, k, ioff;

    MAT3(cc, *nx, *ny, *nz);
    MAT3(fc, *nx, *ny, *nz);
    MAT3( x, *nx, *ny, *nz);
    MAT3(w1, *nx, *ny, *nz);
    MAT3(w2, *nx, *ny, *nz);
    MAT3( r, *nx, *ny, *nz);

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);

    for (*iters=1; *iters<=*itmax; (*iters)++) {

		for (k=2; k<=*nz-1; k++) {
			for (j=2; j<=*ny-1; j++) {
				for (i=2; i<=*nx-1; i++) {
							VAT3(x, i, j, k) = (1-*omega)*VAT3(x, i, j ,k)+ (*omega)*(
									VAT3(fc,   i,  j,  k)
								 +  VAT3(oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
								 +  VAT3(oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
								 +  VAT3(oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
								 +  VAT3(oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
								 + VAT3( uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
								 + VAT3( uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
								 ) / (VAT3(oC, i, j, k) + VAT3(cc, i, j, k));
						}
					}
				}

    }

    if (*iresid == 1)
        Vmresid7_1s(nx, ny, nz, ipc, rpc, oC, cc, fc, oE, oN, uC, x, r);
}

VPUBLIC void Vsor27x(int *nx,int *ny,int *nz,
        int *ipc, double *rpc,
        double  *oC, double  *cc, double  *fc,
        double  *oE, double  *oN, double  *uC, double *oNE, double *oNW,
        double  *uE, double  *uW, double  *uN, double  *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint){

	Vnm_print(2, "SOR:  27x diag not currently supported.\n");

}
