/**
 *  @file    sorparm.c
 *  @ingroup SORparm
 *  @author  Juan M Brandi
 *  @brief   Class SORparm methods
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2014 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "sorparm.h"

VEMBED(rcsid="$Id$")

VPUBLIC SORparm* SORparm_ctor(SORparm_CalcType type) {

    /* Set up the structure */
    SORparm *thee = VNULL;
    thee = (SORparm*)Vmem_malloc(VNULL, 1, sizeof(SORparm));
    VASSERT( thee != VNULL);
    VASSERT( SORparm_ctor2(thee, type) == VRC_SUCCESS );

    return thee;
}

VPUBLIC Vrc_Codes SORparm_ctor2(SORparm *thee, SORparm_CalcType type) {

    int i;

    if (thee == VNULL) return VRC_FAILURE;

    for (i=0; i<3; i++) {
        thee->dime[i] = -1;
    }

    thee->parsed = 0;
    thee->type = type;

    /* *** GENERIC PARAMETERS *** */
    thee->setdime = 0;
    thee->setchgm = 0;

    /* *** TYPE 0 PARAMETERS *** */
    thee->etol = 1.0e-6;
    thee->setetol = 0;
    thee->setgrid = 0;
    thee->setglen = 0;
    thee->setgcent = 0;

    /* *** TYPE 1 & 2 PARAMETERS *** */
    thee->setfglen = 0;

    /* *** Default parameters for TINKER *** */
    thee->chgs = VCM_CHARGE;

    thee->useAqua = 0;
    thee->setUseAqua = 0;

    return VRC_SUCCESS;
}

VPUBLIC void SORparm_dtor(SORparm **thee) {
    if ((*thee) != VNULL) {
        SORparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(SORparm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void SORparm_dtor2(SORparm *thee) { ; }

VPUBLIC int SORparm_copyMGparm(MGparm *mgparm, SORparm *sorparm){

	if(mgparm == VNULL || sorparm == VNULL){
		 Vnm_tprint( 2, "SORparm_copyMGparm: received null object.\n");
		 return 0;
	}

	//check for parallel status
	if(mgparm->type == MCT_PARALLEL){
		/**TODO: SOR doesn't support parallel right now*/
		return 0;
	}

	sorparm->parsed = mgparm->parsed;

	int i;
	for(i=0; i<3; i++){
		sorparm->dime[i] = mgparm->dime[i];
	}

	//generic parameters
	sorparm->setdime = mgparm->setdime;
	sorparm->setchgm = mgparm->setchgm;

	/* *** TYPE 0 PARAMETERS *** */
	sorparm->etol = mgparm->etol;
	sorparm->setetol = mgparm->setetol;
	sorparm->setgrid = mgparm->setgrid;
	sorparm->setglen = mgparm->setglen;
	sorparm->setgcent = mgparm->setgcent;

	/* *** TYPE 1 & 2 PARAMETERS *** */
	sorparm->setfglen = mgparm->setfglen;
	sorparm->fcmeth = mgparm->fcmeth;
	sorparm->fcentmol = mgparm->fcentmol;

	/* *** Default parameters for TINKER *** */
	sorparm->chgs = mgparm->chgs;

	sorparm->useAqua = mgparm->useAqua;
	sorparm->setUseAqua = mgparm->setUseAqua;

	return 1;
}


