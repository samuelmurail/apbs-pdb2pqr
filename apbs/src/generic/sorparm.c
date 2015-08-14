/**
 *  @file    sorparm.c
 *  @ingroup SORparm
 *  @author  Juan Brandi
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

VEMBED(rcsid = "$Id")

#if !definine(VINLINE_SORPARM)

#endif

VPUBLIC SORparm* SORparm_ctor(){

	/*set up structure*/
	SORparm* thee = VNULL;
	thee = (SORparm*)Vmem_malloc(VNULL,1,sizeof(SORparm));
	VASSERT(thee!=VNULL);
	VASSSERT(SORparm_ctor2(thee));

	return thee;
}

VPUBLIC Vrc_Codes SORparm_ctor2(SORparm *thee){

	if(thee==VNULL)
		return VRC_FAILURE;

	int i;
	for(i=0;i<3;i++){
		thee->dim[i]=-1;
	}

	thee->parsed = 0;

	/* *** GENERIC PARAMETERS *** */
	thee->setdime = 0;
	thee->setchgm = 0;

	/* *** TYPE 0 PARAMETERS *** */
	thee->etol = 1.0e-6;
	thee->setetol = 0;
	thee->setgrid = 0;
	thee->setglen = 0;
	thee->setgcent = 0;

	/* *** TYPE 1 PARAMETERS *** */
	thee->setcglen = 0;
	thee->setfglen = 0;

	/* *** DEFAULT PARAMETERS FOR TINKER *** */
	thee->chgm = VCM_CHARGE;
	thee->useAqua = 0;
	thee->setUseAqua = 0;

	return VRC_SUCCESS;
}

VPUBLIC void SORparm_dtor(SORparm **thee){

	if((*thee) != VNULL){
		SORparm_dtor2(*thee);
		Vmem_free(VNULL, 1, sizeof(SORparm), (void **)thee);
		(*thee) = VNULL;
	}
}

VPUBLIC void SORparm_dtor2(SORparm *thee){ ; }

VPUBLIC int SORparm_copyMGparm(NOsh *thee){

	if(thee == VNULL){
		Vnm_print(0,"SORparm_CopyMGparm: received null thee...\n");
		return 0;
	}

	MGparm *mgcalc;
	mgcalc = thee->calc->mgparm;
	SORparm *sorcalc;
	sorcalc = thee->calc->sorparm;

	/* ***Generic Parameters*** */
	int i;
	for(i=0;i<3;i++){
		sorcalc->dim[i] =  (mgcalc->dime)[i];
	}
	sorcalc->setdime = mgcalc->setdime;
	sorcalc->chgm = mgcalc->chgm;
	sorcalc->chgs = mgcalc->chgs;

	/* ***TYPE 0 PARAMETERS*** */
	sorcalc->etol = mgcalc->etol;
	sorcalc->setetol = mgcalc->setetol;
	for(i=0;i<3;i++){
		sorcalc->grid[i] = (mgcalc->grid)[i];
		sorcalc->glen[i] = (mgcalc->glen)[i];
		sorcalc->center[i] = (mgcalc->center)[i];
	}
	sorcalc->setgrid = mgcalc->setgrid;
	sorcalc->setglen = mgcalc->setglen;
	sorcalc->cmeth = mgcalc->cmeth;
	sorcalc->centmol = mgcalc->centmol;
	sorcalc->setgcent = mgcalc->setgcent;

	/* ***TYPE 1 PARAMETERS*** */
	sorcalc->setcglen = mgcalc->setcglen;
	sorcalc->setfglen = mgcalc->setfglen;
	sorcalc->method = mgcalc->method;
	sorcalc->setmethod = mgcalc->setmethod;
	sorcalc->useAqua = mgcalc->useAqua;
	sorcalc->setUseAqua = mgcalc->setUseAqua;


	return 1;
}






























