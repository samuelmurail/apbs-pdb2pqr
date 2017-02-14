/**
 *  @file    gpuparm.c
 *  @ingroup GPUparm
 *  @author  Nathan Baker
 *  @brief   Class GPUparm methods
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

#include "gpuparm.h"

VEMBED(rcsid="$Id")

#if !defined(VINLINE_GPUPARM)

#endif /* if !defined(VINLINE_GPUPARM) */

VPUBLIC GPUparm* GPUparm_ctor(MGparm_CalcType type){

	/* Set up the structure */
	GPUparm *thee = VNULL;
	thee = (GPUparm*)Vmem_malloc(VNULL, 1, sizeof(GPUparm));
	VASSERT(thee != VNULL);
	VASSERT(GPUparm_ctor2(thee, type) == VRC_SUCCESS);

	return thee;
}

VPUBLIC Vrc_Codes GPUparm_ctor2(GPUparm *thee, MGparm_CalcType type){

	//TODO: add code to get number of gpu, devices, sizes, etc..

	thee->mgparm = MGparm_ctor(type);
	return VRC_SUCCESS;
}

VPUBLIC void GPUparm_dtor(GPUparm **thee){
	if((*thee) != VNULL){
		MGparm_dtor(&((*thee)->mgparm));
		GPUparm_dtor2(*thee);
		Vmem_free(VNULL, 1, sizeof(GPUparm), (void**)thee);
		(*thee) = VNULL;
	}
}

VPUBLIC void GPUparm_dtor2(GPUparm *thee){ ; };

VPUBLIC Vrc_Codes GPUparm_parseToken(GPUparm *thee, char tok[VMAX_BUFSIZE], Vio *sock){

	if(thee == VNULL){
		Vnm_print(2, "parseGPU: got NULL thee!\n");
		return VRC_WARNING;
	}
	if(sock == VNULL){
		Vnm_print(2, "parseGPU: got NULL socket!\n");
		return VRC_WARNING;
	}

	Vnm_print(0, "GPUparm_parseToken: trying %s...\n", tok);

	return 1;//MGparm_parseToken(thee->mgparm, tok, sock);
}

VPUBLIC void GPUparm_copy(GPUparm *thee, GPUparm *parm){

	MGparm_copy(thee->mgparm, parm->mgparm);
}
