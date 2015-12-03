/**
 *  @file    cuajcparm.c
 *  @ingroup CUjacparm
 *  @author  Juan M. Brandi
 *  @brief   Class CUjacparm methods
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

#include "cujacparm.h"

VEMBED(rcsid = "$Id")

#if !defined(VINLINE_CUJACPARM)
#endif

VPUBLIC void CUjacparm_setCenterX(CUjacparm *thee, double x){
		VASSERT(thee != VNULL);
		thee->center[0] = x;
}

VPUBLIC void CUjacparm_setCenterY(CUjacparm *thee, double y){
	VASSERT(thee != VNULL);
	thee->center[1] = y;
}

VPUBLIC void CUjacparm_setCenterZ(CUjacparm *thee, double z){
	VASSERT(thee != VNULL);
	thee->center[2] = z;
}

VPUBLIC double CUjacparm_getCenterX(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->center[0];
}

VPUBLIC double CUjacparm_getCenterY(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->center[1];
}

VPUBLIC double CUjacparm_getCetnerZ(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->center[2];
}

VPUBLIC int CUjacparm_getNx(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->dime[0];
}

VPUBLIC int CUjacparm_getNy(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->dime[1];
}

VPUBLIC int CUjacparm_getNz(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->dime[2];
}

VPUBLIC double CUjacparm_getHx(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->grid[0];
}

VPUBLIC double CUjacparm_getHy(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->grid[1];
}

VPUBLIC double CUjacparm_getHz(CUjacparm *thee){
	VASSERT(thee != VNULL);
	return thee->grid[2];
}

VPUBLIC CUjacparm* CUjacparm_ctor(CUjacparm_CalcType type){
	/*set up the structure*/
	CUjacparm *thee = VNULL;
	thee = (CUjacparm*)Vmem_malloc(VNULL, 1, sizeof(CUjacparm));
	VASSERT(thee != VNULL);
	VASSERT(CUjacparm_ctor2(thee, type) == VRC_SUCCESS);

	return thee;
}

VPUBLIC Vrc_Codes CUjacparm_ctor2(CUjacparm *thee, CUjacparm_CalcType type){

	int i;
	if(thee == VNULL)
		return VRC_FAILURE;

	for(i=0; i<3; i++){
			thee->dime[i] = -1;
			thee->pdime[i] = 1;
	}

	thee->parsed = 0;
	thee->type = type;

	/*generic parameters*/
	thee->etol = 1.0e-6;
	thee->setetol = 0;
	thee->setgrid = 0;
	thee->setgcent = 0;
	thee->setglen = 0;

	/*type 1 & 2 parameters*/
	thee->setcglen = 0;
	thee->setfglen = 0;
	thee->setcgcent = 0;
	thee->setfgcent = 0;

	/*type 2 parameters*/
	thee->setpdime = 0;
	thee->setrank = 0;
	thee->setsize = 0;
	thee->setofrac = 0;
	for(i=0;i<6; i++){
		thee->partDisjOwnSide[i] = 0;
	}
	thee->setasync = 0;

	/*default parameters for TINKER*/
	thee->chgs = VCM_CHARGE;

	thee->useAqua = 0;
	thee->setUseAqua = 0;

	return VRC_SUCCESS;

}

VPUBLIC void CUjacparm_dtor(CUjacparm **thee){
		if((*thee) != VNULL){
				CUjacparm_ctor2(*thee);
				Vmem_free(VNULL, 1, sizeof(CUjacparm), (void **)thee);
				(*thee) = VNULL;
		}
}

VPUBLIC void CUjacparm_dtor2(CUjacparm *thee){ ; }


VPUBLIC Vrc_Codes CUjacparm_check(CUjacparm *thee){

	Vrc_Codes rc;
	int i, tdime[3], ti;

	rc = VRC_SUCCESS;

	Vnm_print(0, "CUjacparm_check: checking CUjacparm of tupe %d.\n",
			thee->type);

	/*check to see if the structure is filled*/
	if(!thee->parsed){
		Vnm_print(2, "CUjacparm_check: not filled!\n");
		return VRC_FAILURE;
	}

	/*check for generic settings*/
	if(!thee->setdime){
		Vnm_print(2,"CUjacparm_check: DIME not set!\n");
		rc = VRC_FAILURE;
	}
	if(!thee->setchgm){
		Vnm_print(2, "CUjacparm_check: CHGM not set!");
		rc = VRC_FAILURE;
	}

	/*check sequential and other settings (probalby will need to change*/
	/*some of this checks                                              */
	if((thee->type == CCT_AMGX) || (thee->type == CCT_DUMMY)){
		if((!thee->setgrid) && (!thee->setglen)){
			Vnm_print(2, "CUjacparm_check: Neither GRID nor GLEN set!\n");
			rc = VRC_FAILURE;
		}
		if((thee->setgrid) && (thee->setglen)){
			Vnm_print(2, "CUjacparm_check: Both GRID and GLEN set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setgcent){
			Vnm_print(2, "CUjacparm_check: GCENT not set!\n");
			rc = VRC_FAILURE;
		}
	}

	/*check sequential an parallel auto settings*/
	if((thee->type == CCT_AUTO) || (thee->type == CCT_OTHER)){
		if(!thee->setcglen){
			Vnm_print(2, "CUjacparm_check: CGLEN not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setfglen){
			Vnm_print(2, "CUjacparm_check: FGLEN not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setcgcent){
			Vnm_print(2, "CUjacparm_check: CGCENT not set!\n");
			rc = VRC_FAILURE;
		}
	}

	/*check parallel and focusing settings*/
	if(thee->type == CCT_PARALLEL){
		if(!thee->setpdime){
			Vnm_print(2, "CUjacparm_check: PDIME not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setrank){
			Vnm_print(2, "CUjacparm_check: PROC_RANK not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setofrac){
			Vnm_print(2, "CUjacparm_check: OFRAC not set!\n");
			rc = VRC_FAILURE;
		}
	}

	/*TODO: may need some kind of check that the provided dime settings
	 * make sense. For now, I will just check that they are positive and
	 * greater than 65 if not fix them*/
	if(thee->type != CCT_DUMMY){
		for(i=0; i<3; i++){
			if(thee->dime[i] < 65){
				thee->dime[i] = 65;
			}
		}
	}
	else{ /*this is a dummy calculation, but still need positive numbers*/
		for(i=0; i<3; i++){
			if(thee->dime[i] <= 0){
				Vnm_print(2, "NOsh: Resetting dime[%d] from %d to 3.\n",
						i, thee->dime[i]);
				thee->dime[i] = 3;
			}
		}
	}

	if(!thee->setUseAqua)
		thee->useAqua = 0;

	return rc;
}

VPUBLIC Vrc_Codes CUjacparm_parseToken(CUjacparm *thee, char tok[VMAX_BUFSIZE],
		Vio *sock){

		if(thee == VNULL){
			Vnm_print(2, "parseCU: got NULL thee!\n");
			return VRC_WARNING;
		}
		if(sock == VNULL){
			Vnm_print(2, "parseCU: got Null socket!\n");
			return VRC_WARNING;
		}

		Vnm_print(0, "CUjacparm_parseToken: trying %s...\n", tok);

		if(Vstring_strcasecmp(tok, "dime") == 0){
			return CUjacparm_parseDime(thee, sock);
		} else if (Vstring_strcasecmp(tok, "chgm") == 0){
			return CUjacparm_parseCHGM(thee, sock);
		} else if(Vstring_strcasecmp(tok, "etol") == 0){
			return CUjacparm_parseETOL(thee, sock);
		} else if(Vstring_strcasecmp(tok, "grid") == 0){
			return CUjacparm_parseGRID(thee, sock);
		} else if(Vstring_strcasecmp(tok, "glen") == 0){
			return CUjackparm_parseGLEN(thee, sock);
		} else if(Vstring_strcasecmp(tok, "gcent") == 0){
			return CUjackparm_parseGCENT(thee, sock);
		} else if(Vstring_strcasecmp(tok, "cglen") == 0){
			return CUjackparm_parseCGLEN(thee, sock);
		} else if(Vstring_strcasecmp(tok, "fglen") == 0){
			return CUjackparm_parseFGLEN(thee, sock);
		} else if(Vstring_strcasecmp(tok, "cgcent") == 0){
			return CUjacparm_parseCGCENT(thee, sock);
		} else if(Vstring_strcasecmp(tok, "fgcent") == 0){
			return CUjacparm_parseFGCENT(thee, sock);
		} else if(Vstring_strcasecmp(tok, "pdime") == 0){
			return CUjacparm_parsePDIME(thee, sock);
		} else if(Vstring_strcasecmp(tok, "ofrac") == 0){
			return CUjacparm_parseOFRAC(thee, sock);
		} else if(Vstring_strcasecmp(tok, "async") == 0){
			return CUjacparm_parseASYNC(thee, sock);
		} else if(Vstring_strcasecmp(tok, "gamma") == 0){
			return CUjackparm_parseGAMMA(thee, sock);
		} else if(Vstring_strcasecmp(tok, "useaqua") == 0){
			return CUjacparm_parseUSEAQUA(thee, sock);
		} else {
			Vnm_print(2, "parseCU: Unrecogized keyword (%s)!\n", tok);
			return VRC_WARNING;
		}

		return VRC_FAILURE;
}

VPRIVATE Vrc_Codes CUjacparm_parseDIME(CUjacparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	int ti;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", &ti) == 0){
			Vnm_print(2, "parseCU: Read a non-integer (%s) while parsing DIME keyword!\n", tok);
			return VRC_WARNING;
	}
	else{
		thee->dime[0] = ti;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", &ti) == 0){
		Vnm_print(2, "parseCU: Read non-integer (%s) while parsing DIME keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->dime[1] = ti;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", &ti) == 0){
		Vnm_print(2, "parseCU: Read non-integer (%s) while parsing DIME keyword!\n");
		return VRC_WARNING;
	}
	else{
		thee->dime[2] = ti;
	}

	thee->setdime = 1;
	return VRC_SUCCESS;

}

















