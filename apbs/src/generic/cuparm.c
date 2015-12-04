/**
 *  @file    cuajcparm.c
 *  @ingroup CUparm
 *  @author  Juan M. Brandi
 *  @brief   Class CUparm methods
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

#include "cuparm.h"

VEMBED(rcsid = "$Id")

#if !defined(VINLINE_CUJACPARM)
#endif

VPUBLIC void CUparm_setCenterX(CUparm *thee, double x){
		VASSERT(thee != VNULL);
		thee->center[0] = x;
}

VPUBLIC void CUparm_setCenterY(CUparm *thee, double y){
	VASSERT(thee != VNULL);
	thee->center[1] = y;
}

VPUBLIC void CUparm_setCenterZ(CUparm *thee, double z){
	VASSERT(thee != VNULL);
	thee->center[2] = z;
}

VPUBLIC double CUparm_getCenterX(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->center[0];
}

VPUBLIC double CUparm_getCenterY(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->center[1];
}

VPUBLIC double CUparm_getCetnerZ(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->center[2];
}

VPUBLIC int CUparm_getNx(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->dime[0];
}

VPUBLIC int CUparm_getNy(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->dime[1];
}

VPUBLIC int CUparm_getNz(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->dime[2];
}

VPUBLIC double CUparm_getHx(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->grid[0];
}

VPUBLIC double CUparm_getHy(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->grid[1];
}

VPUBLIC double CUparm_getHz(CUparm *thee){
	VASSERT(thee != VNULL);
	return thee->grid[2];
}

VPUBLIC CUparm* CUparm_ctor(CUparm_CalcType type){
	/*set up the structure*/
	CUparm *thee = VNULL;
	thee = (CUparm*)Vmem_malloc(VNULL, 1, sizeof(CUparm));
	VASSERT(thee != VNULL);
	VASSERT(CUparm_ctor2(thee, type) == VRC_SUCCESS);

	return thee;
}

VPUBLIC Vrc_Codes CUparm_ctor2(CUparm *thee, CUparm_CalcType type){

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

VPUBLIC void CUparm_dtor(CUparm **thee){
		if((*thee) != VNULL){
				CUparm_ctor2(*thee);
				Vmem_free(VNULL, 1, sizeof(CUparm), (void **)thee);
				(*thee) = VNULL;
		}
}

VPUBLIC void CUparm_dtor2(CUparm *thee){ ; }


VPUBLIC Vrc_Codes CUparm_check(CUparm *thee){

	Vrc_Codes rc;
	int i, tdime[3], ti;

	rc = VRC_SUCCESS;

	Vnm_print(0, "CUparm_check: checking CUparm of tupe %d.\n",
			thee->type);

	/*check to see if the structure is filled*/
	if(!thee->parsed){
		Vnm_print(2, "CUparm_check: not filled!\n");
		return VRC_FAILURE;
	}

	/*check for generic settings*/
	if(!thee->setdime){
		Vnm_print(2,"CUparm_check: DIME not set!\n");
		rc = VRC_FAILURE;
	}
	if(!thee->setchgm){
		Vnm_print(2, "CUparm_check: CHGM not set!");
		rc = VRC_FAILURE;
	}

	/*check sequential and other settings (probalby will need to change*/
	/*some of this checks                                              */
	if((thee->type == CCT_AMGX) || (thee->type == CCT_DUMMY)){
		if((!thee->setgrid) && (!thee->setglen)){
			Vnm_print(2, "CUparm_check: Neither GRID nor GLEN set!\n");
			rc = VRC_FAILURE;
		}
		if((thee->setgrid) && (thee->setglen)){
			Vnm_print(2, "CUparm_check: Both GRID and GLEN set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setgcent){
			Vnm_print(2, "CUparm_check: GCENT not set!\n");
			rc = VRC_FAILURE;
		}
	}

	/*check sequential an parallel auto settings*/
	if((thee->type == CCT_AUTO) || (thee->type == CCT_OTHER)){
		if(!thee->setcglen){
			Vnm_print(2, "CUparm_check: CGLEN not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setfglen){
			Vnm_print(2, "CUparm_check: FGLEN not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setcgcent){
			Vnm_print(2, "CUparm_check: CGCENT not set!\n");
			rc = VRC_FAILURE;
		}
	}

	/*check parallel and focusing settings*/
	if(thee->type == CCT_PARALLEL){
		if(!thee->setpdime){
			Vnm_print(2, "CUparm_check: PDIME not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setrank){
			Vnm_print(2, "CUparm_check: PROC_RANK not set!\n");
			rc = VRC_FAILURE;
		}
		if(!thee->setofrac){
			Vnm_print(2, "CUparm_check: OFRAC not set!\n");
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

VPUBLIC Vrc_Codes CUparm_parseToken(CUparm *thee, char tok[VMAX_BUFSIZE],
		Vio *sock){

		if(thee == VNULL){
			Vnm_print(2, "parseCU: got NULL thee!\n");
			return VRC_WARNING;
		}
		if(sock == VNULL){
			Vnm_print(2, "parseCU: got Null socket!\n");
			return VRC_WARNING;
		}

		Vnm_print(0, "CUparm_parseToken: trying %s...\n", tok);

		if(Vstring_strcasecmp(tok, "dime") == 0){
			return CUparm_parseDIME(thee, sock);
		} else if (Vstring_strcasecmp(tok, "chgm") == 0){
			return CUparm_parseCHGM(thee, sock);
		} else if(Vstring_strcasecmp(tok, "etol") == 0){
			return CUparm_parseETOL(thee, sock);
		} else if(Vstring_strcasecmp(tok, "grid") == 0){
			return CUparm_parseGRID(thee, sock);
		} else if(Vstring_strcasecmp(tok, "glen") == 0){
			return CUparm_parseGLEN(thee, sock);
		} else if(Vstring_strcasecmp(tok, "gcent") == 0){
			return CUparm_parseGCENT(thee, sock);
		} else if(Vstring_strcasecmp(tok, "cglen") == 0){
			return CUparm_parseCGLEN(thee, sock);
		} else if(Vstring_strcasecmp(tok, "fglen") == 0){
			return CUparm_parseFGLEN(thee, sock);
		} else if(Vstring_strcasecmp(tok, "cgcent") == 0){
			return CUparm_parseCGCENT(thee, sock);
		} else if(Vstring_strcasecmp(tok, "fgcent") == 0){
			return CUparm_parseFGCENT(thee, sock);
		} else if(Vstring_strcasecmp(tok, "pdime") == 0){
			return CUparm_parsePDIME(thee, sock);
		} else if(Vstring_strcasecmp(tok, "ofrac") == 0){
			return CUparm_parseOFRAC(thee, sock);
		} else if(Vstring_strcasecmp(tok, "async") == 0){
			return CUparm_parseASYNC(thee, sock);
		} else if(Vstring_strcasecmp(tok, "gamma") == 0){
			return CUparm_parseGAMMA(thee, sock);
		} else if(Vstring_strcasecmp(tok, "useaqua") == 0){
			return CUparm_parseUSEAQUA(thee, sock);
		} else {
			Vnm_print(2, "parseCU: Unrecogized keyword (%s)!\n", tok);
			return VRC_WARNING;
		}

		return VRC_FAILURE;
}

VPRIVATE Vrc_Codes CUparm_parseDIME(CUparm *thee, Vio *sock){

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

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;

}

VPRIVATE Vrc_Codes CUparm_parseCHGM(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	Vchrg_Meth ti;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", (int*)(&ti)) == 1){
		thee->chgm = ti;
		thee->setchgm = 1;
		Vnm_print(2, "NOsh: Warning -- parsed deprecated statement \"chgm %d\,\n", ti);
		Vnm_print(2, "Noah please use \"chgm ");
		switch(thee->chgm){
		case VCM_TRIL:
			Vnm_print(2, "spl0");
			break;
		case VCM_BSPL2:
			Vnm_print(2, "spl2");
			break;
		case VCM_BSPL4:
			Vnm_print(2, "spl4");
			break;
		default:
			Vnm_print(2, "UNKNOWN");
			break;
		}
		Vnm_print(2,"\" instead!\n");
		return VRC_SUCCESS;
	} else if(Vstring_strcasecmp(tok, "spl0") == 0){
		thee->chgm = VCM_TRIL;
		thee->setchgm = 1;
		return VRC_SUCCESS;
	} else if(Vstring_strcasecmp(tok, "spl2") == 0){
		thee->chgm = VCM_BSPL2;
		thee->setchgm = 1;
		return VRC_SUCCESS;
	} else if(Vstring_stscasecmp(tok, "spl4") == 0){
		thee->chgm = VCM_BSPL4;
		thee->setchgm = 1;
		return VRC_SUCCESS;
	} else {
		Vnm_print(2, "NOsh: Unrecognized parameter (%s) when parsing chgm!n", tok);
		return VRC_WARNING;
	}

	return VRC_WARNING;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseETOL(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing etol keyword!\n", tok);
		return VRC_WARNING;
	} else if(tf <= 0.0){
		Vnm_print(2, "parseCU: etol must be greater than 0!\n");
	} else {
		thee->etol = tf;
	}

	thee->setetol = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseGRID(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing GRID keyword!\n", tok);
		return VRC_WARNING;
	} else {
		thee->grid[0] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing GRID keyword!\n", tok);
		return VRC_WARNING;
	} else {
		thee->grid[1] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%ld", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing GRID keyword!\n", tok);
		return VRC_WARNING;
	} else {
		thee->grid[2] = tf;
	}

	thee->setgrid = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseGLEN(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing GLEN keyword!\n", tok);
		return VRC_WARNING;
	} else {
		thee->glen[0] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing GLEN keyword!\n", tok);
		return VRC_WARNING;
	} else {
		thee->glen[1] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "Nosh: Read non-float (%d) while parsing GLEN keyword!\n", tok);
		return VRC_WARNING;
	} else {
		thee->glen[2] = tf;
	}

	thee->setglen = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseGAMMA(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];

	VJMPERR1(Vio_scabf(sock, "%s", tok) == 1);
	Vnm_print(2, "parseCU: GAMMA keyword deprecated!\n");
	Vnm_print(2, "parseCU: If you are using PyMol or VMD and still seeing this message,\n");
	Vnm_print(2, "parseCU: please contact the developers of those programs regarding this message.\n");

	return VRC_SUCCESS;
	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseGCENT(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;
	int ti;

	/*if the next token is not a float, it probable means we want to
	 * center on a molecule.
	 */

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		if(Vstring_strcasecmp(tok, "mol") == 0){
			VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
			if(sscanf(tok, "%d", &ti) == 0){
				Vnm_print(2, "NOsh: Read not-int (%s) while parsing GCENT MOL keyword!\n", tok);
				return VRC_WARNING;
			}
			else {
				thee->cmeth = CCT_MOLECULE;
				/*subtract 1 here to convert user numbering into array index*/
				thee->centmol = ti - 1;
			}
		}
		else {
			Vnm_print(2, "NOsh: Read non-int (%s) while parsing GCENT MOL keyword!\n", tok);
			return VRC_WARNING;
		}
	}
	else {
		thee->center[0] = tf;

		VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
		if(sscanf(tok, "%lf", &tf) == 0){
			Vnm_print(2, "NOsh: Read non-float (%s) while parsing GCENT keyword!\n", tok);
			return VRC_WARNING;
		}
		thee->center[1] = tf;

		VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
		if(sscanf(tok, "%lf", &tf) == 1){
			Vnm_print(2, "NOsh: Read non-float (%s) while parsing GCENT keyword!\n", tok);
			return VRC_WARNING;
		}
		thee->center[2] = tf;
	}

	thee->setgcent = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseCGLEN(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing CGLEN keyword!\n", tok);
		return VRC_WARNING;
	}
	else {
		thee->cglen[0] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_printf(2, "NOsh: Read non-float (%s) while parsing CGLEN keyword!\n", tok);
		return VRC_WARNING;
	}
	else {
		thee->cglen[1] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing CGLEN keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->cglen[2] = tf;
	}

	thee->setcglen = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseFGLEN(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing FGLEN keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->fglen[0] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-float (%s) while parsing FGLEN keyword!\n", tok);
		return VRC_WARNING;
	}
	else {
		thee->fglen[1] = tf;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-flaot (%s) while parsing FGLEN keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->fglen[2] = tf;
	}

	thee->setfglen = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseCGCENT(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;
	int ti;

	/*if the next token is not a float, it probably means we want to
	 * center on a molecule
	 */
	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		if(Vstring_strcasecmp(tok, "mol") == 0){
			VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
			if(sscanf(tok, "%d", &ti) == 0){
				Vnm_print(2, "NOsh: Read non-int (%s) while parsing CFCENT MOL keyword!\n", tok);
				return VRC_WARNING;
			}
			else{
				thee->ccmeth = CCT_MOLECULE;
				/*subtract 1 here to convert user numbering into array index*/
				thee->ccentmol = ti -1;
			}
		}
		else {
			Vnm_print(2, "NOsh: Unexpected keyword (%s) while parsing CGCENT!\n", tok);
			return VRC_WARNING;
		}
	}
	else{
		thee->ccenter[0] = tf;

		VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
		if(sscanf(tok, "%lf", &tf) == 0){
			Vnm_print(2, "NOsh: Read non-float (%s) while parsing CGCENT keyword!\n", tok);
			return VRC_WARNING;
		}
		else{
			thee->ccenter[1] = tf;
		}

		VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
		if(sscanf(tok, "%lf", &tf) == 0){
			Vnm_print(2, "NOsh: Read non-float (%s) while parsing CGCENT keyword!\n", tok);
			return VRC_WARNING;
		}
		else{
			thee->ccenter[2] = tf;
		}

	}

	thee->setcgcent = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseFGCENT(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;
	int ti;

	/*if the next token is not a float, it probably means we watn to
	 * center on a molecule.
	 */
	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		if(Vstring_strcasecmp(tok, "mol") == 0){
			VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
			if(sscanf(tok, "%d", &ti) == 0){
				Vnm_print(2, "NOsh: Read non-int (%s) while parsing FGCENT MOL keyword!\n");
				return VRC_WARNING;
			}
			else{
				thee->fcmeth = CCT_MOLECULE;
				/*subtract 1 here to convert user numbering to array index*/
				thee->fcentmol = ti -1;
			}
		}
		else{
			Vnm_print(2, "NOsh: Unexpected keyword (%s) while parsing FGCENT!\n", tok);
			return VRC_WARNING;
		}
	}
	else{
		thee->fcenter[0] = tf;

		VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
		if(sscanf(tok, "%lf", &tf) == 0){
			Vnm_print(2, "NOsh: Read non-float (%s) while parsing FGCENT keyword!\n", tok);
			return VRC_WARNING;
		}
		thee->fcenter[1] = tf;

		VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
		if(sscanf(tok, "%lf", &tf) == 0){
			Vnm_print(2, "NOsh: Read non-float (%s) while parsing FGCENT keyword!\n", tok);
			return VRC_WARNING;
		}
		thee->fcenter[2] = tf;
	}

	thee->setfgcent = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parsePDIME(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	int ti;

	/*read the number of grid points*/
	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", &ti) == 0){
		Vnm_print(2, "NOsh: Read non-integer (%s) while parsing PDIME keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->pdime[0] = ti;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", &ti) == 0){
		Vnm_print(2, "NOsh: Read non-integer (%s) while parsing PDIME keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->pdime[1] = ti;
	}

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%d", &ti) == 0){
		Vnm_print(2, "NOsh: Read non-integer (%s) while parsing PDIME keyword!\n", tok);
		return VRC_WARNING;
	}
	else{
		thee->pdime[2] = ti;
	}

	thee->setpdime = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseOFRAC(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	double tf;

	VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%lf", &tf) == 0){
		Vnm_print(2, "NOsh: Read non-int (%s) while parsing OFRAC keyword!\n", tok);
		return VRC_WARNING;
	}

	thee->ofrac = tf;
	thee->setofrac = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseASYNC(CUparm *thee, Vio *sock){

	char tok[VMAX_BUFSIZE];
	int ti;

	VJEMPERR1(Vio_scanf(sock, "%s", tok) == 1);
	if(sscanf(tok, "%i", &ti) == 0){
		Vnm_print(2, "NOsh: Read non-integer (%s) while parsing ASYNC keyword!\n", tok);
		return VRC_SUCCESS;
	}

	thee->async = ti;
	thee->setasync = 1;
	return VRC_SUCCESS;

	VERROR1:
		Vnm_print(2, "parseCU: ran out of tokens!\n");
		return VRC_WARNING;
}

VPRIVATE Vrc_Codes CUparm_parseUSEAQUA(CUparm *thee, Vio *sock){

	Vnm_print(0, "NOsh: parsed useaqua\n");
	thee->useAqua =1;
	thee->setUseAqua = 1;
	return VRC_SUCCESS;
}
