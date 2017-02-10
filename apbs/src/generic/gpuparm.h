/**@defgroup GPUparm GPUparm class */

/**
 *  @file     gpuparm.h
 *  @ingroup  GPUparm
 *  @brief    Contains declarations for class GPUparm
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *
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


#ifndef _GPUPARM_H_
#define _GPUPARM_H_

/* Generic header files */
#include "maloc/maloc.h"
#include "generic/vhal.h"
#include "generic/vstring.h"
#include "generic/mgparm.h"

/**
 * @ingroup GPUparm
 * @author Juan Brandi (Mostly copied from Nathan Baker)
 * @brief Parameter structure for GPU-specific variables
 */
struct sGPUparm {

	MGparm mgparm; /**< Object to contain mg specific paramters*/


};

/**
 * @typedef GPUparm
 * @ingroup GPUparm
 * @brief Declaration of the GPUparm class as GPUparm structure
 */
typedef struct sGPUparm GPUparm;

/**
 * @brief 	Construct GPUparm object
 * @ingroup GPUparm
 * @author	Juan Brandi (mostly copied from Nathan Baker)
 * @param 	type Type of calculation
 * @returns	Newly allocated and initialized GPUparm object.
 */
VEXTERNC GPUparm* GPUparm_ctor(MGparm_CalcType type);

/**
 * @brief 	Object constructor
 * @ingroup	GPUparm
 * @author 	Juan Brandi (mostly copied from Nathan Baker)
 * @param 	thee Pointer to memory location of GPUparm object
 * @param	type Type of calculation.
 * @return	Success Enumeration.
 */
VEXTERNC Vrc_Codes GPUparm_ctor2(GPUparm *thee, MGparm_CalcType type);

/**
 * @brief	Object destructor
 * @ingroup	GPUparm
 * @author	Juan Brandi (mostly copied from Nathan Baker)
 * @param	thee Pointer to GPUparm object
 */
VEXTERNC void GPUparm_dtor(GPUparm **thee);

/**
 * @brief 	FORTRAN stuv for object destructor
 * @ingroup	GPUparm
 * @author	Juan Brandi (mostly copied from Nathan Baker)
 * @param	thee Pointer to GPUparm object
 */
VEXTERNC void GPUparm_dtor2(MGparm *thee);

/**
 * @brief 	Constistency check for parameters stored in object
 * @ingroup	GPUparm
 * @author	Juan Brandi (mostly copied from Nathan Baker)
 * @param	thee GPUparm object
 * @returns	Success enumeration
 */
VEXTERNC Vrc_Codes GPUparm_check(GPUparm *thee);

/**
 * @brief	Copy GPUparm object into a new one (thee)
 * @author	Juan Brandi (mostly copied from Nathan Baker)
 * @ingroup	GPUparm
 * @param	thee GPUparm object (target for copy);
 * @param	parm GPUparm object (source for copy);
 */
VEXTERNC void GPUparm_copy(GPUparm *thee, GPUparm *parm);

/**
 * @brief 	Parse GPU keyword from input file
 * @ingroup	SORparm
 * @author	Juan Brandi (moslty copied from Nathan Baker)
 * @param	thee GPUparm object
 * @param	tok Token to parse
 * @param	sock Stream for more tokens
 * @returns	Success enumeration (1 if matched and assigned; -1 if matched, but there's some sort of error, (i.e., too few args; 0 if not matched)
 */
VEXTERNC Vrc_Codes GPUparm_parseToken(GPUparm *thee, char tok[VMAX_BUFSIZE]);


#endif /* _GPUPARM_H_ */
