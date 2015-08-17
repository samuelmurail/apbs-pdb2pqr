/** @defgroup SORparm SORparm class
 *  @brief    Parameter which holds useful parameters for generic SOR
 *            calculations
 */

/**
 *  @file     sorparm.h
 *  @ingroup  SORparm
 *  @brief    Contains declarations for class SORparm
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

#ifndef _SORPARM_H_
#define _SORPARM_H_

/*Generic header file*/
#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"

/**
 * @ingroup SORparm
 * @brief Centering Method
 */
enum eSORparm_CentMeth{
	sMCM_POINT=0, /**< Center on a point */
	sMCM_MOLECULE=1,  /**< Center on a molecule */
	sMCM_FOCUS=2  /**< Determined by focusing */
};

/**
 * @ingroup SORparm
 * @brief Declare SORparm_CentMeth type.
 */
typedef enum eSORparm_CentMeth SORparm_CentMeth;

/**
 * @ingroup SORparm
 * @brief Parameter structure for SOR specific variables from input files.
 * @author Juan Brandi
 */
struct sSORparm{

	int parsed; /**< Has the structure benn filled? (0=no, 1=yes) */

	/* ***GENERIC PAREMETERS*** */
	int dim[3]; /**< Grid Dimensions */
	int setdime;
	Vchrg_Meth chgm; /**< Charge discretization method */
	int setchgm;
	Vchrg_Src chgs; /**< Charge source (Charge, Multipole, Induce Dipole, NL Induced */

	/* *** TYPE 0 PARAMETERS *** */
	double etol; /**< User defined error tolerance */
	int setetol;
	double grid[3]; /**< Grid spacings */
	int setgrid;
	double glen[3]; /**< Grid side lengths */
	int setglen;
	SORparm_CentMeth cmeth; /**< Centering Method */
	double center[3]; /**< Grid Center */
	int centmol; /**< Particular molecule on which we want to center the grid.
					*This should be the appropriate index in an array of molecules,
					*not the positive definite integer specified by the user. */
	int setgcent;

	/* ***TYPE 1 PARAMETERS*** */
	double cglen[3]; /**< Coarse grid side lengths */
	int setcglen;
	double fglen[3]; /**< Fine grid side lengths */
	int setfglen[3];

	int method; /** Solver method to be used */
	int setmethod;

	int useAqua; /** Enable the use of lpbe/aqua */
	int setUseAqua;

};

/**
 * @ingroup SORparm
 * @typedef SORparm
 * @brief Declaration of the SORparm class as the SORparm structure.
 */
typedef struct sSORparm SORparm;

/**
 * @ingroup SORparm
 * @brief SOR object constructor
 * @param void
 * @returns Newly created and memory allocated SOR object
 */
VEXTERNC SORparm* SORparm_ctor();

/**
 * @ingroup SORparm
 * @brief SOR object constructor
 * @param pointer to locations of SORparm object
 * @returns Success enumeration
 */
VEXTERNC Vrc_Codes SORparm_ctor2(SORparm *thee);

/**
 * @ingroup SORparm
 * @brief Object destructor
 * @param thee Pointer to memory location of SORparm object
 */
VEXTERNC void SORparm_dtor(SORparm **thee);

/**
 * @ingroup SORparm
 * @brief Object destructor
 * @param thee Pointer to SORparm object
 */
VEXTERNC void SORparm_dtor2(SORparm *thee);

#endif

















