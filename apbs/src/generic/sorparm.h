/** @defgroup SORparm SORparm class
 *  @brief    Parameter which holds useful parameters for generic apbs
 *            calculations
 */

/**
 *  @file     sorparm.h
 *  @ingroup  SORparm
 *  @brief    Contains declarations for class SORparm
 *  @version  $Id$
 *  @author   Juan M. Brandi
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

/* Generic header files */
#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"
#include "generic/mgparm.h"

/**
 * @brief  Calculation type
 * @ingroup SORparm
 */
enum eSORparm_CalcType {
    SCT_AUTO=0,  /**< SOR-auto */
    SCT_PARALLEL=1,  /**< sor-para */
    SCT_DUMMY=2,  /**< sor-dummy */
    SCT_NONE=3  /**< unspecified */
};

/**
 * @brief  Declare MGparm_CalcType type
 * @ingroup  MGparm
 */
typedef enum eSORparm_CalcType SORparm_CalcType;

/**
 * @brief  Centering method
 * @ingroup SORparm
 */
enum eSORparm_CentMeth {
    SCM_POINT=0, /**< Center on a point */
    SCM_MOLECULE=1,  /**< Center on a molecule */
};

/**
 * @brief  Declare SORparm_CentMeth type
 * @ingroup  SORparm
 */
typedef enum eSORparm_CentMeth SORparm_CentMeth;

/**
 *  @ingroup SORparm
 *  @author  Juan Brandi
 *  @brief   Parameter structure for SOR-specific variables from input files
 */
struct sSORparm {

	SORparm_CalcType type;  /**< What type of SOR run? */
    int parsed;  /**< Has this structure been filled? (0 = no, 1 = yes) */

    /* *** GENERIC PARAMETERS *** */
    int dime[3];  /**< Grid dimensions */
    int setdime;  /**< Flag, @see dime */
    Vchrg_Meth chgm;  /**< Charge discretization method */
    int setchgm;  /**< Flag, @see chgm */
    Vchrg_Src  chgs; /**< Charge source (Charge, Multipole, Induced Dipole,
                      * NL Induced */

    /* *** TYPE 0 PARAMETERS (SEQUENTIAL MANUAL) *** */
    double etol;  /**< User-defined error tolerance */
    int setetol;  /**< Flag, @see etol */
    double grid[3];  /**< Grid spacings */
    int setgrid;  /**< Flag, @see grid */
    double glen[3];  /**< Grid side lengths. */
    int setglen;  /**< Flag, @see glen */
    SORparm_CentMeth cmeth;  /**< Centering method */
    double center[3];  /**< Grid center. If ispart = 0, then this is
                        * only meaningful if cmeth = 0.  However, if
                        * ispart = 1 and cmeth = MCM_PNT, then this is the
                        * center of the non-disjoint (overlapping)
                        * partition.  If ispart = 1 and cmeth = MCM_MOL, then
                        * this is the vector that must be added to the
                        * center of the molecule to give the center of
                        * the non-disjoint partition.  */
    int centmol;  /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setgcent;  /**< Flag, @see cmeth */

    /* ******** TYPE 1 & 2 PARAMETERS (SEQUENTIAL & PARALLEL AUTO-FOCUS) *** */
    double fglen[3];  /**< Fine grid side lengths */
    int setfglen;  /**< Flag, @see fglen */
    SORparm_CentMeth fcmeth;
    int fcentmol;
    double fcenter[3];

    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    /**TODO: this sections is probably not going to work as is,
     * and will change later.
     */

    int nonlintype; /**< Linearity Type Method to be used */
    int setnonlintype; /**< Flag, @see nonlintype */

    int useAqua;  /**< Enable use of lpbe/aqua */
    int setUseAqua; /**< Flag, @see useAqua */
};


/** @typedef SORparm
 *  @ingroup SORparm
 *  @brief   Declaration of the SORparm class as the SORparm structure
 */
typedef struct sSORparm SORparm;

VPUBLIC SORparm* SORparm_ctor(SORparm_CalcType type);

VPUBLIC Vrc_Codes SORparm_ctor2(SORparm *thee, SORparm_CalcType type);

VPUBLIC void SORparm_dtor(SORparm **thee);

VPUBLIC void SORparm_dtor2(SORparm *thee);

VPUBLIC int SORparm_copyMGparm(MGparm *mgparm, SORparm *sorparm);

#endif


