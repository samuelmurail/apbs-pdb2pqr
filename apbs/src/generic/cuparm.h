/** @defgroup CUparm CUparm class
 *  @brief Parameter which holds useful parameters for generic
 *  1 grid level gpu calculations
 */

/**
 *  @file     cuparm.h
 *  @ingroup  CUparm
 *  @brief    Contains declarations for class CUparm
 *  @version  $Id$
 *  @author   Juan Brandi
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

#ifndef _CUJACPARM_H_
#define _CUJACPARM_H_

/*generic header files*/
#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"

/**
 * @brief Calculation type
 * @ingroup CUjacParm
 */
enum eCUparm_CalcType{
	CCT_DUMMY = 0, /**< cu-sor*/
	CCT_JAC = 1, /**< cu-jac*/
	CCT_OTHER = 2, /**< this could be use as a hook for user plug-ins later on*/
	CCT_AUTO = 3, /**< cu-auto*/
	CCT_NONE = 4, /**< unspecified*/
	CCT_AMGX = 5, /**< cuda AMGX library.*/
	CCT_PARALLEL = 6 /**< may impement paralles and multithreading one day*/
};

/**
 * @brief Declare CUparm_CalcType type
 * @ingroup CUparm
 */
typedef enum eCUparm_CalcType CUparm_CalcType;

/**
 * @brief Centering Method
 * @ingroup CUparm
 */
enum eCUparm_CentMeth{
	CCT_POINT = 0, /**< Center on a point*/
	CCT_MOLECULE = 1, /**Center on a molecule*/
	CCT_FOCUS = 2 /**Determined by focusing*/
};

/**
 * @brief Declare CUparm_CentMetht type
 * @ingroup CUparm
 */
typedef enum eCUparm_CentMeth CUparm_CentMeth;

/**
 * @ingroup CUparm
 * @author Juan M. Brandi
 * @brief Parameter structure for different gpu calculations
 */
struct sCUparm{

	CUparm_CalcType type; /**< What type of gpu solver*/
	int parsed; /**< Has this struct been filled (0=no, 1=yes)*/

	/*generic parameters*/
	int dime[3]; /**< Grid dimensions*/
	int setdime; /**< Flag, @see dime*/
	Vchrg_Meth chgm; /**< Charge discretization methos*/
	int setchgm; /**< Flag, @see chgm*/
	Vchrg_Src chgs; /**< Charge source (Charge, Multipole, Induced
					 * Dipole, NL Induced*/

	/* Type 0 Parameters (Sequential Manual) */
	double etol; /**< User-defined error tolereance*/
	int setetol; /**< Flag, @see etol*/
	double grid[3]; /**< Grid spacings*/
	int setgrid; /**< Flag, @see grid*/
	double glen[3]; /**< Grid side lengths*/
	int setglen; /**< Flag, @see glen*/
	CUparm_CentMeth cmeth; /**< Centering Method*/
	double center[3]; /**< Grid center. This probably will have to be
	 	 	 	 	   * modified since it has to do with parrallel
	 	 	 	 	   * stuff and modified accordingly.
	 	 	 	 	   **/
	int centmol; /**< Particular molecule on which we want to center the
	              * grid. This number should be the index in the array of
	              * molecules, not the positive integer specified by the
	              * user
	              * */
	int setgcent; /**< Flag, @see centmol*/

	/* Type 1 & 2 parameters (sequential auto-focus) */
	/* some of this parameters may not be used or will have to be */
	/*modified later to maybe allow for multicore and gpu */
    double cglen[3];  /**< Coarse grid side lengths */
    int setcglen;  /**< Flag, @see cglen */
    double fglen[3];  /**< Fine grid side lengths */
    int setfglen;  /**< Flag, @see fglen */
    CUparm_CentMeth ccmeth;  /**< Coarse grid centering method */
    double ccenter[3];  /**< Coarse grid center.  */
    int ccentmol;  /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive integer specified by the user. */
    int setcgcent;  /**< Flag, @see ccmeth */
    CUparm_CentMeth fcmeth;  /**< Fine grid centering method */
    double fcenter[3];  /**< Fine grid center.  */
    int fcentmol; /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive integer specified by the user. */
    int setfgcent;  /**< Flag, @see fcmeth */

    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    /* same as above some of this will not be used right now but maybe*/
    /*later.                                                          */
    double partDisjCenter[3];  /**< This gives the center
                                     of the disjoint partitions */
    double partDisjLength[3];  /**< This gives the lengths of the disjoint
                                * partitions */
    int partDisjOwnSide[6];  /**< Tells whether the boundary points are ours
                              * (1) or not (0) */

    int pdime[3];  /**< Grid of processors to be used in calculation */
    int setpdime;  /**< Flag, @see pdime */
    int proc_rank;  /**< Rank of this processor */
    int setrank;  /**< Flag, @see proc_rank */
    int proc_size;  /**< Total number of processors */
    int setsize;  /**< Flag, @see proc_size */
    double ofrac;  /**< Overlap fraction between procs */
    int setofrac;  /**< Flag, @see ofrac */
    int async; /**< Processor ID for asynchronous calculation */
    int setasync; /**< Flag, @see asynch */

    int nonlintype; /**< Linearity Type Method to be used */
    int setnonlintype; /**< Flag, @see nonlintype */

    int method;		/**< Solver Method */
    int setmethod; /**< Flag, @see method */

    int useAqua;  /**< Enable use of lpbe/aqua */
    int setUseAqua; /**< Flag, @see useAqua */

};

/** @typedef CUparm
 *  @ingroup CUparm
 *  @brief Declaration of the CUparm class as the CUparm structure.
 */
typedef struct sCUparm CUparm;

/**@brief Get the number of grid points in the x direction
 * @author Juan M. Brandi (code copied mainly from Nathan Baker)
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return Number of points in the x direction.
 */
VEXTERNC int CUparm_getNx(CUparm *thee);

/**@brief Get the number of grid points in the y direction
 * @author Juan M. Brandi (code copied mainly from Nathan Baker)
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return Number of points in the y direction.
 */
VEXTERNC int CUparm_getNy(CUparm *thee);

/**@brief Get the number of grid points in the z direction.
 * @author Juan M. Brandi (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return NUmber of points in the z direction.
 */
VEXTERNC int CUparm_getNz(CUparm *thee);

/**@brief Get grid spacing in the x direction
 * @author Juan M. Brandi (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object
 * @return Grid spacing in the x direction
 */
VEXTERNC double CUparm_getHx(CUparm *thee);

/**@brief Get grid spacing in the y direction.
 * @author Juan M. Brandi (code mainly copied from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object
 * @return Grid spacing in the y direction.
 */
VEXTERNC double CUparm_getHy(CUparm *thee);

/**@brief Get grid spacing in the z direction.
 * @author Juan M. Brandi (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return Grid spacing in the z direction
 */
VEXTERNC double CUparm_getHz(CUparm *thee);

/**@brief Set center x-coordinate.
 * @author Juan M. Brandi (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @param x x-coordinate
 */
VEXTERNC void CUparm_setCenterX(CUparm *thee, double x);

/**@brief Set center y-coordinate.
 * @author Juan M. Brandi (code copied mainly from Nathan. Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @param y y-coordinate.
 */
VEXTERNC void CUparm_setCenterY(CUparm *thee, double y);

/**@brief Set center z-coordinate.
 * @author Juan M. Brandi. (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm
 * @param z z-coordinate
 */
VEXTERNC void CUparm_setCenterZ(CUparm *thee, double z);

/**@brief Get center x-coordinate.
 * @author Juan M. Brandi. (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object
 * @return x-coordinate
 */
VEXTERNC double CUparm_getCenterX(CUparm *thee);

/**@brief Get center y-coordinate.
 * @author Juan M. Brandi. (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return y-coordinate
 */
VEXTERNC double CUparm_getCenterY(CUparm *thee);

/**@brief Get center z-coordinate.
 * @author Juan M. Brandi. (code copied mainly from Natha Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return z-coordinate.
 */
VEXTERNC double CUparm_getCenterZ(CUparm *thee);

/**@brief Construct CUparm object
 * @author Juan M.Brandi (code copied mainly from Nathan Baker)
 * @ingroup CUparm
 * @param type Type of CuSolver
 * @return Newly allocated and intialized CUparm object.
 */
VEXTERNC CUparm* CUparm_ctor(CUparm_CalcType type);

/**
 * @brief FORTRAN stub to construct CUparm object.
 * @author Juan M. Brandi (code copied from Nathan Baker and Todd Dolinsky.
 * @ingroup CUparm
 * @param thee Space for CUparm object.
 * @param type Type of CuSolver
 * @return Success enumeration.
 */
VEXTERNC Vrc_Codes CUparm_ctor2(CUparm *thee, CUparm_CalcType type);

/**@brief object destructor.
 * @author Juan M. Brandi (code copied from Nathan Baker).
 * @ingroup CUparm
 * @param the Pointer to Cujacparm object
 */
VEXTERNC void CUparm_dtor(CUparm **thee);

/**@brief FORTRAN stub for object destructor
 * @author Juan M. Brandi (code copied mainly from Nathan Baker)
 * @ingroup CUparm
 * @param thee CUparm object
 * @return Success enumeration
 */
VEXTERNC Vrc_Codes CUparm_ctor2(CUparm *thee);

/**@brief Consistency check for parameter values stored in object.
 * @author Juan M. Brandi (code copied mainly from Nathan Baker).
 * @ingroup CUparm
 * @param thee CUparm object.
 * @return Success enumemration.
 */
VEXTERNC Vrc_Codes CUparm_check(CUparm *thee);

/**@brief Parse a CU keyword from an input file.
 * @author Juan M. Brandi (code mostly copied from Nathan Baker).
 * @ingroup CUjackparm
 * @param thee CUparm object
 * @param tok Token to parse
 * @param sock Stream for more tokens
 * @return Success enumeration code (1 if matched and assigned; -1 if
 * 			matched, but there's some sort of error (i.e. too few args);
 * 			0 if not matched)
 */
VEXTERNC Vrc_Codes CUparm_parseToken(CUparm *thee, char tok[VMAX_BUFSIZE],
		Vio *sock);

#endif


