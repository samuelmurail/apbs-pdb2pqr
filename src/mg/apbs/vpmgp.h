/** @defgroup Vpmgp Vpmgp class
 *  @brief  Parameter structure for Mike Holst's PMGP code
 *  @note   Variables and many default values taken directly from PMG
 */

/**
 *  @file     vpmgp.h
 *  @ingroup  Vpmgp
 *  @brief    Contains declarations for class Vpmgp
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *  @note     Variables and many default values taken directly from PMG
  *  @attention
 *  @verbatim
 *

 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002.  Washington University in St. Louis.
 * All Rights Reserved.
 *
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
 * AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS".  THE AUTHORS HAVE NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

 * @endverbatim
 */


#ifndef _VPMGP_H_
#define _VPMGP_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/valist.h"
#include "apbs/vunit.h"

/**
 *  @struct  Vpmgp
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpmgp class/module
 */
struct Vpmgp {

    /* ********** USER-SPECIFIED PARAMETERS ********** */
    int nx;              /**< Grid x dimensions [no default]  */
    int ny;              /**< Grid y dimensions [no default]  */
    int nz;              /**< Grid z dimensions [no default]  */
    int nlev;            /**< Number of mesh levels [no default] */
    double hx;           /**< Grid x spacings [no default]  */
    double hy;           /**< Grid y spacings [no default]  */
    double hzed;         /**< Grid z spacings [no default]  */
    int nonlin;          /**< Problem type [no default]
                              \li 0: linear
                              \li 1: nonlinear
                              \li 2: linear then nonlinear */

    /* ********** DERIVED PARAMETERS ********** */
    int nxc;             /**< Coarse level grid x dimensions */
    int nyc;             /**< Coarse level grid y dimensions */
    int nzc;             /**< Coarse level grid z dimensions */
    int nf;              /**< Number of fine grid unknowns */
    int nc;              /**< Number of coarse grid unknowns */
    int narrc;           /**< Size of vector on coarse level */
    int n_rpc;           /**< Real info work array required storage */
    int n_iz;            /**< Integer storage parameter (index max) */
    int n_ipc;           /**< Integer info work array required storage */

    int nrwk;            /**< Real work storage */
    int niwk;            /**< Integer work storage */
    int narr;            /**< Array work storage */
    int ipkey;           /**< Toggles nonlinearity (set by nonlin)
                              \li  -1: Linearized PBE
                              \li   0: Nonlinear PBE with capped sinh 
                                     term [default]
                              \li  >1: Polynomial approximation to sinh, 
                                     note that ipkey must be odd  */

    /* ********** PARAMETERS WITH DEFAULT VALUES ********** */
    double xcent;        /**< Grid x center [0]  */
    double ycent;        /**< Grid y center [0]  */
    double zcent;        /**< Grid z center [0]  */
    double errtol;       /**< Desired error tolerance [default = 1e-9] */
    int itmax;           /**< Maximum number of iters [default = 100] */
    int istop;           /**< Stopping criterion [default = 1]
                              \li 0: residual
                              \li 1: relative residual
                              \li 2: diff
                              \li 3: errc
                              \li 4: errd
                              \li 5: aerrd */
    int iinfo;           /**< Runtime status messages [default = 1]
                              \li 0: none
                              \li 1: some
                              \li 2: lots
                              \li 3: more */
    int bcfl;            /**< Boundary condition method [default = 1]
                              \li   0: zero boundary conditions
                              \li   1: boundary condition approximated by
                                      single Debye-Huckel sphere for
                                      entire molecule
                              \li   2: boundary condition approximated by 
                                      Debye-Huckel spheres for each atom 
                              \li   4: boundary condition determined from 
                                      previously calculated solution 
                                      (i.e., for focusing calculations) */
    int key;             /**< Print solution to file [default = 0] 
                              \li   0: no
                              \li   1: yes */
    int iperf;           /**< Analysis of the operator [default = 0]
                              \li   0: no
                              \li   1: condition number
                              \li   2: spectral radius
                              \li   3: cond. number & spectral radius */
    int meth;            /**< Solution method [default = 2]
                              \li   0: conjugate gradient multigrid
                              \li   1: newton
                              \li   2: multigrid
                              \li   3: conjugate gradient
                              \li   4: sucessive overrelaxation
                              \li   5: red-black gauss-seidel
                              \li   6: weighted jacobi
                              \li   7: richardson */
    int mgkey;           /**< Multigrid method [default = 0]
                              \li   0: variable v-cycle
                              \li   1: nested iteration */
    int nu1;             /**< Number of pre-smoothings [default = 2] */
    int nu2;             /**< Number of post-smoothings [default = 2] */
    int mgsmoo;          /**< Smoothing method [default = 1]
                              \li   0: weighted jacobi
                              \li   1: gauss-seidel
                              \li   2: SOR
                              \li   3: richardson
                              \li   4: cghs */
    int mgprol;          /**< Prolongation method [default = 0]
                              \li   0: trilinear
                              \li   1: operator-based
                              \li   2: mod. operator-based */
    int mgcoar;          /**< Coarsening method [default = 2]
                              \li   0: standard
                              \li   1: harmonic
                              \li   2: galerkin */
    int mgsolv;          /**< Coarse equation solve method [default = 1]
                              \li   0: cghs
                              \li   1: banded linpack */
    int mgdisc;          /**< Discretization method [default = 0]
                              \li   0: finite volume
                              \li   1: finite element */
    double omegal;       /**< Linear relax parameter [default = 8e-1] */
    double omegan;       /**< Nonlin relax parameter [default = 9e-1] */
    int irite;           /**< FORTRAN output unit [default = 8] */
    int ipcon;           /**< Preconditioning method [default = 3]
                              \li   0: diagonal
                              \li   1: ICCG 
                              \li   2: ICCGDW
                              \li   3: MICCGDW
                              \li   4: none */
    double xlen;        /**< Domain x length */
    double ylen;        /**< Domain y length */
    double zlen;        /**< Domain z length */
    double xmin;        /**< Domain lower x corner */
    double ymin;        /**< Domain lower y corner */
    double zmin;        /**< Domain lower z corner */
    double xmax;        /**< Domain upper x corner */
    double ymax;        /**< Domain upper y corner */
    double zmax;        /**< Domain upper z corner */
};

/** @typedef Vpmgp
 *  @ingroup Vpmgp
 *  @brief   Declaration of the Vpmgp class as the Vpmgp structure
 */
typedef struct Vpmgp Vpmgp;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Inlineable methods (vpmgp.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPMGP)
#else /* if defined(VINLINE_VPMGP) */
#endif /* if !defined(VINLINE_VPMGP) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Non-Inlineable methods (vpmgp.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct PMG parameter object and initialize to default values
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   nx    Number of x grid points
 *  @param   ny    Number of y grid points
 *  @param   nz    Number of z grid points
 *  @param   nlev  Number of levels in multigrid hierarchy
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   nonlin  Nonlinearity flag
 *                   \li 0: Linearized PBE
 *                   \li 1: Nonlinear PBE
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC Vpmgp* Vpmgp_ctor(int nx, int ny, int nz, int nlev, 
  double hx, double hy, double hzed, int nonlin);

/** @brief   FORTRAN stub to construct PMG parameter object and initialize to
 *           default values 
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   thee  Newly allocated PMG object
 *  @param   nx    Number of x grid points
 *  @param   ny    Number of y grid points
 *  @param   nz    Number of z grid points
 *  @param   nlev  Number of levels in multigrid hierarchy
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   nonlin  Nonlinearity flag
 *                   \li 0: Linearized PBE
 *                   \li 1: Nonlinear PBE
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vpmgp_ctor2(Vpmgp *thee, int nx, int ny, int nz, int nlev, 
  double hx, double hy, double hzed, int nonlin);

/** @brief   Object destructor
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location for Vpmgp object
 */
VEXTERNC void Vpmgp_dtor(Vpmgp **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   thee  Pointer to Vpmgp object
 */
VEXTERNC void Vpmgp_dtor2(Vpmgp *thee);

#endif    /* ifndef _VPMGP_H_ */
