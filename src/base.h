/**
 *    @file  base.h
 *   @brief  Base header file.
 *
 *  The base header file including essential macros that should be defined globally
 *  in the entire package. They can be redefined at compile time.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Dec 06, 2016  22:37
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_BASE_H__
#define PSI_BASE_H__

// Define/Undefine NDEBUG according to GREM_DEBUG macro value specified at compile time.
#undef NDEBUG
#if !GREM_DEBUG
#define NDEBUG
#else
#define GREM_DEBUG_ENABLED
#endif

#undef GIT_VERSION
#ifdef GREM_GIT_VERSION
#define GIT_VERSION GREM_GIT_VERSION
#else
#define GIT_VERSION VERSION
#endif

#undef UPDATE_DATE
#ifdef GREM_UPDATE_DATE
#define UPDATE_DATE GREM_UPDATE_DATE
#else
#define UPDATE_DATE __DATE__
#endif

#endif  /* --- #ifndef PSI_BASE_H__ --- */
