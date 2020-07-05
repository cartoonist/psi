/**
 *    @file  base.hpp
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

#ifndef PSI_BASE_HPP__
#define PSI_BASE_HPP__

#include "config.hpp"

/* Define/Undefine NDEBUG according to PSI_DEBUG macro value specified at compile time.
 * This is necessary because SeqAn manipulates NDEBUG. */
#undef NDEBUG
#if !PSI_DEBUG
#define NDEBUG
#else
#define PSI_DEBUG_ENABLED
#endif

#endif  /* --- #ifndef PSI_BASE_HPP__ --- */
