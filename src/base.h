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

#ifndef BASE_H__
#define BASE_H__

// Define/Undefine NDEBUG according to GREM_DEBUG macro value specified at compile time.
#undef NDEBUG
#if !GREM_DEBUG
#define NDEBUG
#endif

// TODO: Move to a class specifically aimed to store performance statistics.
// Add a checkpoint for performance tracker when this number of loci is traversed.
#define TRAVERSE_CHECKPOINT_LOCI_NO 1000000

#endif  // BASE_H__
