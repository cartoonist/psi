/*
 * =====================================================================================
 *
 * Filename: base.h
 *
 * Created: Tue Dec 06, 2016  22:37
 * Last modified: Tue Feb 28, 2017  17:43
 *
 * Description: grem base definitions.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
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
// Report the average path length when no. of traversed paths reach this number.
#define NOF_PATHLEN_SAMPLES 1000000

#endif  // BASE_H__
