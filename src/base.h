/*
 * =====================================================================================
 *
 * Filename: base.h
 *
 * Created: Tue Dec 06, 2016  22:37
 * Last modified: Fri Feb 03, 2017  00:05
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

// Easylogging++ performance tracker in microseconds.
#define ELPP_PERFORMANCE_MICROSECONDS
// Prevent creation of default empty log file during pre-processing.
#define ELPP_NO_DEFAULT_LOG_FILE

// Define/Undefine NDEBUG according to GREM_DEBUG value.
#undef NDEBUG
#if !GREM_DEBUG
#define NDEBUG
// Disable Easylogging++ performance tracking logs.
#define ELPP_DISABLE_PERFORMANCE_TRACKING
#endif

// Add a checkpoint for performance tracker when this number of loci is traversed.
#define TRAVERSE_CHECKPOINT_LOCI_NO 1000000
// Report the average path length when no. of traversed paths reach this number.
#define NOF_PATHLEN_SAMPLES 1000000

#endif  // BASE_H__
