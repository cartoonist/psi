/*
 * =====================================================================================
 *
 * Filename: base.h
 *
 * Created: Tue Dec 06, 2016  22:37
 * Last modified: Mon Dec 19, 2016  13:06
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

// Define/Undefine NDEBUG according to GREM_DEBUG value.
#undef NDEBUG
#if !GREM_DEBUG
#define NDEBUG
// Disable Easylogging++ performance tracking logs.
#define ELPP_DISABLE_PERFORMANCE_TRACKING
#endif

// Add a checkpoint for performance tracker when this number of loci is traversed.
#define TRAVERSE_CHECKPOINT_LOCI_NO 100000
// Report the number of seed hits when it exceeds this number.
#define SEEDHITS_REPORT_BUF 100000
// Report the average number of go downs when it exceeds this number of samples.
#define AVG_GODOWNS_SAMPLES 10000

#endif  // BASE_H__
