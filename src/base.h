/*
 * =====================================================================================
 *
 * Filename: base.h
 *
 * Created: Tue Dec 06, 2016  22:37
 * Last modified: Wed Dec 07, 2016  02:07
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

#endif  // BASE_H__
