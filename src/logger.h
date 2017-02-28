/**
 *    @file  logger.h
 *   @brief  Logger header file.
 *
 *  This header file wraps the external logging library and sets up the general
 *  configuration of the library required for the entire package.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Feb 26, 2017  12:59
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef LOGGER_H__
#define LOGGER_H__

#include "base.h"

// Setup Easylogging++ performance tracker to report in microseconds.
#define ELPP_PERFORMANCE_MICROSECONDS
// Prevent creation of default empty log file during pre-processing.
#define ELPP_NO_DEFAULT_LOG_FILE

#ifdef NDEBUG
// Disable Easylogging++ performance tracking logs.
#define ELPP_DISABLE_PERFORMANCE_TRACKING
#endif  // NDEBUG

#include <easyloggingpp/src/easylogging++.h>

#endif  // LOGGER_H__
