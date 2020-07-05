/**
 *    @file  test_base.hpp
 *   @brief  Test base header file.
 *
 *  This header file includes essential macros such as test data directory which can be
 *  redefine during compile time and definitions here are just default values.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Feb 26, 2017  17:36
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_TEST_BASE_HPP__
#define PSI_TEST_BASE_HPP__

#include <string>

#include "catch2/catch.hpp"

#include "test_config.hpp"

#define TEMPLATE_SCENARIO TEMPLATE_TEST_CASE

#ifndef TEST_DATA_DIR
#define TEST_DATA_DIR PROJECT_SOURCE_DIR "/test/data"
#endif

static const std::string test_data_dir( TEST_DATA_DIR );

#endif  /* --- #ifndef PSI_TEST_BASE_HPP__ --- */
