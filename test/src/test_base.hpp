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

#ifndef TEST_BASE_HPP__
#define TEST_BASE_HPP__

#include <string>

#include "catch2/catch.hpp"

#ifndef TESTDIR
#define TESTDIR ".."  // Test directory path relative to test binary path.
#endif

std::string _testdir(TESTDIR);

#endif  // TEST_BASE_HPP__
