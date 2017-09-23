/**
 *    @file  tests_stat.cc
 *   @brief  Test stat module.
 *
 *  This test contains test scenarios for stat module.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Sep 24, 2017  02:48
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <iostream>
#include <chrono>
#include <thread>

#include "tests_base.h"
#include "stat.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;
using namespace std::chrono_literals;


SCENARIO ( "Test the Timer", "[stat]" )
{
  GIVEN( "A Timer" )
  {
    WHEN( "It measures a time period" )
    {
      Timer::clock_type::time_point start, end;
      {
        auto timer = Timer( "test-timer" );
        start = Timer::clock_type::now();
        std::this_thread::sleep_for( 678912us );
      }
      end = Timer::clock_type::now();
      THEN( "It should get the correct duration in microseconds" )
      {
        auto d1 = std::chrono::duration_cast< std::chrono::microseconds >( end - start );
        auto d2 = Timer::get_duration( "test-timer" );
        REQUIRE( d2.count() - d1.count() < 100 );
      }
    }
  }
}
