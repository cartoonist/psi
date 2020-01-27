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


using namespace grem;
using namespace std::chrono_literals;

SCENARIO ( "Test the Timer", "[stat]" )
{
  GIVEN( "A CPU clock Timer" )
  {
    WHEN( "It measures a short time period" )
    {
      {
        auto timer = Timer< >( "test-timer" );
        std::this_thread::sleep_for( 678912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< >::get_duration_rep( "test-timer" );
        REQUIRE( d == Approx( 0 ).margin( 0.0001 ) );
      }
    }

    WHEN( "It measures longer time period" )
    {
      {
        auto timer = Timer< >( "test-timer" );
        std::this_thread::sleep_for( 1278912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< >::get_duration_rep( "test-timer" );
        REQUIRE( d == Approx( 0 ).margin( 0.0001 ) );
      }
    }
  }

  GIVEN( "A wall clock Timer" )
  {
    WHEN( "It measures a short time period" )
    {
      {
        auto timer = Timer< SteadyClock >( "test-timer" );
        std::this_thread::sleep_for( 678912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< SteadyClock >::get_duration_rep( "test-timer" );
        REQUIRE( static_cast< float >( d ) == Approx( 678912 ).epsilon( 0.001 ) );
      }
    }

    WHEN( "It measures longer time period" )
    {
      {
        auto timer = Timer< SteadyClock >( "test-timer" );
        std::this_thread::sleep_for( 1278912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< SteadyClock >::get_duration_rep( "test-timer" );
        REQUIRE( static_cast< float >( d ) == Approx( 1278912 ).epsilon( 0.001 ) );
      }
    }
  }
}
