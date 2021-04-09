/**
 *    @file  test_stats.cpp
 *   @brief  Test stats module.
 *
 *  This test contains test scenarios for stats module.
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

#include <psi/stats.hpp>

#include "test_base.hpp"


using namespace psi;
using namespace std::chrono_literals;

SCENARIO ( "Test the Timer", "[stats]" )
{
  GIVEN( "A CPU clock Timer" )
  {
    WHEN( "It measures a short time period" )
    {
      {
        auto timer = Timer< CpuClock >( "timer-short" );
        std::this_thread::sleep_for( 678912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< CpuClock >::get_duration_rep( "timer-short" );
        REQUIRE( d == Approx( 0 ).margin( 0.0001 ) );
      }
    }

    WHEN( "It measures longer time period" )
    {
      {
        auto timer = Timer< CpuClock >( "timer-long" );
        std::this_thread::sleep_for( 1278912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< CpuClock >::get_duration_rep( "timer-long" );
        REQUIRE( d == Approx( 0 ).margin( 0.0001 ) );
      }
    }


    WHEN( "Resume an existing timer" )
    {
      auto pd = Timer< CpuClock >::get_duration_rep( "timer-long" );
      {
        auto timer = Timer< CpuClock >( "timer-long" );
        std::this_thread::sleep_for( 1000000us );
      }
      THEN( "It should get the accumulated duration" )
      {
        auto d = Timer< CpuClock >::get_duration_rep( "timer-long" );
        REQUIRE( d == Approx( pd ).margin( 0.0001 ) );
      }
    }
  }

  GIVEN( "A wall clock Timer" )
  {
    WHEN( "It measures a short time period" )
    {
      {
        auto timer = Timer< SteadyClock >( "timer-short" );
        std::this_thread::sleep_for( 678912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< SteadyClock >::get_duration_rep( "timer-short" );
        REQUIRE( static_cast< float >( d ) == Approx( 678912 ).epsilon( 0.001 ) );
      }
    }

    WHEN( "It measures longer time period" )
    {
      {
        auto timer = Timer< SteadyClock >( "timer-long" );
        std::this_thread::sleep_for( 1278912us );
      }
      THEN( "It should get the correct duration" )
      {
        auto d = Timer< SteadyClock >::get_duration_rep( "timer-long" );
        REQUIRE( static_cast< float >( d ) == Approx( 1278912 ).epsilon( 0.001 ) );
      }
    }

    WHEN( "Resume an existing timer" )
    {
      {
        auto timer = Timer< SteadyClock >( "timer-long" );
        std::this_thread::sleep_for( 1000000us );
      }
      THEN( "It should get the accumulated duration" )
      {
        auto d = Timer< SteadyClock >::get_duration_rep( "timer-long" );
        REQUIRE( static_cast< float >( d ) == Approx( 2278912 ).epsilon( 0.001 ) );
      }
    }
  }
}
