/**
 *    @file  test_main.cpp
 *   @brief  Unit test framework (Catch2) implementation file.
 *
 *  In order to compile the unit testing framework once, not for each tranlation unit
 *  which includes the header file, we put the `CATCH_CONFIG_MAIN` macro in a separated
 *  source file while the header file is included. This file should be linked with other
 *  test implementations.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Feb 13, 2017  00:44
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <iostream>

#include <Kokkos_Core.hpp>

#include "catch2/catch_session.hpp"
#include "catch2/catch_get_random_seed.hpp"
#include "catch2/reporters/catch_reporter_event_listener.hpp"
#include "catch2/reporters/catch_reporter_registrars.hpp"

#include "test_base.hpp"

class MyTestEventListener : public Catch::EventListenerBase {
public:
    using Catch::EventListenerBase::EventListenerBase;

    void set_rnd_seed() {
      auto seed = Catch::getSeed();
      if ( seed != 0 ) {
        std::cout << "Setting random generator seed to " << seed << "..." << std::endl;
        rnd::set_seed( seed );
      }
    }

    void testRunStarting( Catch::TestRunInfo const& ) override {
      this->set_rnd_seed();
    }
};

CATCH_REGISTER_LISTENER(MyTestEventListener)

int main( int argc, char* argv[] )
{
  Kokkos::initialize( argc, argv );

  int result = Catch::Session().run( argc, argv );

  Kokkos::finalize();

  return result;
}
