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

#include<Kokkos_Core.hpp>

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include "test_base.hpp"

class MyTestEventListener : public Catch::TestEventListenerBase {
public:
    using Catch::TestEventListenerBase::TestEventListenerBase;

    void set_rnd_seed() {
      auto seed = Catch::rngSeed();
      if ( seed != 0 ) {
        std::cout << "Setting random generator seed to " << seed << "..." << std::endl;
        rnd::set_seed( seed );
      }
    }

    void testRunStarting( Catch::TestRunInfo const& ) override {
      this->set_rnd_seed();
      /*
      // For debugging
      Kokkos::InitializationSettings args;
      args.set_num_threads( 1 );
      Kokkos::initialize( args );
      */
      Kokkos::initialize();
    }

    void testRunEnded( Catch::TestRunStats const& testRunStats ) override {
      Kokkos::finalize();
    }
};

CATCH_REGISTER_LISTENER(MyTestEventListener)
