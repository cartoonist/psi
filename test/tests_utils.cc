/**
 *    @file  tests_utils.cc
 *   @brief  utils test cases.
 *
 *  Contains test cases for `utils.h` header file.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Aug 20, 2017  16:38
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include "tests_base.h"
#include "utils.h"


using namespace grem;

SCENARIO( "Two strings can be checked for suffix match", "[utils]" )
{
  GIVEN( "A string" )
  {
    std::string str( "mississipi" );

    THEN( "Suffix strings should be matched" )
    {
      REQUIRE( ends_with( str, "pi" ) );
      REQUIRE( ends_with( str, "issipi" ) );
      REQUIRE( ends_with( str, "" ) );
      REQUIRE( ends_with( str, "mississipi" ) );
    }

    THEN( "Non-suffix strings should not be matched")
    {
      REQUIRE( !ends_with( str, "m" ) );
      REQUIRE( !ends_with( str, "missi" ) );
      REQUIRE( !ends_with( str, "issi" ) );
      REQUIRE( !ends_with( str, "MISSISSIPI" ) );
      REQUIRE( !ends_with( str, "I" ) );
      REQUIRE( !ends_with( str, "arizona" ) );
    }
  }
}
