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

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>

#include <seqan/sequence.h>

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

  GIVEN( "A seqan string" )
  {
    seqan::CharString str( "mississipi" );
    seqan::CharString pattern;

    THEN( "Suffix strings should be matched" )
    {
      pattern = "pi";
      REQUIRE( ends_with( str, pattern ) );
      pattern = "issipi";
      REQUIRE( ends_with( str, pattern ) );
      pattern = "";
      REQUIRE( ends_with( str, pattern ) );
      pattern = "mississipi";
      REQUIRE( ends_with( str, pattern ) );
    }

    THEN( "Non-suffix strings should not be matched")
    {
      pattern = "m";
      REQUIRE( !ends_with( str, pattern ) );
      pattern = "missi";
      REQUIRE( !ends_with( str, pattern ) );
      pattern = "issi";
      REQUIRE( !ends_with( str, pattern ) );
      pattern = "MISSISSIPI";
      REQUIRE( !ends_with( str, pattern ) );
      pattern = "I";
      REQUIRE( !ends_with( str, pattern ) );
      pattern = "arizona";
      REQUIRE( !ends_with( str, pattern ) );
    }
  }
}

SCENARIO( "Two strings can be checked for prefix match", "[utils]" )
{
  GIVEN( "A string" )
  {
    std::string str( "mississipi" );

    THEN( "Prefix strings should be matched" )
    {
      REQUIRE( starts_with( str, "mi" ) );
      REQUIRE( starts_with( str, "missis" ) );
      REQUIRE( starts_with( str, "" ) );
      REQUIRE( starts_with( str, "mississipi" ) );
    }

    THEN( "Non-prefix strings should not be matched")
    {
      REQUIRE( !starts_with( str, "i" ) );
      REQUIRE( !starts_with( str, "ssipi" ) );
      REQUIRE( !starts_with( str, "issi" ) );
      REQUIRE( !starts_with( str, "MISSISSIPI" ) );
      REQUIRE( !starts_with( str, "I" ) );
      REQUIRE( !starts_with( str, "arizona" ) );
    }
  }

  GIVEN( "A seqan string" )
  {
    seqan::CharString str( "mississipi" );
    seqan::CharString pattern;

    THEN( "Prefix strings should be matched" )
    {
      pattern = "mi";
      REQUIRE( starts_with( str, pattern ) );
      pattern = "missis";
      REQUIRE( starts_with( str, pattern ) );
      pattern = "";
      REQUIRE( starts_with( str, pattern ) );
      pattern = "mississipi";
      REQUIRE( starts_with( str, pattern ) );
    }

    THEN( "Non-prefix strings should not be matched")
    {
      pattern = "i";
      REQUIRE( !starts_with( str, pattern ) );
      pattern = "ssipi";
      REQUIRE( !starts_with( str, pattern ) );
      pattern = "issi";
      REQUIRE( !starts_with( str, pattern ) );
      pattern = "MISSISSIPI";
      REQUIRE( !starts_with( str, pattern ) );
      pattern = "I";
      REQUIRE( !starts_with( str, pattern ) );
      pattern = "arizona";
      REQUIRE( !starts_with( str, pattern ) );
    }
  }
}

SCENARIO( "Serialize and deserialize a vector", "[utils]" )
{
  std::string file_name_prefix = _testdir + "/test_";

  size_t size = 20;
  GIVEN( "A vector of integer with size " + std::to_string( size ) )
  {
    std::vector< int > v;
    for ( unsigned int i = 0; i < size; ++i ) {
      v.push_back( i * 2 );
    }

    WHEN( "it is serialized to a file" )
    {
      std::string file_name = file_name_prefix + "1";
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      serialize( ofs, v.begin(), v.end() );
      ofs.close();
      THEN( "should be deserialized correctly" )
      {
        std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
        std::vector< int > w;
        deserialize( ifs, w, std::back_inserter( w ) );

        for ( unsigned int i = 0; i < w.size(); ++i ) {
          REQUIRE( w[i] == i * 2 );
        }
      }

      std::remove( file_name.c_str() );
    }
  }

  GIVEN( "An empty string" )
  {
    std::vector< int > v;

    WHEN( "it is serialized to a file" )
    {
      std::string file_name = file_name_prefix + "2";
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      serialize( ofs, v.begin(), v.end() );
      ofs.close();
      THEN( "should be deserialized correctly" )
      {
        std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
        std::vector< int > w;
        deserialize( ifs, w, std::back_inserter( w ) );
        REQUIRE( w.size() == 0 );
      }

      std::remove( file_name.c_str() );
    }
  }

  struct position {
    int i;
    double d;
    char c;
    char s[10];
  };

  size = 10;
  GIVEN( "A vector of a structure with size " + std::to_string( size ) )
  {
    std::vector< position > v;
    for ( unsigned int i = 0; i < size; ++i ) {
      v.push_back( {
          static_cast<int>( i + 10 ),
          i / 3.0,
          static_cast<char>( i + 65 ),
          { static_cast<char>( i + 65 ), static_cast<char>( i + 97 ) }
          } );
    }

    WHEN( "it is serialized to a file" )
    {
      std::string file_name = file_name_prefix + "3";
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      serialize( ofs, v.begin(), v.end() );
      ofs.close();
      THEN( "should be deserialized correctly" )
      {
        std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
        std::vector< position > w;
        deserialize( ifs, w, std::back_inserter( w ) );

        for ( unsigned int i = 0; i < w.size(); ++i ) {
          REQUIRE( w[i].i == i + 10 );
          REQUIRE( w[i].d == i / 3.0 );
          REQUIRE( w[i].c == i + 65 );
          const char s[10] = { static_cast<char>( i + 65 ), static_cast<char>( i + 97 ) };
          REQUIRE( std::strcmp( w[i].s, s ) == 0 );
        }
      }

      std::remove( file_name.c_str() );
    }
  }
}

SCENARIO( "Check if a file path is readable", "[utils]" )
{
  GIVEN( "A file name that exists and is readable" )
  {
    std::string tmpfpath = SEQAN_TEMP_FILENAME();
    std::ofstream ofs( tmpfpath );
    ofs.close();

    THEN( "It should be readable" )
    {
      REQUIRE( readable( tmpfpath ) );
    }
  }

  GIVEN( "A file name that does not exist" )
  {
    std::string tmpfpath = SEQAN_TEMP_FILENAME();

    THEN( "It should not be readable" )
    {
      REQUIRE( ! readable( tmpfpath ) );
    }
  }

  GIVEN( "A file name that exists, but is not readable" )
  {
    std::string filepath = "/root/.Xauthority";

    THEN( "It should not be readable" )
    {
      REQUIRE( ! readable( filepath ) );
    }
  }
}

SCENARIO( "Check if a file path is writable", "[utils]" )
{
  GIVEN( "A file name that exists and is writable" )
  {
    std::string tmpfpath = SEQAN_TEMP_FILENAME();
    std::ofstream ofs( tmpfpath );
    ofs.close();

    THEN( "It should be writable" )
    {
      REQUIRE( writable( tmpfpath ) );
    }
  }

  GIVEN( "A file name that does not exist" )
  {
    std::string tmpfpath = SEQAN_TEMP_FILENAME();

    THEN( "It should be writable" )
    {
      REQUIRE( writable( tmpfpath ) );
      REQUIRE( ! readable( tmpfpath ) );
    }
  }

  GIVEN( "A file name that exists, but is not writable" )
  {
    std::string filepath = "/root/.Xauthority";

    THEN( "It should not be writable" )
    {
      REQUIRE( ! writable( filepath ) );
    }
  }
}

SCENARIO( "Check file appendability", "[utils]" )
{
  GIVEN( "A file name that exists and is appendable" )
  {
    std::string tmpfpath = SEQAN_TEMP_FILENAME();
    std::ofstream ofs( tmpfpath );
    ofs.close();

    THEN( "It should be appendable" )
    {
      REQUIRE( appendable( tmpfpath ) );
    }
  }

  GIVEN( "A file name that does not exist" )
  {
    std::string tmpfpath = SEQAN_TEMP_FILENAME();

    THEN( "It should not be appendable" )
    {
      REQUIRE( ! appendable( tmpfpath ) );
    }
  }

  GIVEN( "A file name that exists, but is not appendable" )
  {
    std::string filepath = "/root/.Xauthority";

    THEN( "It should not be appendable" )
    {
      REQUIRE( ! appendable( filepath ) );
    }
  }
}
