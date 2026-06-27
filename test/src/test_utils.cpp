/**
 *    @file  test_utils.cpp
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

#include <psi/utils.hpp>
#include <sdsl/bit_vectors.hpp>
#include <seqan/sequence.h>

#include "test_base.hpp"


using namespace psi;

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
    seqan2::CharString str( "mississipi" );
    seqan2::CharString pattern;

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
    seqan2::CharString str( "mississipi" );
    seqan2::CharString pattern;

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

SCENARIO( "Rounding up to the next closest power of 2", "[utils]" )
{
  GIVEN( "A 32-bit number" )
  {
    WHEN( "Rounding up a number which is NOT a power of 2" )
    {
      std::vector< uint32_t > numbers = { 0, 3, 15, 243, 65336, 65539, 2147483543, 2147483651, 4294967295 };
      std::vector< uint32_t > truth = { 1, 4, 16, 256, 65536, 131072, 2147483648, 0, 0 };
      THEN( "It should yield the next closest power of 2" )
      {
        for ( unsigned int i = 0; i < numbers.size(); ++i ) {
          REQUIRE( roundup32( numbers[ i ] ) == truth[ i ] );
        }
      }
    }

    WHEN( "Rounding up a number which is a power of 2" )
    {
      std::vector< uint32_t > numbers = { 1, 2, 4, 16, 256, 65536, 131072, 2147483648 };
      THEN( "The next closest power of 2 is itself" )
      {
        for ( unsigned int i = 0; i < numbers.size(); ++i ) {
          REQUIRE( roundup32( numbers[ i ] ) == numbers[ i ] );
        }
      }
    }
  }

  GIVEN( "A 64-bit number" )
  {
    WHEN( "Rounding up a number which is NOT a power of 2" )
    {
      std::vector< uint64_t > numbers = { 0, 3, 15, 243, 65336, 65539, 2147483543, 2147483651, 4611686018427387915u, 9223372036854775809u, 18446744073709551614u };
      std::vector< uint64_t > truth = { 1, 4, 16, 256, 65536, 131072, 2147483648, 4294967296, 9223372036854775808u, 0, 0 };
      THEN( "It should yield the next closest power of 2" )
      {
        for ( unsigned int i = 0; i < numbers.size(); ++i ) {
          REQUIRE( roundup64( numbers[ i ] ) == truth[ i ] );
        }
      }
    }

    WHEN( "Rounding up a number which is a power of 2" )
    {
      std::vector< uint64_t > numbers = { 1, 2, 4, 16, 256, 65536, 131072, 2147483648, 4611686018427387904u, 9223372036854775808u, 9223372036854775808u };
      THEN( "The next closest power of 2 is itself" )
      {
        for ( unsigned int i = 0; i < numbers.size(); ++i ) {
          REQUIRE( roundup64( numbers[ i ] ) == numbers[ i ] );
        }
      }
    }
  }
}

SCENARIO( "Serialize and deserialize a vector", "[utils]" )
{
  std::string file_name_prefix = test_data_dir + "/test_";

  size_t size = 20;
  GIVEN( "A vector of integer with size " + std::to_string( size ) )
  {
    std::vector< unsigned int > v;
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
        std::vector< unsigned int > w;
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
    std::vector< unsigned int > v;

    WHEN( "it is serialized to a file" )
    {
      std::string file_name = file_name_prefix + "2";
      std::ofstream ofs( file_name, std::ofstream::out | std::ofstream::binary );
      serialize( ofs, v.begin(), v.end() );
      ofs.close();
      THEN( "should be deserialized correctly" )
      {
        std::ifstream ifs( file_name, std::ifstream::in | std::ifstream::binary );
        std::vector< unsigned int > w;
        deserialize( ifs, w, std::back_inserter( w ) );
        REQUIRE( w.size() == 0 );
      }

      std::remove( file_name.c_str() );
    }
  }

  struct position {
    unsigned int i;
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
          i + 10,
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
          REQUIRE( w[i].c == static_cast< char >( i + 65 ) );
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

SCENARIO( "Find a value in an integer vector by reversal iteration", "[utils]" )
{
  GIVEN( "A high-dimensional `enc_vector` storing incremental integers" )
  {
    size_t len = 1000000;
    std::vector< unsigned int > v(len);
    for ( std::size_t i = 0; i < v.size(); ++i ) v[ i ] = i;
    sdsl::enc_vector< sdsl::coder::elias_delta<> > cv( v );

    WHEN( "The last item is searched" )
    {
      auto lc = rfind( cv, len - 1);

      THEN( "It should return an iterator pointing to the last item in the vector" )
      {
        REQUIRE( lc == cv.end() );
      }
    }

    WHEN( "Another item is searched" )
    {
      auto lc = rfind( cv, len - 10 );

      THEN( "It should return an iterator pointing to the that item in the vector" )
      {
        REQUIRE( lc == cv.end() - 9 );
      }
    }

    WHEN( "Non-existing item is searched" )
    {
      auto lc = rfind( cv, len );

      THEN( "It should return an iterator pointing to the that item in the vector" )
      {
        REQUIRE( lc == cv.begin() );
      }
    }
  }
}

SCENARIO( "Check equality of two vectors by reversal iteration", "[utils]" )
{
  GIVEN( "A `enc_vector` storing incremental integers and an existing query" )
  {
    size_t len = 1000000;
    std::vector< unsigned int > v(len);
    for ( std::size_t i = 0; i < v.size(); ++i ) v[ i ] = i;
    sdsl::enc_vector< sdsl::coder::elias_delta<> > cv( v );
    std::vector< unsigned int > query = { 999980, 999981, 999982, 999983, 999984, 999985 };
    auto lc = rfind( cv, *query.rbegin() );

    WHEN( "Check the equality by using reversed iteration" )
    {
      bool equal = requal( query.rbegin(), query.rend(), lc, cv.begin() );

      THEN( "The query should be found" )
      {
        REQUIRE( equal );
      }
    }
  }

  GIVEN( "A `enc_vector` storing incremental integers and an non-existing query" )
  {
    std::size_t len = 10;
    std::vector< unsigned int > v(len);
    for ( std::size_t i = 0; i < v.size(); ++i ) v[ i ] = i + 2;
    sdsl::enc_vector< sdsl::coder::elias_delta<> > cv( v );
    std::vector< unsigned int > query = { 0, 1, 2 };
    auto lc = rfind( cv, *query.rbegin() );

    WHEN( "Check the equality by using reversed iteration" )
    {
      bool equal = requal( query.rbegin(), query.rend(), lc, cv.begin() );

      THEN( "The query should not be found" )
      {
        REQUIRE( !equal );
      }
    }
  }
}

SCENARIO( "Word-wise range copy for bit-vectors", "[utils]" )
{
  GIVEN( "A bit vector whose size is less than a word length (64-bit)" )
  {
    sdsl::bit_vector sbv( 12, 0 );
    sbv[ 0 ] = 1;
    sbv[ 5 ] = 1;
    sbv[ 11 ] = 1;

    WHEN( "The whole bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 30, 1 );
      bv_icopy( sbv, dbv );

      THEN( "All first bits of length of source bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xfffffffffffff821 );
      }
    }

    WHEN( "The bit vector is partially copied to another" )
    {
      sdsl::bit_vector dbv( 30, 1 );
      bv_icopy( sbv, dbv, 6, 1 );

      THEN( "That part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffbf );
      }
    }

    WHEN( "A suffix of the bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 30, 1 );
      bv_icopy( sbv, dbv, 5 );

      THEN( "The suffix part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xfffffffffffff83f );
      }
    }
  }

  GIVEN( "A bit vector whose size is a multiplier of a word length (64-bit)" )
  {
    sdsl::bit_vector sbv( 7872, 0 );
    sbv.set_int( 542, 0x900000000fafabcd );
    sbv.set_int( 7808, 0xaaaaaaaaaaaaaaaa );

    WHEN( "The whole bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 7872, 1 );
      bv_icopy( sbv, dbv );

      THEN( "All first bits of length of source bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0x0 );
        REQUIRE( dbv.get_int( 100 ) == 0x0 );
        REQUIRE( dbv.get_int( 478 ) == 0x0 );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 893 ) == 0x0 );
        REQUIRE( dbv.get_int( 7744 ) == 0x0 );
        REQUIRE( dbv.get_int( 7808 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }

    WHEN( "The bit vector is partially copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bv_icopy( sbv, dbv, 542, 64 );

      THEN( "That part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
      }
    }

    WHEN( "A suffix of the bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bv_icopy( sbv, dbv, 7808 );

      THEN( "The suffix part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7744 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7808 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }
  }

  GIVEN( "A bit vector whose size is larger than word length (64-bit)" )
  {
    sdsl::bit_vector sbv( 7800, 0 );
    sbv.set_int( 0, 0xdddddddddddddddd );
    sbv.set_int( 542, 0x900000000fafabcd );
    sbv.set_int( 7736, 0xaaaaaaaaaaaaaaaa );

    WHEN( "The whole bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 7872, 1 );
      bv_icopy( sbv, dbv );

      THEN( "All first bits of length of source bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xdddddddddddddddd );
        REQUIRE( dbv.get_int( 100 ) == 0x0 );
        REQUIRE( dbv.get_int( 478 ) == 0x0 );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 893 ) == 0x0 );
        REQUIRE( dbv.get_int( 7672 ) == 0x0 );
        REQUIRE( dbv.get_int( 7736 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }

    WHEN( "The bit vector is partially copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bv_icopy( sbv, dbv, 542, 74 );

      THEN( "That part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 606 ) == 0xfffffffffffffc00 );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
      }
    }

    WHEN( "A suffix of the bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bv_icopy( sbv, dbv, 7736 );

      THEN( "The suffix part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7672 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7736 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }
  }

  GIVEN( "A bit vector whose size is less than a word length (64-bit)" )
  {
    sdsl::bit_vector sbv( 12, 0 );
    sbv[ 0 ] = 1;
    sbv[ 5 ] = 1;
    sbv[ 11 ] = 1;

    WHEN( "The whole bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 30, 1 );
      bvcopy( sbv, dbv );

      THEN( "All first bits of length of source bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xfffffffffffff821 );
      }
    }

    WHEN( "The bit vector is partially copied to another" )
    {
      sdsl::bit_vector dbv( 30, 1 );
      bvcopy( sbv, dbv, 6, 1, 6 );

      THEN( "That part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffbf );
      }
    }

    WHEN( "A suffix of the bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 30, 1 );
      bvcopy( sbv, dbv, 5, 7, 5 );

      THEN( "The suffix part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xfffffffffffff83f );
      }
    }
  }

  GIVEN( "A bit vector whose size is a multiplier of a word length (64-bit)" )
  {
    sdsl::bit_vector sbv( 7872, 0 );
    sbv.set_int( 542, 0x900000000fafabcd );
    sbv.set_int( 7808, 0xaaaaaaaaaaaaaaaa );

    WHEN( "The whole bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 7872, 1 );
      bvcopy( sbv, dbv );

      THEN( "All first bits of length of source bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0x0 );
        REQUIRE( dbv.get_int( 100 ) == 0x0 );
        REQUIRE( dbv.get_int( 478 ) == 0x0 );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 893 ) == 0x0 );
        REQUIRE( dbv.get_int( 7744 ) == 0x0 );
        REQUIRE( dbv.get_int( 7808 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }

    WHEN( "The bit vector is partially copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bvcopy( sbv, dbv, 542, 64, 542 );

      THEN( "That part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
      }
    }

    WHEN( "A suffix of the bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bvcopy( sbv, dbv, 7808, 64, 7808 );

      THEN( "The suffix part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7744 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7808 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }
  }

  GIVEN( "A bit vector whose size is larger than word length (64-bit)" )
  {
    sdsl::bit_vector sbv( 7800, 0 );
    sbv.set_int( 0, 0xdddddddddddddddd );
    sbv.set_int( 542, 0x900000000fafabcd );
    sbv.set_int( 7736, 0xaaaaaaaaaaaaaaaa );

    WHEN( "The whole bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 7872, 1 );
      bvcopy( sbv, dbv );

      THEN( "All first bits of length of source bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xdddddddddddddddd );
        REQUIRE( dbv.get_int( 100 ) == 0x0 );
        REQUIRE( dbv.get_int( 478 ) == 0x0 );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 893 ) == 0x0 );
        REQUIRE( dbv.get_int( 7672 ) == 0x0 );
        REQUIRE( dbv.get_int( 7736 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }

    WHEN( "The bit vector is partially copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bvcopy( sbv, dbv, 542, 74, 542 );

      THEN( "That part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0x900000000fafabcd );
        REQUIRE( dbv.get_int( 606 ) == 0xfffffffffffffc00 );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
      }
    }

    WHEN( "A suffix of the bit vector is copied to another" )
    {
      sdsl::bit_vector dbv( 8000, 1 );
      bvcopy( sbv, dbv, 7736, 64, 7736 );

      THEN( "The suffix part of destination bit vector should be equal to the source" )
      {
        REQUIRE( dbv.get_int( 0 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 100 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 478 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 542 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 893 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 6936 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7672 ) == 0xffffffffffffffff );
        REQUIRE( dbv.get_int( 7736 ) == 0xaaaaaaaaaaaaaaaa );
      }
    }
  }
}

SCENARIO( "Generate random real numbers", "[utils][random]" )
{
  WHEN( "No range is given" )
  {
    int n = 100000;

    THEN( "Generated random real number should be in range [0, 1)" )
    {
      for ( int i = 0; i < n; ++i ) {
        auto rnum = random::random_real< double >();
        REQUIRE( ( 0 <= rnum && rnum <= 1 ) );
      }
    }
  }

  WHEN( "A small range is given" )
  {
    float low = 10;
    float high = 100;
    int n = 100000;

    THEN( "Generated random real numbers should be in the range" )
    {
      for ( int i = 0; i < n; ++i ) {
        auto rnum = random::random_real< float >( low, high );
        REQUIRE( ( low <= rnum && rnum <= high ) );
      }
    }
  }

  WHEN( "A large range is given" )
  {
    double low = std::numeric_limits< double >::min();
    double high = std::numeric_limits< double >::max();
    int n = 100000;

    THEN( "Generated random real numbers should be in the range" )
    {
      for ( int i = 0; i < n; ++i ) {
        auto rnum = random::random_real( low, high );
        REQUIRE( ( low <= rnum && rnum <= high ) );
      }
    }
  }
}

SCENARIO( "Generate random integers", "[utils][random]" )
{
  WHEN( "A range is given" )
  {
    int low = random::random_integer< int >();
    int high = random::random_integer< int >( low + 1 );
    int n = 1000000;

    THEN( "Generated random integers should be in the range" )
    {
      for ( int i = 0; i < n; ++i ) {
        auto rnum = random::random_integer( low, high );
        REQUIRE( ( low <= rnum && rnum <= high ) );
      }
    }
  }

  WHEN( "A length is given" )
  {
    unsigned int len = random::random_integer< unsigned int >();
    int n = 1000000;

    THEN( "Generated random index should be in the range" )
    {
      for ( int i = 0; i < n; ++i ) {
        auto rnum = random::random_index( len );
        REQUIRE( ( 0 <= rnum && rnum <= len - 1 ) );
      }
    }
  }
}

SCENARIO( "Generate random strings", "[utils][random]" )
{
  WHEN( "A length is given" )
  {
    unsigned int len = 6060;
    int n = 1000;

    THEN( "Generated random string should have the given length" )
    {
      for ( int i = 0; i < n; ++i ) {
        REQUIRE( random::random_string( len ).length() == len );
      }
    }
  }
}

SCENARIO( "Compute the average of a long stream of integers in parallel", "[utils]" )
{
  using value_type = uint16_t;

  uint8_t nofthreads = 8;
  value_type lbound = 20;
  value_type ubound = 50;
  GIVEN( std::to_string( nofthreads ) + " threads contributing to a shared sum variable" )
  {
    std::atomic< value_type > sum( 0 );
    std::atomic< value_type > total( 0 );
    std::atomic< double > avg( -1 );
    RWSpinLock rws_lock;
    std::atomic_ullong real_sum( 0 );
    std::atomic_ullong real_tot( 0 );
    std::vector< std::thread > threads;

    auto update =
        [&sum, &total, &avg]() -> void {
          auto partial_sum = sum.load();
          auto partial_total = total.load();
          sum.store( 0 );
          total.store( 0 );
          double pre_avg = avg.load();
          double new_avg = partial_sum / static_cast< double >( partial_total );
          if ( pre_avg == -1 ) avg.store( new_avg );
          else avg.store( ( new_avg + pre_avg ) / 2 );
        };

    auto run =
        [&sum, &total, &rws_lock, &real_sum, &real_tot, update, lbound, ubound]() -> void {
          constexpr const unsigned int RETRY_THRESHOLD = 4;
          constexpr const unsigned int NOF_VALUES = 2000;

          unsigned int retry = RETRY_THRESHOLD;
          for ( std::size_t i = 0; i < NOF_VALUES; ++i ) {
            value_type value = random::random_integer( lbound, ubound );
            while ( true ) {
              auto peek_sum = sum.load();
              if ( peek_sum >= std::numeric_limits< value_type >::max() - value ) {
                UniqWriterLock reducer( rws_lock );
                if ( reducer && peek_sum == sum.load() ) {
                  REQUIRE( !rws_lock.acquire_writer_weak() );
                  REQUIRE( sum.load() >= std::numeric_limits< value_type >::max() - value );
                  update();
                }
                continue;
              }
              {
                ReaderLock adder( rws_lock );
                if ( sum.compare_exchange_weak( peek_sum, peek_sum+value,
                                                std::memory_order_release, std::memory_order_relaxed ) ) {
                  total.fetch_add( 1 );
                  break;
                }
              }

              if ( --retry == 0 ) {
                retry = RETRY_THRESHOLD;
                std::this_thread::yield();
              }
            }
            real_sum.fetch_add( value );
            real_tot.fetch_add( 1 );
          }
        };

    WHEN( "The average is computed" )
    {
      for ( std::size_t tidx = 0; tidx < nofthreads; ++tidx ) threads.emplace_back( run );
      THEN( "It should be approximately correct" )
      {
        for ( std::size_t tidx = 0; tidx < nofthreads; ++tidx ) threads[ tidx ].join();
        update();
        double diff = real_sum / static_cast< double >( real_tot ) - avg;
        REQUIRE( diff == Approx( 0 ).margin( 0.5 ) );
      }
    }
  }
}
