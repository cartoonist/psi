/**
 *    @file  test_sequence.cpp
 *   @brief  Test sequence module.
 *
 *  Test scenarios for sequence module.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Oct 04, 2017  02:13
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <string>

#include <seqan/seq_io.h>

#include "sequence.h"

#include "test_base.hpp"


using namespace grem;


SCENARIO( "Subsetting a reads chunk from a reads set", "[sequence]" )
{
  unsigned int reads_num = 10;
  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = test_data_dir + "/small/reads_n10l10e0i0.fastq";
    seqan::SeqFileIn infile;
    if ( !open( infile, fqpath.c_str() ) ) {
      throw std::runtime_error( "cannot open file '" + fqpath + "'" );
    }
    Records< Dna5QStringSet<> > reads;
    readRecords( reads, infile );
    Records< Dna5QStringSet< grem::Dependent > > reads_chunk;
    unsigned int subset_len = 4;
    unsigned int offset = 2;

    WHEN( "A " + std::to_string( subset_len ) + "-reads subset is selected" )
    {
      load_chunk( reads_chunk, reads, subset_len, offset );
      THEN( "Each read's ID in the chunk should be its position in the original set" )
      {
        REQUIRE( length( reads_chunk ) == subset_len );
        for ( unsigned int i = 0; i < subset_len; ++i ) {
          REQUIRE( position_to_id( reads_chunk, i ) == offset + i );
        }
      }

      WHEN( "Next " + std::to_string( subset_len ) + "-reads subset is selected" )
      {
        auto last_id = position_to_id( reads_chunk, length( reads_chunk ) - 1 );
        load_chunk( reads_chunk, reads, subset_len, last_id + 1 );
        THEN( "Each read's ID in the chunk should be its position in the original set" )
        {
          REQUIRE( length( reads_chunk ) == subset_len );
          for ( unsigned int i = 0; i < subset_len; ++i ) {
            REQUIRE( position_to_id( reads_chunk, i ) == offset + subset_len + i );
          }
        }
      }
    }
  }
}

SCENARIO( "Load reads to an owner Records with non-zero offset", "[sequence]" )
{
  unsigned int reads_num = 10;
  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = test_data_dir + "/small/reads_n10l10e0i0.fastq";
    klibpp::SeqStreamIn iss( fqpath.c_str() );
    Records< seqan::StringSet< MemString > > records;
    unsigned int subset_len = 4;
    unsigned int offset = 2;
    readRecords( records, iss, offset );

    WHEN( std::to_string( subset_len ) + " reads are parsed from file" )
    {
      readRecords( records, iss, subset_len );

      THEN( "Each read's ID in the chunk should be its position in the original set" )
      {
        REQUIRE( length( records ) == subset_len );
        for ( unsigned int i = 0; i < subset_len; ++i ) {
          REQUIRE( position_to_id( records, i ) == offset + i );
        }
      }

      WHEN( "Next " + std::to_string( subset_len ) + " reads are parsed" )
      {
        readRecords( records, iss, subset_len );

        THEN( "Each read's ID in the chunk should be its position in the original set" )
        {
          REQUIRE( length( records.str ) == subset_len );
          for ( unsigned int i = 0; i < subset_len; ++i ) {
            REQUIRE( position_to_id( records, i ) == offset + subset_len + i );
          }
        }
      }
    }
  }
}

SCENARIO( "Constructing a DiskString", "[sequence]" )
{
  auto check_content =
    []( DiskString& d, const std::string& data ) {
      std::ifstream in( d.get_file_path() );
      std::string content;
      std::string buffer;
      while ( std::getline( in, buffer ) ) content += buffer;
      REQUIRE( data == content );
    };

  GIVEN( "A simple text" )
  {
    std::string text( "a mississippian lazy fox sits on a pie" );

    WHEN( "DiskString is constructed by providing data to constructor" )
    {
      DiskString dstr( text );

      THEN( "It should contain the data" )
      {
        check_content( dstr, text );
        REQUIRE( length( dstr ) == 38 );
      }
    }

    WHEN( "DiskString is constructed by providing data to constructor with C-style string" )
    {
      DiskString dstr( text.c_str() );

      THEN( "It should contain the data" )
      {
        check_content( dstr, text );
        REQUIRE( dstr.length() == 38 );
      }
    }

    WHEN( "Some strings are appended to DiskString (random file)")
    {
      DiskString dstr;
      dstr.reserve( 38 );
      dstr = text.substr( 0, 15 );
      dstr += text.substr( 15, 10 );
      dstr += text.substr( 25 ).c_str();

      THEN( "The value should be the concatenation of the strings" )
      {
        check_content( dstr, text );
        REQUIRE( dstr.length() == 38 );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstr, another_tmpfile );
        clear( dstr );
        open( dstr, another_tmpfile );

        THEN( "The value should be the concatenation of the strings" )
        {
          check_content( dstr, text );
          REQUIRE( dstr.length() == 38 );
          REQUIRE( readable( another_tmpfile ) );
        }
      }
    }

    WHEN( "Some strings are appended to DiskString (specific file)")
    {
      std::string tmpfile = SEQAN_TEMP_FILENAME();
      DiskString dstr( "", tmpfile );
      dstr.reserve( 38 );
      dstr = text.substr( 0, 15 );
      dstr += text.substr( 15, 10 );
      dstr += text.substr( 25 ).c_str();

      THEN( "The value should be the concatenation of the strings" )
      {
        check_content( dstr, text );
        REQUIRE( dstr.length() == 38 );
        REQUIRE( readable( tmpfile ) );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstr, another_tmpfile );
        clear( dstr );
        open( dstr, another_tmpfile );

        THEN( "The value should be the concatenation of the strings" )
        {
          check_content( dstr, text );
          REQUIRE( dstr.length() == 38 );
          REQUIRE( readable( another_tmpfile ) );
        }
      }
    }

    WHEN( "Some strings are appended to DiskString (specific file with C-Style string)")
    {
      std::string tmpfile = SEQAN_TEMP_FILENAME();
      DiskString dstr( text, tmpfile.c_str() );

      THEN( "The value should be the concatenation of the strings" )
      {
        check_content( dstr, text );
        REQUIRE( dstr.length() == 38 );
        REQUIRE( readable( tmpfile ) );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstr, another_tmpfile );
        clear( dstr );
        open( dstr, another_tmpfile );

        THEN( "The value should be the concatenation of the strings" )
        {
          check_content( dstr, text );
          REQUIRE( dstr.length() == 38 );
          REQUIRE( readable( another_tmpfile ) );
        }
      }
    }

    GIVEN( "An open DiskString containing the text" )
    {
      DiskString dstr( text );

      WHEN( "It is moved to another DiskString" )
      {
        DiskString another = std::move( dstr );

        THEN( "The new owner should get the same content" )
        {
          check_content( another, text );
          REQUIRE( length( dstr ) == 38 );
        }

        AND_WHEN( "It is appended by another string after move" )
        {
          std::string appending_text = " suddenly!";
          another += appending_text;

          THEN( "The new owner should get the same content appended by that string" )
          {
            auto new_text = text + appending_text;
            check_content( another, new_text );
            REQUIRE( length( another ) == 48 );
          }

          AND_WHEN( "It is written to file and loaded again" )
          {
            std::string another_tmpfile = SEQAN_TEMP_FILENAME();
            save( another, another_tmpfile );
            clear( another );
            open( another, another_tmpfile );

            THEN( "The value should be the concatenation of the strings" )
            {
              auto new_text = text + appending_text;
              check_content( another, new_text );
              REQUIRE( another.length() == 48 );
              REQUIRE( readable( another_tmpfile ) );
            }
          }
        }
      }

      WHEN( "It is closed by getting the file path" )
      {
        REQUIRE( dstr.is_open() );
        auto fpath = dstr.get_file_path();
        REQUIRE( !dstr.is_open() );

        THEN( "Adding another string raises an exception" )
        {
          REQUIRE_THROWS_AS( dstr.operator+=( text ), std::runtime_error );
        }
      }

      WHEN( "It is cleared then new text is added" )
      {
        std::string new_text = "another brazilian cute beaver builds a dam";
        clear( dstr );
        REQUIRE( length( dstr ) == 0 );
        dstr += new_text;

        THEN( "It should hold the new text (overwriting the old one)" )
        {
          check_content( dstr, new_text );
          REQUIRE( dstr.length() == 42 );
        }
      }
    }
  }
}

SCENARIO( "Constructing a set of disk-based string", "[sequence]" )
{
  auto check_content =
    []( seqan::StringSet< DiskString >& d, const std::string& data ) {
      std::ifstream in( d.get_file_path() );
      std::string content;
      std::string buffer;
      while ( std::getline( in, buffer ) ) content += buffer;
      REQUIRE( data == content );
    };

  GIVEN( "A set of strings" )
  {
    std::string str1 = "a mississippian lazy fox sits on a pie";
    std::string str2 = "another brazilian cute beaver builds a dam";
    std::string str3 = "some african stupid chimps eat banana";
    std::string raw_total =
      str1 + std::string( 1, SEQUENCE_DEFAULT_SENTINEL_CHAR ) +
      str2 + std::string( 1, SEQUENCE_DEFAULT_SENTINEL_CHAR ) + str3;

    WHEN( "The strings are added to string set by const reference" )
    {
      seqan::StringSet< DiskString, seqan::Owner<> > dstrset;
      dstrset.reserve( 3 );
      appendValue( dstrset, str1 );
      push_back( dstrset, str2 );
      dstrset.push_back( str3 );

      THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
      {
        check_content( dstrset, raw_total );
        REQUIRE( dstrset.length() == 3 );
        REQUIRE( length( dstrset[ 0 ] ) == 38 );
        REQUIRE( length( dstrset[ 1 ] ) == 42 );
        REQUIRE( length( dstrset[ 2 ] ) == 37 );
        REQUIRE( dstrset.get_id( 0 ) == 0 );
        REQUIRE( dstrset.get_offset( 0 ) == 0 );
        REQUIRE( dstrset.get_id( 25 ) == 0 );
        REQUIRE( dstrset.get_offset( 25 ) == 25 );
        REQUIRE( dstrset.get_id( 37 ) == 0 );
        REQUIRE( dstrset.get_offset( 37 ) == 37 );
        REQUIRE( dstrset.get_id( 39 ) == 1 );
        REQUIRE( dstrset.get_offset( 39 ) == 0 );
        REQUIRE( dstrset.get_id( 51 ) == 1 );
        REQUIRE( dstrset.get_offset( 51 ) == 12 );
        REQUIRE( dstrset.get_id( 80 ) == 1 );
        REQUIRE( dstrset.get_offset( 80 ) == 41 );
        REQUIRE( dstrset.get_id( 82 ) == 2 );
        REQUIRE( dstrset.get_offset( 82 ) == 0 );
        REQUIRE( dstrset.get_id( 100 ) == 2 );
        REQUIRE( dstrset.get_offset( 100 ) == 18 );
        REQUIRE( dstrset.get_id( 118 ) == 2 );
        REQUIRE( dstrset.get_offset( 118 ) == 36 );
      }
    }

    WHEN( "The strings are added to string set by const reference (specific file)" )
    {
      std::string tmpfile = SEQAN_TEMP_FILENAME();
      seqan::StringSet< DiskString, seqan::Owner<> > dstrset( tmpfile );
      dstrset.reserve( 3 );
      appendValue( dstrset, str1 );
      push_back( dstrset, str2 );
      dstrset.push_back( str3 );

      THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
      {
        check_content( dstrset, raw_total );
        REQUIRE( dstrset.length() == 3 );
        REQUIRE( length( dstrset[ 0 ] ) == 38 );
        REQUIRE( length( dstrset[ 1 ] ) == 42 );
        REQUIRE( length( dstrset[ 2 ] ) == 37 );
        REQUIRE( dstrset.get_id( 0 ) == 0 );
        REQUIRE( dstrset.get_offset( 0 ) == 0 );
        REQUIRE( dstrset.get_id( 25 ) == 0 );
        REQUIRE( dstrset.get_offset( 25 ) == 25 );
        REQUIRE( dstrset.get_id( 37 ) == 0 );
        REQUIRE( dstrset.get_offset( 37 ) == 37 );
        REQUIRE( dstrset.get_id( 39 ) == 1 );
        REQUIRE( dstrset.get_offset( 39 ) == 0 );
        REQUIRE( dstrset.get_id( 51 ) == 1 );
        REQUIRE( dstrset.get_offset( 51 ) == 12 );
        REQUIRE( dstrset.get_id( 80 ) == 1 );
        REQUIRE( dstrset.get_offset( 80 ) == 41 );
        REQUIRE( dstrset.get_id( 82 ) == 2 );
        REQUIRE( dstrset.get_offset( 82 ) == 0 );
        REQUIRE( dstrset.get_id( 100 ) == 2 );
        REQUIRE( dstrset.get_offset( 100 ) == 18 );
        REQUIRE( dstrset.get_id( 118 ) == 2 );
        REQUIRE( dstrset.get_offset( 118 ) == 36 );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstrset, another_tmpfile );
        clear( dstrset );
        open( dstrset, another_tmpfile );

        THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
        {
          check_content( dstrset, raw_total );
          REQUIRE( dstrset.length() == 3 );
          REQUIRE( length( dstrset[ 0 ] ) == 38 );
          REQUIRE( length( dstrset[ 1 ] ) == 42 );
          REQUIRE( length( dstrset[ 2 ] ) == 37 );
          REQUIRE( dstrset.get_id( 0 ) == 0 );
          REQUIRE( dstrset.get_offset( 0 ) == 0 );
          REQUIRE( dstrset.get_id( 25 ) == 0 );
          REQUIRE( dstrset.get_offset( 25 ) == 25 );
          REQUIRE( dstrset.get_id( 37 ) == 0 );
          REQUIRE( dstrset.get_offset( 37 ) == 37 );
          REQUIRE( dstrset.get_id( 39 ) == 1 );
          REQUIRE( dstrset.get_offset( 39 ) == 0 );
          REQUIRE( dstrset.get_id( 51 ) == 1 );
          REQUIRE( dstrset.get_offset( 51 ) == 12 );
          REQUIRE( dstrset.get_id( 80 ) == 1 );
          REQUIRE( dstrset.get_offset( 80 ) == 41 );
          REQUIRE( dstrset.get_id( 82 ) == 2 );
          REQUIRE( dstrset.get_offset( 82 ) == 0 );
          REQUIRE( dstrset.get_id( 100 ) == 2 );
          REQUIRE( dstrset.get_offset( 100 ) == 18 );
          REQUIRE( dstrset.get_id( 118 ) == 2 );
          REQUIRE( dstrset.get_offset( 118 ) == 36 );
        }
      }
    }

    WHEN( "The strings are added to string set by const reference (specific file with C-Style string)" )
    {
      std::string tmpfile = SEQAN_TEMP_FILENAME();
      seqan::StringSet< DiskString, seqan::Owner<> > dstrset( tmpfile.c_str() );
      dstrset.reserve( 3 );
      appendValue( dstrset, str1 );
      push_back( dstrset, str2 );
      dstrset.push_back( str3 );

      THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
      {
        check_content( dstrset, raw_total );
        REQUIRE( dstrset.length() == 3 );
        REQUIRE( length( dstrset[ 0 ] ) == 38 );
        REQUIRE( length( dstrset[ 1 ] ) == 42 );
        REQUIRE( length( dstrset[ 2 ] ) == 37 );
        REQUIRE( dstrset.get_id( 0 ) == 0 );
        REQUIRE( dstrset.get_offset( 0 ) == 0 );
        REQUIRE( dstrset.get_id( 25 ) == 0 );
        REQUIRE( dstrset.get_offset( 25 ) == 25 );
        REQUIRE( dstrset.get_id( 37 ) == 0 );
        REQUIRE( dstrset.get_offset( 37 ) == 37 );
        REQUIRE( dstrset.get_id( 39 ) == 1 );
        REQUIRE( dstrset.get_offset( 39 ) == 0 );
        REQUIRE( dstrset.get_id( 51 ) == 1 );
        REQUIRE( dstrset.get_offset( 51 ) == 12 );
        REQUIRE( dstrset.get_id( 80 ) == 1 );
        REQUIRE( dstrset.get_offset( 80 ) == 41 );
        REQUIRE( dstrset.get_id( 82 ) == 2 );
        REQUIRE( dstrset.get_offset( 82 ) == 0 );
        REQUIRE( dstrset.get_id( 100 ) == 2 );
        REQUIRE( dstrset.get_offset( 100 ) == 18 );
        REQUIRE( dstrset.get_id( 118 ) == 2 );
        REQUIRE( dstrset.get_offset( 118 ) == 36 );
      }
    }

    WHEN( "The strings are added to string set by r-value reference" )
    {
      seqan::StringSet< DiskString, seqan::Owner<> > dstrset;
      reserve( dstrset, 3 );
      appendValue( dstrset, std::move( str1 ) );
      push_back( dstrset, std::move( str2 ) );
      dstrset.push_back( std::move( str3 ) );

      THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
      {
        check_content( dstrset, raw_total );
        REQUIRE( length( dstrset ) == 3 );
        REQUIRE( length( dstrset[ 0 ] ) == 38 );
        REQUIRE( length( dstrset[ 1 ] ) == 42 );
        REQUIRE( length( dstrset[ 2 ] ) == 37 );
        REQUIRE( dstrset.get_id( 0 ) == 0 );
        REQUIRE( dstrset.get_offset( 0 ) == 0 );
        REQUIRE( dstrset.get_id( 25 ) == 0 );
        REQUIRE( dstrset.get_offset( 25 ) == 25 );
        REQUIRE( dstrset.get_id( 37 ) == 0 );
        REQUIRE( dstrset.get_offset( 37 ) == 37 );
        REQUIRE( dstrset.get_id( 39 ) == 1 );
        REQUIRE( dstrset.get_offset( 39 ) == 0 );
        REQUIRE( dstrset.get_id( 51 ) == 1 );
        REQUIRE( dstrset.get_offset( 51 ) == 12 );
        REQUIRE( dstrset.get_id( 80 ) == 1 );
        REQUIRE( dstrset.get_offset( 80 ) == 41 );
        REQUIRE( dstrset.get_id( 82 ) == 2 );
        REQUIRE( dstrset.get_offset( 82 ) == 0 );
        REQUIRE( dstrset.get_id( 100 ) == 2 );
        REQUIRE( dstrset.get_offset( 100 ) == 18 );
        REQUIRE( dstrset.get_id( 118 ) == 2 );
        REQUIRE( dstrset.get_offset( 118 ) == 36 );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstrset, another_tmpfile );
        clear( dstrset );
        open( dstrset, another_tmpfile );

        THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
        {
          check_content( dstrset, raw_total );
          REQUIRE( dstrset.length() == 3 );
          REQUIRE( length( dstrset[ 0 ] ) == 38 );
          REQUIRE( length( dstrset[ 1 ] ) == 42 );
          REQUIRE( length( dstrset[ 2 ] ) == 37 );
          REQUIRE( dstrset.get_id( 0 ) == 0 );
          REQUIRE( dstrset.get_offset( 0 ) == 0 );
          REQUIRE( dstrset.get_id( 25 ) == 0 );
          REQUIRE( dstrset.get_offset( 25 ) == 25 );
          REQUIRE( dstrset.get_id( 37 ) == 0 );
          REQUIRE( dstrset.get_offset( 37 ) == 37 );
          REQUIRE( dstrset.get_id( 39 ) == 1 );
          REQUIRE( dstrset.get_offset( 39 ) == 0 );
          REQUIRE( dstrset.get_id( 51 ) == 1 );
          REQUIRE( dstrset.get_offset( 51 ) == 12 );
          REQUIRE( dstrset.get_id( 80 ) == 1 );
          REQUIRE( dstrset.get_offset( 80 ) == 41 );
          REQUIRE( dstrset.get_id( 82 ) == 2 );
          REQUIRE( dstrset.get_offset( 82 ) == 0 );
          REQUIRE( dstrset.get_id( 100 ) == 2 );
          REQUIRE( dstrset.get_offset( 100 ) == 18 );
          REQUIRE( dstrset.get_id( 118 ) == 2 );
          REQUIRE( dstrset.get_offset( 118 ) == 36 );
        }
      }
    }

    GIVEN( "A disk-based string set containing two of those strings" )
    {
      seqan::StringSet< DiskString > dstrset;
      dstrset.push_back( str1 );
      push_back( dstrset, str2 );

      WHEN( "It is moved to another string set" )
      {
        seqan::StringSet< DiskString > another_strset = std::move( dstrset );

        THEN( "It should contains the same strings" )
        {
          check_content( another_strset, raw_total.substr( 0, 81 ) );
          REQUIRE( another_strset.length() == 2 );
          REQUIRE( length( another_strset[ 0 ] ) == 38 );
          REQUIRE( length( another_strset[ 1 ] ) == 42 );
          REQUIRE( another_strset.get_id( 0 ) == 0 );
          REQUIRE( another_strset.get_offset( 0 ) == 0 );
          REQUIRE( another_strset.get_id( 25 ) == 0 );
          REQUIRE( another_strset.get_offset( 25 ) == 25 );
          REQUIRE( another_strset.get_id( 37 ) == 0 );
          REQUIRE( another_strset.get_offset( 37 ) == 37 );
          REQUIRE( another_strset.get_id( 39 ) == 1 );
          REQUIRE( another_strset.get_offset( 39 ) == 0 );
          REQUIRE( another_strset.get_id( 51 ) == 1 );
          REQUIRE( another_strset.get_offset( 51 ) == 12 );
          REQUIRE( another_strset.get_id( 80 ) == 1 );
          REQUIRE( another_strset.get_offset( 80 ) == 41 );
        }

        AND_WHEN( "Another string is added" )
        {
          another_strset.push_back( str3 );

          THEN( "It should contains all strings" )
          {
            check_content( another_strset, raw_total );
            REQUIRE( another_strset.length() == 3 );
            REQUIRE( length( another_strset[ 0 ] ) == 38 );
            REQUIRE( length( another_strset[ 1 ] ) == 42 );
            REQUIRE( length( another_strset[ 2 ] ) == 37 );
            REQUIRE( another_strset.get_id( 0 ) == 0 );
            REQUIRE( another_strset.get_offset( 0 ) == 0 );
            REQUIRE( another_strset.get_id( 25 ) == 0 );
            REQUIRE( another_strset.get_offset( 25 ) == 25 );
            REQUIRE( another_strset.get_id( 37 ) == 0 );
            REQUIRE( another_strset.get_offset( 37 ) == 37 );
            REQUIRE( another_strset.get_id( 39 ) == 1 );
            REQUIRE( another_strset.get_offset( 39 ) == 0 );
            REQUIRE( another_strset.get_id( 51 ) == 1 );
            REQUIRE( another_strset.get_offset( 51 ) == 12 );
            REQUIRE( another_strset.get_id( 80 ) == 1 );
            REQUIRE( another_strset.get_offset( 80 ) == 41 );
            REQUIRE( another_strset.get_id( 82 ) == 2 );
            REQUIRE( another_strset.get_offset( 82 ) == 0 );
            REQUIRE( another_strset.get_id( 100 ) == 2 );
            REQUIRE( another_strset.get_offset( 100 ) == 18 );
            REQUIRE( another_strset.get_id( 118 ) == 2 );
            REQUIRE( another_strset.get_offset( 118 ) == 36 );
          }

          AND_WHEN( "It is written to file and loaded again" )
          {
            std::string another_tmpfile = SEQAN_TEMP_FILENAME();
            save( another_strset, another_tmpfile );
            clear( another_strset );
            open( another_strset, another_tmpfile );

            THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
            {
              check_content( another_strset, raw_total );
              REQUIRE( another_strset.length() == 3 );
              REQUIRE( length( another_strset[ 0 ] ) == 38 );
              REQUIRE( length( another_strset[ 1 ] ) == 42 );
              REQUIRE( length( another_strset[ 2 ] ) == 37 );
              REQUIRE( another_strset.get_id( 0 ) == 0 );
              REQUIRE( another_strset.get_offset( 0 ) == 0 );
              REQUIRE( another_strset.get_id( 25 ) == 0 );
              REQUIRE( another_strset.get_offset( 25 ) == 25 );
              REQUIRE( another_strset.get_id( 37 ) == 0 );
              REQUIRE( another_strset.get_offset( 37 ) == 37 );
              REQUIRE( another_strset.get_id( 39 ) == 1 );
              REQUIRE( another_strset.get_offset( 39 ) == 0 );
              REQUIRE( another_strset.get_id( 51 ) == 1 );
              REQUIRE( another_strset.get_offset( 51 ) == 12 );
              REQUIRE( another_strset.get_id( 80 ) == 1 );
              REQUIRE( another_strset.get_offset( 80 ) == 41 );
              REQUIRE( another_strset.get_id( 82 ) == 2 );
              REQUIRE( another_strset.get_offset( 82 ) == 0 );
              REQUIRE( another_strset.get_id( 100 ) == 2 );
              REQUIRE( another_strset.get_offset( 100 ) == 18 );
              REQUIRE( another_strset.get_id( 118 ) == 2 );
              REQUIRE( another_strset.get_offset( 118 ) == 36 );
            }
          }
        }
      }

      WHEN( "It is cleared" )
      {
        clear( dstrset );

        THEN( "Its length should be zero" )
        {
          REQUIRE( length( dstrset ) == 0 );
        }

        AND_WHEN( "A string is added after clearance" )
        {
          push_back( dstrset, str3 );

          THEN( "It should contain that string" )
          {
            REQUIRE( dstrset.length() == 1 );
            REQUIRE( length( dstrset[ 0 ] ) == 37 );
            check_content( dstrset, str3 );
            REQUIRE( dstrset.get_id( 0 ) == 0 );
            REQUIRE( dstrset.get_offset( 0 ) == 0 );
            REQUIRE( dstrset.get_id( 25 ) == 0 );
            REQUIRE( dstrset.get_offset( 25 ) == 25 );
            REQUIRE( dstrset.get_id( 36 ) == 0 );
            REQUIRE( dstrset.get_offset( 36 ) == 36 );
          }
        }
      }
    }
  }
}

SCENARIO( "Constructing a set of in-memory string", "[sequence]" )
{
  auto check_content =
    []( seqan::StringSet< MemString >& d, const std::string& data ) {
      REQUIRE( ( data == d ) );
    };

  GIVEN( "A set of strings" )
  {
    std::string str1 = "a mississippian lazy fox sits on a pie";
    std::string str2 = "another brazilian cute beaver builds a dam";
    std::string str3 = "some african stupid chimps eat banana";
    std::string raw_total =
      str1 + std::string( 1, SEQUENCE_DEFAULT_SENTINEL_CHAR ) +
      str2 + std::string( 1, SEQUENCE_DEFAULT_SENTINEL_CHAR ) + str3;

    WHEN( "The strings are added to string set by const reference" )
    {
      seqan::StringSet< MemString, seqan::Owner<> > dstrset;
      dstrset.reserve( 3 );
      appendValue( dstrset, str1 );
      push_back( dstrset, str2 );
      dstrset.push_back( str3 );

      THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
      {
        check_content( dstrset, raw_total );
        REQUIRE( dstrset.length() == 3 );
        REQUIRE( length( dstrset[ 0 ] ) == 38 );
        REQUIRE( length( dstrset[ 1 ] ) == 42 );
        REQUIRE( length( dstrset[ 2 ] ) == 37 );
        REQUIRE( dstrset.get_id( 0 ) == 0 );
        REQUIRE( dstrset.get_offset( 0 ) == 0 );
        REQUIRE( dstrset.get_id( 25 ) == 0 );
        REQUIRE( dstrset.get_offset( 25 ) == 25 );
        REQUIRE( dstrset.get_id( 37 ) == 0 );
        REQUIRE( dstrset.get_offset( 37 ) == 37 );
        REQUIRE( dstrset.get_id( 39 ) == 1 );
        REQUIRE( dstrset.get_offset( 39 ) == 0 );
        REQUIRE( dstrset.get_id( 51 ) == 1 );
        REQUIRE( dstrset.get_offset( 51 ) == 12 );
        REQUIRE( dstrset.get_id( 80 ) == 1 );
        REQUIRE( dstrset.get_offset( 80 ) == 41 );
        REQUIRE( dstrset.get_id( 82 ) == 2 );
        REQUIRE( dstrset.get_offset( 82 ) == 0 );
        REQUIRE( dstrset.get_id( 100 ) == 2 );
        REQUIRE( dstrset.get_offset( 100 ) == 18 );
        REQUIRE( dstrset.get_id( 118 ) == 2 );
        REQUIRE( dstrset.get_offset( 118 ) == 36 );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstrset, another_tmpfile );
        clear( dstrset );
        open( dstrset, another_tmpfile );

        THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
        {
          check_content( dstrset, raw_total );
          REQUIRE( dstrset.length() == 3 );
          REQUIRE( length( dstrset[ 0 ] ) == 38 );
          REQUIRE( length( dstrset[ 1 ] ) == 42 );
          REQUIRE( length( dstrset[ 2 ] ) == 37 );
          REQUIRE( dstrset.get_id( 0 ) == 0 );
          REQUIRE( dstrset.get_offset( 0 ) == 0 );
          REQUIRE( dstrset.get_id( 25 ) == 0 );
          REQUIRE( dstrset.get_offset( 25 ) == 25 );
          REQUIRE( dstrset.get_id( 37 ) == 0 );
          REQUIRE( dstrset.get_offset( 37 ) == 37 );
          REQUIRE( dstrset.get_id( 39 ) == 1 );
          REQUIRE( dstrset.get_offset( 39 ) == 0 );
          REQUIRE( dstrset.get_id( 51 ) == 1 );
          REQUIRE( dstrset.get_offset( 51 ) == 12 );
          REQUIRE( dstrset.get_id( 80 ) == 1 );
          REQUIRE( dstrset.get_offset( 80 ) == 41 );
          REQUIRE( dstrset.get_id( 82 ) == 2 );
          REQUIRE( dstrset.get_offset( 82 ) == 0 );
          REQUIRE( dstrset.get_id( 100 ) == 2 );
          REQUIRE( dstrset.get_offset( 100 ) == 18 );
          REQUIRE( dstrset.get_id( 118 ) == 2 );
          REQUIRE( dstrset.get_offset( 118 ) == 36 );
        }
      }
    }

    WHEN( "The strings are added to string set by r-value reference" )
    {
      seqan::StringSet< MemString, seqan::Owner<> > dstrset;
      reserve( dstrset, 3 );
      appendValue( dstrset, std::move( str1 ) );
      push_back( dstrset, std::move( str2 ) );
      dstrset.push_back( std::move( str3 ) );

      THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
      {
        check_content( dstrset, raw_total );
        REQUIRE( length( dstrset ) == 3 );
        REQUIRE( length( dstrset[ 0 ] ) == 38 );
        REQUIRE( length( dstrset[ 1 ] ) == 42 );
        REQUIRE( length( dstrset[ 2 ] ) == 37 );
        REQUIRE( dstrset.get_id( 0 ) == 0 );
        REQUIRE( dstrset.get_offset( 0 ) == 0 );
        REQUIRE( dstrset.get_id( 25 ) == 0 );
        REQUIRE( dstrset.get_offset( 25 ) == 25 );
        REQUIRE( dstrset.get_id( 37 ) == 0 );
        REQUIRE( dstrset.get_offset( 37 ) == 37 );
        REQUIRE( dstrset.get_id( 39 ) == 1 );
        REQUIRE( dstrset.get_offset( 39 ) == 0 );
        REQUIRE( dstrset.get_id( 51 ) == 1 );
        REQUIRE( dstrset.get_offset( 51 ) == 12 );
        REQUIRE( dstrset.get_id( 80 ) == 1 );
        REQUIRE( dstrset.get_offset( 80 ) == 41 );
        REQUIRE( dstrset.get_id( 82 ) == 2 );
        REQUIRE( dstrset.get_offset( 82 ) == 0 );
        REQUIRE( dstrset.get_id( 100 ) == 2 );
        REQUIRE( dstrset.get_offset( 100 ) == 18 );
        REQUIRE( dstrset.get_id( 118 ) == 2 );
        REQUIRE( dstrset.get_offset( 118 ) == 36 );
      }

      AND_WHEN( "It is written to file and loaded again" )
      {
        std::string another_tmpfile = SEQAN_TEMP_FILENAME();
        save( dstrset, another_tmpfile );
        clear( dstrset );
        open( dstrset, another_tmpfile );

        THEN( "It should contain the concatenation of all strings delimited by SENTINEL" )
        {
          check_content( dstrset, raw_total );
          REQUIRE( dstrset.length() == 3 );
          REQUIRE( length( dstrset[ 0 ] ) == 38 );
          REQUIRE( length( dstrset[ 1 ] ) == 42 );
          REQUIRE( length( dstrset[ 2 ] ) == 37 );
          REQUIRE( dstrset.get_id( 0 ) == 0 );
          REQUIRE( dstrset.get_offset( 0 ) == 0 );
          REQUIRE( dstrset.get_id( 25 ) == 0 );
          REQUIRE( dstrset.get_offset( 25 ) == 25 );
          REQUIRE( dstrset.get_id( 37 ) == 0 );
          REQUIRE( dstrset.get_offset( 37 ) == 37 );
          REQUIRE( dstrset.get_id( 39 ) == 1 );
          REQUIRE( dstrset.get_offset( 39 ) == 0 );
          REQUIRE( dstrset.get_id( 51 ) == 1 );
          REQUIRE( dstrset.get_offset( 51 ) == 12 );
          REQUIRE( dstrset.get_id( 80 ) == 1 );
          REQUIRE( dstrset.get_offset( 80 ) == 41 );
          REQUIRE( dstrset.get_id( 82 ) == 2 );
          REQUIRE( dstrset.get_offset( 82 ) == 0 );
          REQUIRE( dstrset.get_id( 100 ) == 2 );
          REQUIRE( dstrset.get_offset( 100 ) == 18 );
          REQUIRE( dstrset.get_id( 118 ) == 2 );
          REQUIRE( dstrset.get_offset( 118 ) == 36 );
        }
      }
    }

    GIVEN( "A in-memory string set containing two of those strings" )
    {
      seqan::StringSet< MemString > dstrset;
      dstrset.push_back( str1 );
      push_back( dstrset, str2 );

      WHEN( "It is moved to another string set" )
      {
        seqan::StringSet< MemString > another_strset = std::move( dstrset );

        THEN( "It should contains the same strings" )
        {
          check_content( another_strset, raw_total.substr( 0, 81 ) );
          REQUIRE( another_strset.length() == 2 );
          REQUIRE( length( another_strset[ 0 ] ) == 38 );
          REQUIRE( length( another_strset[ 1 ] ) == 42 );
          REQUIRE( another_strset.get_id( 0 ) == 0 );
          REQUIRE( another_strset.get_offset( 0 ) == 0 );
          REQUIRE( another_strset.get_id( 25 ) == 0 );
          REQUIRE( another_strset.get_offset( 25 ) == 25 );
          REQUIRE( another_strset.get_id( 37 ) == 0 );
          REQUIRE( another_strset.get_offset( 37 ) == 37 );
          REQUIRE( another_strset.get_id( 39 ) == 1 );
          REQUIRE( another_strset.get_offset( 39 ) == 0 );
          REQUIRE( another_strset.get_id( 51 ) == 1 );
          REQUIRE( another_strset.get_offset( 51 ) == 12 );
          REQUIRE( another_strset.get_id( 80 ) == 1 );
          REQUIRE( another_strset.get_offset( 80 ) == 41 );
        }

        AND_WHEN( "Another string is added" )
        {
          another_strset.push_back( str3 );

          THEN( "It should contains all strings" )
          {
            check_content( another_strset, raw_total );
            REQUIRE( another_strset.length() == 3 );
            REQUIRE( length( another_strset[ 0 ] ) == 38 );
            REQUIRE( length( another_strset[ 1 ] ) == 42 );
            REQUIRE( length( another_strset[ 2 ] ) == 37 );
            REQUIRE( another_strset.get_id( 0 ) == 0 );
            REQUIRE( another_strset.get_offset( 0 ) == 0 );
            REQUIRE( another_strset.get_id( 25 ) == 0 );
            REQUIRE( another_strset.get_offset( 25 ) == 25 );
            REQUIRE( another_strset.get_id( 37 ) == 0 );
            REQUIRE( another_strset.get_offset( 37 ) == 37 );
            REQUIRE( another_strset.get_id( 39 ) == 1 );
            REQUIRE( another_strset.get_offset( 39 ) == 0 );
            REQUIRE( another_strset.get_id( 51 ) == 1 );
            REQUIRE( another_strset.get_offset( 51 ) == 12 );
            REQUIRE( another_strset.get_id( 80 ) == 1 );
            REQUIRE( another_strset.get_offset( 80 ) == 41 );
            REQUIRE( another_strset.get_id( 82 ) == 2 );
            REQUIRE( another_strset.get_offset( 82 ) == 0 );
            REQUIRE( another_strset.get_id( 100 ) == 2 );
            REQUIRE( another_strset.get_offset( 100 ) == 18 );
            REQUIRE( another_strset.get_id( 118 ) == 2 );
            REQUIRE( another_strset.get_offset( 118 ) == 36 );
          }

          AND_WHEN( "It is written to file and loaded again" )
          {
            std::string another_tmpfile = SEQAN_TEMP_FILENAME();
            save( another_strset, another_tmpfile );
            clear( another_strset );
            open( another_strset, another_tmpfile );

            THEN( "It should contains all strings" )
            {
              check_content( another_strset, raw_total );
              REQUIRE( another_strset.length() == 3 );
              REQUIRE( length( another_strset[ 0 ] ) == 38 );
              REQUIRE( length( another_strset[ 1 ] ) == 42 );
              REQUIRE( length( another_strset[ 2 ] ) == 37 );
              REQUIRE( another_strset.get_id( 0 ) == 0 );
              REQUIRE( another_strset.get_offset( 0 ) == 0 );
              REQUIRE( another_strset.get_id( 25 ) == 0 );
              REQUIRE( another_strset.get_offset( 25 ) == 25 );
              REQUIRE( another_strset.get_id( 37 ) == 0 );
              REQUIRE( another_strset.get_offset( 37 ) == 37 );
              REQUIRE( another_strset.get_id( 39 ) == 1 );
              REQUIRE( another_strset.get_offset( 39 ) == 0 );
              REQUIRE( another_strset.get_id( 51 ) == 1 );
              REQUIRE( another_strset.get_offset( 51 ) == 12 );
              REQUIRE( another_strset.get_id( 80 ) == 1 );
              REQUIRE( another_strset.get_offset( 80 ) == 41 );
              REQUIRE( another_strset.get_id( 82 ) == 2 );
              REQUIRE( another_strset.get_offset( 82 ) == 0 );
              REQUIRE( another_strset.get_id( 100 ) == 2 );
              REQUIRE( another_strset.get_offset( 100 ) == 18 );
              REQUIRE( another_strset.get_id( 118 ) == 2 );
              REQUIRE( another_strset.get_offset( 118 ) == 36 );
            }
          }
        }
      }

      WHEN( "It is cleared" )
      {
        clear( dstrset );

        THEN( "Its length should be zero" )
        {
          REQUIRE( length( dstrset ) == 0 );
        }

        AND_WHEN( "A string is added after clearance" )
        {
          push_back( dstrset, str3 );

          THEN( "It should contain that string" )
          {
            REQUIRE( dstrset.length() == 1 );
            check_content( dstrset, str3 );
            REQUIRE( length( dstrset[ 0 ] ) == 37 );
            REQUIRE( dstrset.get_id( 0 ) == 0 );
            REQUIRE( dstrset.get_offset( 0 ) == 0 );
            REQUIRE( dstrset.get_id( 25 ) == 0 );
            REQUIRE( dstrset.get_offset( 25 ) == 25 );
            REQUIRE( dstrset.get_id( 36 ) == 0 );
            REQUIRE( dstrset.get_offset( 36 ) == 36 );
          }
        }
      }
    }
  }
}

SCENARIO( "Constructing a mutable YaPair", "[sequence]" )
{
  GIVEN( "A YaPair object" )
  {
    YaPair< int, int > a( 3, 4 );

    WHEN( "The `first` and `second` member are modified" )
    {
      a.first = 2;
      a.second = 5;

      THEN( "The `i1` and `i2` should get the updated value" )
      {
        REQUIRE( a.i1 == 2 );
        REQUIRE( a.i2 == 5 );
        REQUIRE( a.first == a.i1 );
        REQUIRE( a.second == a.i2 );
      }
    }

    WHEN( "The `i1` and `i2` member are modified" )
    {
      a.i1 = 2;
      a.i2 = 5;

      THEN( "The `i1` and `i2` should get the updated value" )
      {
        REQUIRE( a.first == 2 );
        REQUIRE( a.second == 5 );
        REQUIRE( a.first == a.i1 );
        REQUIRE( a.second == a.i2 );
      }
    }

    WHEN( "It is assigned to another YaPair object" )
    {
      YaPair< int, int > b;
      REQUIRE( b.first == 0 );
      REQUIRE( b.second == 0 );
      REQUIRE( b.i1 == 0 );
      REQUIRE( b.i2 == 0 );

      b = a;

      THEN( "It should get the assigned value" )
      {
        REQUIRE( b.i1 == 3 );
        REQUIRE( b.i2 == 4 );
        REQUIRE( b.first == b.i1 );
        REQUIRE( b.second == b.i2 );
      }
    }

    WHEN( "It is assigned to another moved YaPair object" )
    {
      YaPair< int, int > b;
      REQUIRE( b.first == 0 );
      REQUIRE( b.second == 0 );
      REQUIRE( b.i1 == 0 );
      REQUIRE( b.i2 == 0 );

      b = std::move( a );

      THEN( "It should get the assigned value" )
      {
        REQUIRE( b.i1 == 3 );
        REQUIRE( b.i2 == 4 );
        REQUIRE( b.first == b.i1 );
        REQUIRE( b.second == b.i2 );
      }
    }

    WHEN( "A YaPair object constructed by another YaPair object" )
    {
      YaPair< int, int > b( a );

      THEN( "It should get the assigned value" )
      {
        REQUIRE( b.i1 == 3 );
        REQUIRE( b.i2 == 4 );
        REQUIRE( b.first == b.i1 );
        REQUIRE( b.second == b.i2 );
      }
    }

    WHEN( "A YaPair object constructed by another moved YaPair object" )
    {
      YaPair< int, int > b( std::move( a ) );

      THEN( "It should get the assigned value" )
      {
        REQUIRE( b.i1 == 3 );
        REQUIRE( b.i2 == 4 );
        REQUIRE( b.first == b.i1 );
        REQUIRE( b.second == b.i2 );
      }
    }
  }
}

SCENARIO( "Enumerate k-mers in a Records using RecordsIter class", "[sequence]" )
{
  GIVEN( "A records containing four sequences with different lengths" )
  {
    Records< Dna5QStringSet<> > reads;
    appendValue( reads.str, "aaaaaattttttcccccc" );
    appendValue( reads.str, "acgtttacgtttacg" );
    appendValue( reads.str, "acgtttacgtttacgtttacgttt" );
    appendValue( reads.str, "acgtttacgtttacgtttacgtttaaaaaattttttc" );

    unsigned int k = 6;
    WHEN( "Enumerating non-overlapping " + std::to_string( k ) + "-mers in the records" )
    {
      typename seqan::Iterator< decltype( reads ), NonOverlapping >::Type iter( &reads, k );

      THEN( "It should yield all non-overlapping " + std::to_string( k ) + "-mers in the records" )
      {
        REQUIRE( *iter++ == "aaaaaa" );
        REQUIRE( *iter++ == "tttttt" );
        REQUIRE( *iter++ == "cccccc" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "aaaaaa" );
        REQUIRE( *iter++ == "tttttt" );
        REQUIRE( at_end( iter ) == true );
      }
    }
  }

  GIVEN( "A records containing four sequences with same lengths" )
  {
    Records< Dna5QStringSet<> > reads;
    appendValue( reads.str, "aaaaaattttttcccccc" );
    appendValue( reads.str, "acgtttacgtttacg" );

    unsigned int k = 6;
    WHEN( "Enumerating overlapping " + std::to_string( k ) + "-mers in the records" )
    {
      typename seqan::Iterator< decltype( reads ), GreedyOverlapping >::Type iter( &reads, k );

      THEN( "It should yield all overlapping " + std::to_string( k ) + "-mers in the records" )
      {
        REQUIRE( *iter++ == "aaaaaa" );
        REQUIRE( *iter++ == "aaaaat" );
        REQUIRE( *iter++ == "aaaatt" );
        REQUIRE( *iter++ == "aaattt" );
        REQUIRE( *iter++ == "aatttt" );
        REQUIRE( *iter++ == "attttt" );
        REQUIRE( *iter++ == "tttttt" );
        REQUIRE( *iter++ == "tttttc" );
        REQUIRE( *iter++ == "ttttcc" );
        REQUIRE( *iter++ == "tttccc" );
        REQUIRE( *iter++ == "ttcccc" );
        REQUIRE( *iter++ == "tccccc" );
        REQUIRE( *iter++ == "cccccc" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "cgttta" );
        REQUIRE( *iter++ == "gtttac" );
        REQUIRE( *iter++ == "tttacg" );
        REQUIRE( *iter++ == "ttacgt" );
        REQUIRE( *iter++ == "tacgtt" );
        REQUIRE( *iter++ == "acgttt" );
        REQUIRE( *iter++ == "cgttta" );
        REQUIRE( *iter++ == "gtttac" );
        REQUIRE( *iter++ == "tttacg" );
        REQUIRE( at_end( iter ) == true );
      }
    }
  }

  unsigned int reads_num = 10;
  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = test_data_dir + "/small/reads_n10l10e0i0.fastq";
    klibpp::SeqStreamIn iss( fqpath.c_str() );
    if ( !iss ) throw std::runtime_error( "cannot open file '" + fqpath + "'" );

    typedef seqan::StringSet< seqan::DnaQString > TStringSet;
    typedef Records< TStringSet > TRecords;

    TRecords reads;
    readRecords( reads, iss );
    unsigned int k = 4;

    unsigned int step = 2;
    WHEN( "Seeding by non-greedy overlapping strategy with length "
        + std::to_string( k ) + " and step size " + std::to_string( step ) )
    {
      typedef typename seqan::Iterator< TRecords, Overlapping >::Type TIter;
      TIter reads_itr( &reads, k, step );

      THEN( "Seed should be correct" )
      {
        std::string truth[70] =
        { "CAAA", "AATA", "TAAG", "AGAT",
          "AAAT", "ATAA", "AAGA", "GACT",
          "TTTC", "TCTG", "TGGA", "GAGT",
          "ATAA", "AATA", "TATT", "TTCC",
          "TTCC", "CCTG", "TGGT", "GTTG",
          "GTCC", "CCTG", "TGGT", "GTTG",
          "TGCT", "CTAT", "ATGT", "GTGT",
          "TGTT", "TTGG", "GGGC", "GCTT",
          "CTTT", "TTTT", "TTTC", "TCTT",
          "CTTC", "TCTT", "TTCC", "CCTT" };
        unsigned int i = 0;
        while ( !at_end( reads_itr ) ) {
          REQUIRE( *reads_itr == truth[i] );
          REQUIRE( get_position( reads_itr ).i1 == i/4 );
          REQUIRE( get_position( reads_itr ).i2 == i%4*step );
          ++i;
          ++reads_itr;
        }
      }
    }
  }
}

SCENARIO( "Increment a k-mer lexicographically", "[sequence]" )
{
  unsigned int k = 20;
  GIVEN( "A k-mer of length " + std::to_string( k ) + " with all 'A'" )
  {
    seqan::DnaString kmer;
    for ( unsigned int i = 0; i < k; ++i ) appendValue( kmer, 'A' );

    WHEN( "It is incremented" )
    {
      unsigned int s = increment_kmer( kmer );
      REQUIRE( s == length( kmer ) - 1 );
      THEN( "It should be the next lexicographical kmer" )
      {
        REQUIRE( kmer == "AAAAAAAAAAAAAAAAAAAC" );
      }
    }

    WHEN( "It is incremented at two positions in the middle of the string" )
    {
      unsigned int s = increment_kmer( kmer, 12 );
      REQUIRE( s == 11 );
      s = increment_kmer( kmer, 17 );
      REQUIRE( s == 16 );
      THEN( "It should be the next lexicographical kmer at those positions" )
      {
        REQUIRE( kmer == "AAAAAAAAAAACAAAACAAA" );
      }
    }

    WHEN( "It is incremented at a position out of range" )
    {
      unsigned int s = increment_kmer( kmer, 32 );
      REQUIRE( s == length( kmer ) - 1 );
      THEN( "It should be the next lexicographical kmer at last character" )
      {
        REQUIRE( kmer == "AAAAAAAAAAAAAAAAAAAC" );
      }
    }

    WHEN( "It is incremented at a position out of range" )
    {
      unsigned int s = increment_kmer( kmer, -1 );
      REQUIRE( s == length( kmer ) - 1 );
      THEN( "It should be the next lexicographical kmer at last character" )
      {
        REQUIRE( kmer == "AAAAAAAAAAAAAAAAAAAC" );
      }
    }
  }

  GIVEN( "A k-mer of length " + std::to_string( k ) + " with all 'T'" )
  {
    seqan::DnaString kmer;
    for ( unsigned int i = 0; i < k; ++i ) appendValue( kmer, 'T' );

    WHEN( "It is incremented" )
    {
      unsigned int s = increment_kmer( kmer );
      REQUIRE( s == -1 );
      THEN( "It should not be changed" )
      {
        REQUIRE( kmer == "TTTTTTTTTTTTTTTTTTTT" );
      }
    }
  }
}

SCENARIO( "Seeding", "[seeding][sequence]" )
{
  unsigned int reads_num = 10;
  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = test_data_dir + "/small/reads_n10l10e0i0.fastq";
    seqan::SeqFileIn infile;
    if ( !open( infile, fqpath.c_str() ) ) {
      throw std::runtime_error( "cannot open file '" + fqpath + "'" );
    }
    Records< Dna5QStringSet<> > reads;
    readRecords( reads, infile );
    unsigned int k = 4;

    WHEN( "Seeding by greedy non-overlapping strategy with length " + std::to_string( k ) )
    {
      Dna5QStringSet<> seeds;
      seeding( seeds, reads.str, k, GreedyNonOverlapping() );

      THEN( "Seed should be correct" )
      {
        std::string truth[30] =
        { "CAAA", "TAAG", "AGAT",
          "AAAT", "AAGA", "GACT",
          "TTTC", "TGGA", "GAGT",
          "ATAA", "TATT", "TTCC",
          "TTCC", "TGGT", "GTTG",
          "GTCC", "TGGT", "GTTG",
          "TGCT", "ATGT", "GTGT",
          "TGTT", "GGGC", "GCTT",
          "CTTT", "TTTC", "TCTT",
          "CTTC", "TTCC", "CCTT" };
        for ( unsigned int i = 0; i < length( seeds ); ++i ) {
          REQUIRE( seeds[i] == truth[i] );
        }
      }
    }

    WHEN( "Seeding by non-overlapping strategy with length " + std::to_string( k ) )
    {
      Dna5QStringSet<> seeds;
      seeding( seeds, reads.str, k, NonOverlapping() );

      THEN( "Seed should be correct" )
      {
        std::string truth[20] =
        { "CAAA", "TAAG",
          "AAAT", "AAGA",
          "TTTC", "TGGA",
          "ATAA", "TATT",
          "TTCC", "TGGT",
          "GTCC", "TGGT",
          "TGCT", "ATGT",
          "TGTT", "GGGC",
          "CTTT", "TTTC",
          "CTTC", "TTCC" };
        for ( unsigned int i = 0; i < length( seeds ); ++i ) {
          REQUIRE( seeds[i] == truth[i] );
        }
      }
    }

    WHEN( "Seeding by overlapping strategy with length " + std::to_string( k ) )
    {
      Dna5QStringSet<> seeds;
      seeding( seeds, reads.str, k, GreedyOverlapping() );

      THEN( "Seed should be correct" )
      {
        std::string truth[70] =
        { "CAAA", "AAAT", "AATA", "ATAA", "TAAG", "AAGA", "AGAT",
          "AAAT", "AATA", "ATAA", "TAAG", "AAGA", "AGAC", "GACT",
          "TTTC", "TTCT", "TCTG", "CTGG", "TGGA", "GGAG", "GAGT",
          "ATAA", "TAAT", "AATA", "ATAT", "TATT", "ATTC", "TTCC",
          "TTCC", "TCCT", "CCTG", "CTGG", "TGGT", "GGTT", "GTTG",
          "GTCC", "TCCT", "CCTG", "CTGG", "TGGT", "GGTT", "GTTG",
          "TGCT", "GCTA", "CTAT", "TATG", "ATGT", "TGTG", "GTGT",
          "TGTT", "GTTG", "TTGG", "TGGG", "GGGC", "GGCT", "GCTT",
          "CTTT", "TTTT", "TTTT", "TTTT", "TTTC", "TTCT", "TCTT",
          "CTTC", "TTCT", "TCTT", "CTTC", "TTCC", "TCCT", "CCTT" };
        for ( unsigned int i = 0; i < length( seeds ); ++i ) {
          REQUIRE( seeds[i] == truth[i] );
        }
      }
    }
  }

  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = test_data_dir + "/small/reads_n10l10e0i0.fastq";
    klibpp::SeqStreamIn iss( fqpath.c_str() );
    if ( !iss ) throw std::runtime_error( "cannot open file '" + fqpath + "'" );

    typedef seqan::StringSet< std::string > TStringSet;

    Records< TStringSet > reads;
    readRecords( reads, iss );
    unsigned int k = 4;

    WHEN( "Seeding by non-overlapping strategy with length " + std::to_string( k ) )
    {
      Records< TStringSet > seeds;
      seeding( seeds, reads, k, NonOverlapping() );

      THEN( "Seed should be correct" )
      {
        std::string truth[20] =
        { "CAAA", "TAAG",
          "AAAT", "AAGA",
          "TTTC", "TGGA",
          "ATAA", "TATT",
          "TTCC", "TGGT",
          "GTCC", "TGGT",
          "TGCT", "ATGT",
          "TGTT", "GGGC",
          "CTTT", "TTTC",
          "CTTC", "TTCC" };
        for ( unsigned int i = 0; i < length( seeds ); ++i ) {
          REQUIRE( seeds.str[i] == truth[i] );
          for ( unsigned int j = 0; j < k; ++j ) {
            REQUIRE( position_to_id( seeds, { i, j } ) == i/2 );
            REQUIRE( position_to_offset( seeds, { i, j } ) == (i%2)*k+j );
          }
        }
      }
    }

    WHEN( "Seeding by overlapping strategy with length " + std::to_string( k ) )
    {
      Records< TStringSet > seeds;
      seeding( seeds, reads, k, GreedyOverlapping() );

      THEN( "Seed should be correct" )
      {
        std::string truth[70] =
        { "CAAA", "AAAT", "AATA", "ATAA", "TAAG", "AAGA", "AGAT",
          "AAAT", "AATA", "ATAA", "TAAG", "AAGA", "AGAC", "GACT",
          "TTTC", "TTCT", "TCTG", "CTGG", "TGGA", "GGAG", "GAGT",
          "ATAA", "TAAT", "AATA", "ATAT", "TATT", "ATTC", "TTCC",
          "TTCC", "TCCT", "CCTG", "CTGG", "TGGT", "GGTT", "GTTG",
          "GTCC", "TCCT", "CCTG", "CTGG", "TGGT", "GGTT", "GTTG",
          "TGCT", "GCTA", "CTAT", "TATG", "ATGT", "TGTG", "GTGT",
          "TGTT", "GTTG", "TTGG", "TGGG", "GGGC", "GGCT", "GCTT",
          "CTTT", "TTTT", "TTTT", "TTTT", "TTTC", "TTCT", "TCTT",
          "CTTC", "TTCT", "TCTT", "CTTC", "TTCC", "TCCT", "CCTT" };
        for ( unsigned int i = 0; i < length( seeds ); ++i ) {
          REQUIRE( seeds.str[i] == truth[i] );
          for ( unsigned int j = 0; j < k; ++j ) {
            REQUIRE( position_to_id( seeds, { i, j } ) == i/7 );
            REQUIRE( position_to_offset( seeds, { i, j } ) == (i%7)+j );
          }
        }
      }
    }
  }
}
