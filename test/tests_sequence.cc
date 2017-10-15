/**
 *    @file  tests_sequence.cc
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

#include "tests_base.h"
#include "sequence.h"
#include "logger.h"


using namespace grem;

SCENARIO( "Subsetting a reads chunk from a reads set", "[sequence]" )
{
  unsigned int reads_num = 10;
  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = _testdir + "/data/small/reads_n10l10e0i0.fastq";
    seqan::SeqFileIn infile;
    if ( !open( infile, fqpath.c_str() ) ) {
      throw std::runtime_error( "cannot open file '" + fqpath + "'" );
    }
    Records< Dna5QStringSet<> > reads;
    readRecords( reads, infile );
    Records< Dna5QStringSet< grem::Dependent > > reads_chunk;
    unsigned int subset_len = 4;

    WHEN( "A " + std::to_string( subset_len ) + "-reads subset is selected" )
    {
      load_chunk( reads_chunk, reads, subset_len, 2 );
      THEN( "Each read's ID in the chunk should be its position in the original set" )
      {
        REQUIRE( length( reads_chunk ) == subset_len );
        for ( unsigned int i = 0; i < subset_len; ++i ) {
          auto read_id = position_to_id( reads_chunk, i );
          REQUIRE( reads_chunk.str[ i ] == reads.str[ read_id ] );
        }
      }

      WHEN( "Next " + std::to_string( subset_len ) + "-reads subset is selected" )
      {
        load_chunk( reads_chunk, reads, subset_len );
        THEN( "Each read's ID in the chunk should be its position in the original set" )
        {
          REQUIRE( length( reads_chunk ) == subset_len );
          for ( unsigned int i = 0; i < subset_len; ++i ) {
            auto read_id = position_to_id( reads_chunk, i );
            REQUIRE( reads_chunk.str[ i ] == reads.str[ read_id ] );
          }
        }
      }

      subset_len = 8;
      WHEN( std::to_string( subset_len ) + "-reads subset is selected after clear call" )
      {
        clear( reads_chunk );
        load_chunk( reads_chunk, reads, subset_len );
        THEN( "Each read's ID in the chunk should be its position in the original set" )
        {
          REQUIRE( length( reads_chunk ) == subset_len );
          for ( unsigned int i = 0; i < subset_len; ++i ) {
            auto read_id = position_to_id( reads_chunk, i );
            REQUIRE( reads_chunk.str[ i ] == reads.str[ read_id ] );
          }
        }

        WHEN( "Next " + std::to_string( subset_len ) + "-reads subset is selected" )
        {
          load_chunk( reads_chunk, reads, subset_len );
          THEN( "Each read's ID in the chunk should be its position in the original set" )
          {
            REQUIRE( length( reads_chunk ) == reads_num - subset_len );
            for ( unsigned int i = 0; i < reads_num - subset_len; ++i ) {
              auto read_id = position_to_id( reads_chunk, i );
              REQUIRE( reads_chunk.str[ i ] == reads.str[ read_id ] );
            }
          }
        }
      }
    }
  }
}

SCENARIO( "Seeding", "[seeding][sequence]" )
{
  unsigned int reads_num = 10;
  GIVEN( "Read records from a file containing " + std::to_string( reads_num ) + " reads" )
  {
    std::string fqpath = _testdir + "/data/small/reads_n10l10e0i0.fastq";
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
}
