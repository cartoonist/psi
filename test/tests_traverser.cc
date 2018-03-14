/**
 *    @file  tests_traverser.cc
 *   @brief  Test traverser module.
 *
 *  Test scenarios for traverser module.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Mar 17, 2017  00:44
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <vector>
#include <string>

#include <seqan/seq_io.h>

#include "tests_base.h"
#include "sequence.h"
#include "index.h"
#include "vargraph.h"
#include "traverser.h"
#include "utils.h"
#include "logger.h"


using namespace grem;

SCENARIO ( "Find reads in the graph using a Traverser (exact)", "[traverser]" )
{
  GIVEN ( "A small variation graph and a set of reads" )
  {
    typedef seqan::IndexWotd<> TIndexSpec;
    typedef seqan::Index< Dna5QStringSet<>, TIndexSpec > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;

    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );
    std::string readspath = _testdir + "/data/small/reads_n10l10e0i0.fastq";
    seqan::SeqFileIn reads_file;
    if ( !open( reads_file, readspath.c_str() ) ) {
      throw std::runtime_error( "cannot open file " + readspath );
    }

    Records< Dna5QStringSet<> > reads;
    readRecords( reads, reads_file, 10 );
    seqan::Index< Dna5QStringSet<>, TIndexSpec > reads_index( reads.str );

    unsigned int seed_len = 10;

    WHEN ( "Run a traverser on all loci with seed length " + std::to_string( seed_len ) )
    {
      TTraverser traverser( &vargraph, &reads, &reads_index, seed_len );

      unsigned int counter = 0;
      std::size_t truth[10][2] = { {1, 0}, {1, 1}, {9, 4}, {9, 17}, {16, 0}, {17, 0},
        {20, 0}, {20, 31}, {20, 38}, {20, 38} };

      std::function< void( typename TTraverser::output_type const& ) > count_hits =
        [&counter, &truth]( typename TTraverser::output_type const& hit ) {
          if ( counter < 10 ) {
            REQUIRE( hit.node_id == truth[ counter ][0] );
            REQUIRE( hit.node_offset == truth[ counter ][1] );
            REQUIRE( hit.read_id == counter );
            REQUIRE( hit.read_offset == 0 );
          }
          else {
            assert( false );  // shouldn't be reached.
          }
          ++counter;
      };
      THEN ( "It should find all reads in the graph" )
      {
        for ( std::size_t r = 1; r <= vargraph.max_node_rank(); ++r ) {
          const auto& node_id = vargraph.rank_to_id( r );
          VarGraph::offset_type seqlen = vargraph.node_length( node_id );
          for ( VarGraph::offset_type f = 0; f < seqlen; ++f ) {
            traverser.set_start_locus( node_id, f );
            traverser.run( count_hits );
          }
        }
      }
    }
  }
}
