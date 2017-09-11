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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/seeds.h>

#include "tests_base.h"
#include "sequence.h"
#include "traverser.h"
#include "mapper.h"
#include "vargraph.h"
#include "utils.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;


SCENARIO ( "Find reads in the graph using a Traverser (exact)", "[traverser]" )
{
  GIVEN ( "A small variation graph and a set of reads" )
  {
    typedef seqan::IndexWotd<> TIndexSpec;
    typedef typename Traverser< TIndexSpec, BFS, ExactMatching >::Type TTraverser;

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

    Dna5QRecords reads;
    readRecords( reads, reads_file, 10 );
    Dna5QStringSetIndex< TIndexSpec > reads_index( reads.str );

    unsigned int seed_len = 10;

    WHEN ( "Run a traverser on all loci with seed length " + std::to_string( seed_len ) )
    {
      TTraverser traverser( &vargraph, &reads_index, seed_len );

      unsigned int counter = 0;
      std::size_t truth[10][2] = { {1, 0}, {1, 1}, {9, 4}, {9, 17}, {16, 0}, {17, 0},
        {20, 0}, {20, 31}, {20, 38}, {20, 38} };

      std::function< void( typename TTraverser::output_type const& ) > count_hits =
        [&counter, &truth]( typename TTraverser::output_type const& hit ) {
          if ( counter < 10 ) {
            REQUIRE( beginPositionH( hit ) == truth[ counter ][0] );
            REQUIRE( endPositionH( hit ) == truth[ counter ][1] );
            REQUIRE( beginPositionV( hit ) == counter );
            REQUIRE( endPositionV( hit ) == 0 );
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
          VarGraph::offset_type seqlen = vargraph.node_sequence( node_id ).length();
          for ( VarGraph::offset_type f = 0; f < seqlen; ++f ) {
            traverser.set_start_locus( node_id, f );
            traverser.run( count_hits );
          }
        }
      }
    }
  }
}

SCENARIO ( "Serialize/deserialize paths nodes coverage into/from the file", "[traverser]" )
{
  GIVEN ( "Nodes coverage of two paths" )
  {
    std::vector< VarGraph::NodeCoverage > paths_node_coverage;
    VarGraph::NodeCoverage covered_nodes;
    unsigned int paths_num = 2;
    std::string file_path = "/tmp/test_path_coverage";

    for ( unsigned int i = 0; i < paths_num; ++i ) {
      for ( VarGraph::nodeid_type j = 100*(i+1); j < 100*(i+1) + 12; ++j ) {
        covered_nodes.insert ( j );
      }
      paths_node_coverage.push_back ( covered_nodes );
      covered_nodes.clear();
    }

    WHEN ( "Serialize it to the file" )
    {
      save_paths_coverage ( paths_node_coverage, file_path );

      THEN ( "Deserializing should yield the same coverage" )
      {
        std::vector< VarGraph::NodeCoverage > paths_coverage_reloaded;
        load_paths_coverage ( paths_coverage_reloaded, file_path, paths_num );

        std::vector< VarGraph::nodeid_type > sorted_node_ids;
        REQUIRE ( paths_coverage_reloaded.size() == paths_num );
        for ( unsigned int i = 0; i < paths_num; ++i ) {
          REQUIRE ( paths_coverage_reloaded[i].size() == 12 );
          for ( const auto &node_id : paths_coverage_reloaded[i] ) {
            sorted_node_ids.push_back ( node_id );
          }
          unsigned int j = 0;
          std::sort ( sorted_node_ids.begin(), sorted_node_ids.end() );
          for ( const auto &node_id : sorted_node_ids ) {
            REQUIRE ( node_id == 100*(i+1) + j );
            ++j;
          }
          sorted_node_ids.clear();
        }
      }
    }
  }
}
