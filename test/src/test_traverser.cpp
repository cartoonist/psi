/**
 *    @file  test_traverser.cpp
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

#include <psi/sequence.hpp>
#include <psi/index.hpp>
#include <psi/graph.hpp>
#include <psi/traverser.hpp>
#include <psi/utils.hpp>
#include <gum/seqgraph.hpp>
#include <gum/io_utils.hpp>
#include <seqan/seq_io.h>

#include "test_base.hpp"


using namespace psi;

SCENARIO ( "Find reads in the graph using a Traverser (exact)", "[traverser]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;
  typedef graph_type::offset_type offset_type;

  GIVEN ( "A small graph and a set of reads" )
  {
    typedef seqan::IndexWotd<> TIndexSpec;
    typedef seqan::Index< Dna5QStringSet<>, TIndexSpec > TIndex;

    std::string vgpath = test_data_dir + "/small/x.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );
    std::string readspath = test_data_dir + "/small/reads_n10l10e0i0.fastq";
    seqan::SeqFileIn reads_file;
    if ( !open( reads_file, readspath.c_str() ) ) {
      throw std::runtime_error( "cannot open file " + readspath );
    }

    Records< Dna5QStringSet<> > reads;
    readRecords( reads, reads_file, 10 );
    seqan::Index< Dna5QStringSet<>, TIndexSpec > reads_index( reads.str );

    unsigned int seed_len = 10;

    WHEN ( "Run a DFS traverser on all loci with seed length " + std::to_string( seed_len ) )
    {
      typedef typename Traverser< graph_type, TIndex, DFS, ExactMatching >::Type TTraverser;

      TTraverser traverser( &graph, &reads, &reads_index, seed_len );

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
        for ( std::size_t r = 1; r <= graph.get_node_count(); ++r ) {
          const auto& node_id = graph.rank_to_id( r );
          offset_type seqlen = graph.node_length( node_id );
          for ( offset_type f = 0; f < seqlen; ++f ) {
            traverser.add_locus( node_id, f );
            traverser.run( count_hits );
          }
        }
      }
    }

    WHEN ( "Run a BFS traverser on all loci with seed length " + std::to_string( seed_len ) )
    {
      typedef typename Traverser< graph_type, TIndex, BFS, ExactMatching >::Type TTraverser;

      TTraverser traverser( &graph, &reads, &reads_index, seed_len );

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
        for ( std::size_t r = 1; r <= graph.get_node_count(); ++r ) {
          const auto& node_id = graph.rank_to_id( r );
          offset_type seqlen = graph.node_length( node_id );
          for ( offset_type f = 0; f < seqlen; ++f ) {
            traverser.add_locus( node_id, f );
            traverser.run( count_hits );
          }
        }
      }
    }
  }
}
