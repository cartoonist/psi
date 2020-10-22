/**
 *    @file  test_graph.cpp
 *   @brief  Graph interface function test scenarios.
 *
 *  Contains test cases for graph interface functions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Oct 20, 2020  22:49
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <iostream>

#include <psi/graph.hpp>
#include <gum/io_utils.hpp>
#include <pairg/spgemm_utility.hpp>

#include "test_base.hpp"


using namespace psi;

SCENARIO( "Build adjacency matrix of a character graph", "[graph][interface]" )
{
  typedef gum::SeqGraph< gum::Succinct > graph_type;
  typedef typename graph_type::id_type id_type;
  typedef typename graph_type::rank_type rank_type;
  typedef typename graph_type::offset_type offset_type;
  typedef typename graph_type::linktype_type linktype_type;
  typedef pairg::matrixOps traits_type;

  Kokkos::initialize();

  GIVEN( "A tiny variation graph" )
  {
    std::string vgpath = test_data_dir + "/tiny/tiny.gfa";
    graph_type graph;
    gum::util::load( graph, vgpath );
    auto nof_nodes = util::total_nof_loci( graph );
    auto nof_edges = nof_nodes - graph.get_node_count() + graph.get_edge_count();

    WHEN( "The adjacency matrix of its corresponding character graph is build" )
    {
      auto matrix = util::adjacency_matrix( graph, traits_type() );

      THEN( "The number of columns/rows/row map should be equal to the numbers nodes" )
      {
        REQUIRE( matrix.numCols() == nof_nodes );
        REQUIRE( matrix.numRows() == nof_nodes );
        REQUIRE( matrix.graph.row_map.extent( 0 ) == nof_nodes + 1 );
      }

      THEN( "The number of entries/values should be equal to the numbers edges" )
      {
        REQUIRE( matrix.graph.entries.extent( 0 ) == nof_edges );
        REQUIRE( matrix.values.extent( 0 ) == nof_edges );
      }

      THEN( "The matrix entries should be valid" )
      {
        graph.for_each_node(
            [&graph, &matrix, nof_nodes]( rank_type rank, id_type id ) {
              offset_type char_id = gum::util::id_to_charorder( graph, id );
              for ( offset_type i = 1; i < graph.node_length( id ); ++i, ++char_id ) {
                for ( offset_type col = 0; col < nof_nodes; ++col ) {
                  if ( col == char_id + 1 ) {
                    REQUIRE( traits_type::queryValue( matrix, char_id, col ) );
                  }
                  else {
                    REQUIRE( !traits_type::queryValue( matrix, char_id, col ) );
                  }
                }
              }
              offset_type prev = 0;
              graph.for_each_edges_out(
                  id,
                  [&graph, &matrix, &prev, nof_nodes, char_id]( id_type to, linktype_type ) {
                    auto to_char_id = gum::util::id_to_charorder( graph, to );
                    REQUIRE( traits_type::queryValue( matrix, char_id, to_char_id ) );
                    for ( ; prev < to_char_id; ++prev ) {
                      REQUIRE( !traits_type::queryValue( matrix, char_id, prev ) );
                    }
                    ++prev;
                    return true;
                  } );
              return true;
            } );
      }
    }
  }

  Kokkos::finalize();
}
