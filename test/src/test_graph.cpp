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

#include "vg/vg.pb.h"
#include "vg/stream.hpp"

#include "test_base.hpp"


using namespace psi;

static const gum::ExternalLoader< vg::Graph > vg_loader { []( std::istream& in ) -> vg::Graph {
    vg::Graph merged;
    std::function< void( vg::Graph& ) > handle_chunks =
      [&]( vg::Graph& other ) {
        gum::util::merge_vg( merged, static_cast< vg::Graph const& >( other ) );
      };
    stream::for_each( in, handle_chunks );
    return merged;
  } };

TEMPLATE_SCENARIO( "Get graph statistics", "[graph][interface]", gum::Dynamic, gum::Succinct )
{
  typedef TestType spec_type;
  typedef gum::SeqGraph< spec_type > graph_type;

  GIVEN( "A tiny variation graph" )
  {
    std::string vgpath = test_data_dir + "/tiny/tiny.gfa";
    graph_type graph;
    gum::util::load( graph, vgpath, vg_loader );

    WHEN( "Total number of nodes in a subgraph is counted" )
    {
      auto nof_nodes = util::node_count( graph, 5, 10 );
      THEN( "It should be equal to the number of nodes in the subgraph" )
      {
        REQUIRE( nof_nodes == 5 );
      }
    }

    WHEN( "Total number of edges in a component is counted" )
    {
      auto nof_edges = util::edge_count( graph, 5 );
      THEN( "It should be equal to the number of edges in the component" )
      {
        REQUIRE( nof_edges == 13 );
      }
    }
  }
}
