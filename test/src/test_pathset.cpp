/**
 *    @file  test_pathset.cpp
 *   @brief  Test PathSet class.
 *
 *  Test scenarios for PathSet class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Apr 02, 2018  00:02
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <string>

#include <gum/seqgraph.hpp>
#include <gum/io_utils.hpp>
#include <psi/graph.hpp>
#include <psi/pathset.hpp>

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

SCENARIO( "PathSet provides an interface similar to a conventional container", "[pathset]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;

  auto basic_tests = []( auto&& set, const graph_type& graph ) {
    REQUIRE( set.size() == 4 );
    REQUIRE( (*set.begin()).get_nodes().size() == 100 );
    REQUIRE( *(*set.begin()).get_nodes().begin() == 1 );
    REQUIRE( *((*set.begin()).get_nodes().end() - 1) == 100 );
    REQUIRE( set[ 1 ].get_nodes().size() == 12 );
    REQUIRE( *set[ 1 ].get_nodes().begin() == 43 );
    REQUIRE( *( set[ 1 ].get_nodes().end() - 1 ) == 54 );
    REQUIRE( set[ 2 ].get_nodes().size() == 200 );
    REQUIRE( *set[ 2 ].get_nodes().begin() == 1 );
    REQUIRE( *( set[ 2 ].get_nodes().end() - 1 ) == 200 );
    REQUIRE( (*(set.end()-1)).get_nodes().size() == 11 );
    REQUIRE( *(*(set.end()-1)).get_nodes().begin() == 200 );
    REQUIRE( *((*(set.end()-1)).get_nodes().end() - 1 ) == 210 );

    typedef Path< graph_type > TPath;
    typedef typename TPath::nodes_type::value_type TNodeID;

    std::vector< TNodeID > nodes( 100 );
    std::iota( nodes.begin(), nodes.end(), 1 );
    Path< graph_type > path( &graph, nodes );
    REQUIRE( set.found( path ) );

    nodes = std::vector< TNodeID >( 12 );
    std::iota( nodes.begin(), nodes.end(), 94 );
    path = Path< graph_type >( &graph, nodes );
    REQUIRE( set.found( path ) );

    nodes = std::vector< TNodeID >( 12 );
    std::iota( nodes.begin(), nodes.end(), 194 );
    path = Path< graph_type >( &graph, nodes );
    REQUIRE( !covered_by( path, set ) );

    nodes = std::vector< TNodeID >( 1 );
    std::iota( nodes.begin(), nodes.end(), 210 );
    path = Path< graph_type >( &graph, nodes );
    REQUIRE( covered_by( path, set ) );
  };

  GIVEN( "A small graph" )
  {
    std::string vgpath = test_data_dir + "/small/x.gfa";
    graph_type graph;
    gum::util::extend( graph, vgpath, vg_loader );

    GIVEN( "A PathSet containing some paths" )
    {
      typedef Path< graph_type, Compact > TPath;
      typedef typename TPath::nodes_type::value_type TNodeID;

      PathSet< TPath > set( graph );

      WHEN( "The paths are added" )
      {
        std::vector< TNodeID > nodes( 100 );
        std::iota( nodes.begin(), nodes.end(), 1 );
        TPath path( &graph, std::move( nodes ) );
        set.push_back( path );
        nodes = std::vector< TNodeID >( 12 );
        std::iota( nodes.begin(), nodes.end(), 43 );
        path = TPath( &graph, std::move( nodes ) );
        set.push_back( path );
        nodes = std::vector< TNodeID >( 200 );
        std::iota( nodes.begin(), nodes.end(), 1 );
        path = TPath( &graph, std::move( nodes ) );
        set.push_back( path );
        nodes = std::vector< TNodeID >( 11 );
        std::iota( nodes.begin(), nodes.end(), 200 );
        path = TPath( &graph, std::move( nodes ) );
        set.push_back( path );
        set.initialize();

        THEN( "It should pass the basic tests" )
        {
          basic_tests( set, graph );
        }

        AND_WHEN( "The PathSet is moved by assignment" )
        {
          PathSet< TPath > another_set( graph );
          another_set = std::move( set );

          THEN( "The moved PathSet should pass the basic tests" )
          {
            basic_tests( another_set, graph );
          }
        }

        AND_WHEN( "Another PathSet is constructed by moving" )
        {
          PathSet< TPath > another_set( std::move( set ) );

          THEN( "The moved PathSet should pass the basic tests" )
          {
            basic_tests( another_set, graph );
          }
        }

        AND_WHEN( "The PathSet is serialised to a output stream" )
        {
          std::string tmpfpath = get_tmpfile();
          save( set, tmpfpath );
          PathSet< TPath > another_set( graph );
          open( another_set, tmpfpath );

          THEN( "The loaded PathSet should pass the basic tests" )
          {
            basic_tests( another_set, graph );
          }
        }

        AND_WHEN( "It is cleared" )
        {
          set.clear();

          THEN( "The size of the PathSet should be zero" )
          {
            REQUIRE( set.size() == 0 );
          }
        }
      }
    }
  }
}
