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

#include "vargraph.h"
#include "pathset.h"
#include "logger.h"

#include "test_base.hpp"


using namespace grem;

SCENARIO( "PathSet provides an interface similar to a conventional container", "[pathindex]" )
{
  auto basic_tests = []( auto&& set, const VarGraph& vargraph ) {
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

    typedef Path< VarGraph > TPath;
    typedef typename TPath::nodes_type::value_type TNodeID;

    std::vector< TNodeID > nodes( 100 );
    std::iota( nodes.begin(), nodes.end(), 1 );
    Path< VarGraph > path( &vargraph, nodes );
    REQUIRE( set.found( path ) );

    nodes = std::vector< TNodeID >( 12 );
    std::iota( nodes.begin(), nodes.end(), 94 );
    path = Path< VarGraph >( &vargraph, nodes );
    REQUIRE( set.found( path ) );

    nodes = std::vector< TNodeID >( 12 );
    std::iota( nodes.begin(), nodes.end(), 194 );
    path = Path< VarGraph >( &vargraph, nodes );
    REQUIRE( !covered_by( path, set ) );

    nodes = std::vector< TNodeID >( 1 );
    std::iota( nodes.begin(), nodes.end(), 210 );
    path = Path< VarGraph >( &vargraph, nodes );
    REQUIRE( covered_by( path, set ) );
  };

  GIVEN( "A small VarGraph" )
  {
    std::string vgpath = test_data_dir + "/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph;
    vargraph.load( gifs );

    GIVEN( "A PathSet containing some paths" )
    {
      typedef Path< VarGraph, Compact > TPath;
      typedef typename TPath::nodes_type::value_type TNodeID;

      PathSet< TPath > set;

      WHEN( "The paths are added" )
      {
        std::vector< TNodeID > nodes( 100 );
        std::iota( nodes.begin(), nodes.end(), 1 );
        TPath path( &vargraph, std::move( nodes ) );
        set.push_back( path );
        nodes = std::vector< TNodeID >( 12 );
        std::iota( nodes.begin(), nodes.end(), 43 );
        path = TPath( &vargraph, std::move( nodes ) );
        set.push_back( path );
        nodes = std::vector< TNodeID >( 200 );
        std::iota( nodes.begin(), nodes.end(), 1 );
        path = TPath( &vargraph, std::move( nodes ) );
        set.push_back( path );
        nodes = std::vector< TNodeID >( 11 );
        std::iota( nodes.begin(), nodes.end(), 200 );
        path = TPath( &vargraph, std::move( nodes ) );
        set.push_back( path );

        THEN( "It should pass the basic tests" )
        {
          basic_tests( set, vargraph );
        }

        AND_WHEN( "The PathSet is moved by assignment" )
        {
          PathSet< TPath > another_set;
          another_set = std::move( set );

          THEN( "The moved PathSet should pass the basic tests" )
          {
            basic_tests( another_set, vargraph );
          }
        }

        AND_WHEN( "Another PathSet is constructed by moving" )
        {
          PathSet< TPath > another_set( std::move( set ) );

          THEN( "The moved PathSet should pass the basic tests" )
          {
            basic_tests( another_set, vargraph );
          }
        }

        AND_WHEN( "The PathSet is serialised to a output stream" )
        {
          std::string tmpfpath = get_tmpfile();
          save( set, tmpfpath );
          PathSet< TPath > another_set;
          open( another_set, &vargraph, tmpfpath );

          THEN( "The loaded PathSet should pass the basic tests" )
          {
            basic_tests( another_set, vargraph );
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
