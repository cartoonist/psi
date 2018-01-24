/**
 *    @file  tests_pathset.cc
 *   @brief  Test PathSet class.
 *
 *  Test scenarios for PathSet class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Sep 17, 2017  22:41
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <string>

#include "tests_base.h"
#include "vargraph.h"
#include "pathset.h"
#include "index.h"
#include "logger.h"


using namespace grem;

SCENARIO ( "Serialize/deserialize paths set into/from the file", "[pathset]" )
{
  GIVEN ( "Two paths from a small graph" )
  {
    typedef seqan::IndexEsa<> TIndexSpec;

    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );
    Dna5QPathSet< VarGraph, TIndexSpec > paths_set;

    unsigned int paths_num = 2;
    std::string file_path = "/tmp/test_pathset";

    paths_set.reserve( paths_num );
    for ( unsigned int i = 0; i < paths_num; ++i ) {
      Path< VarGraph > path( &vargraph );
      for ( VarGraph::nodeid_type j = 3+i; j <= 210; j+=(i+1)*4 ) {
        add_node( path, j );
      }
      paths_set.add_path( std::move( path ) );
    }

    WHEN ( "Serialize it to the file" )
    {
      paths_set.save( file_path );

      THEN ( "Deserializing should yield the same paths" )
      {
        Dna5QPathSet< VarGraph, TIndexSpec > loaded_paths;
        loaded_paths.load( file_path, &vargraph );
        REQUIRE ( loaded_paths.size() == paths_num );
        REQUIRE( length( loaded_paths.paths_set.at(0) ) == 52 );
        REQUIRE( length( loaded_paths.paths_set.at(1) ) == 26 );

        for ( unsigned int i = 0; i < paths_num; ++i ) {
          unsigned int j = 0;
          for ( const auto &node_id : loaded_paths.paths_set.at(i).get_nodes() ) {
            REQUIRE ( node_id == 3+i+(i+1)*4*j );
            ++j;
          }
        }
      }
    }
  }
}

SCENARIO( "Get node ID/offset by position in the PathSet", "[pathset]" )
{
  GIVEN( "A VarGraph and PathSet within the graph" )
  {
    typedef seqan::IndexEsa<> TIndexSpec;

    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    Dna5QPathSet< VarGraph, TIndexSpec > paths_set;
    Path< VarGraph > path( &vargraph, { 205, 207, 209, 210 } );

    REQUIRE( path.get_sequence_len() == 54 );
    REQUIRE( length( path ) == 4 );

    paths_set.add_path( std::move( path ) );

    WHEN( "Query the positions in the path sequence" )
    {
      THEN( "The corresponding node ID/offset should be returned" )
      {
        REQUIRE( position_to_id( paths_set, { 0, 0 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 0 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 14 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 14 } ) == 14 );
        REQUIRE( position_to_id( paths_set, { 0, 26 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 26 } ) == 26 );
        REQUIRE( position_to_id( paths_set, { 0, 27 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 27 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 30 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 30 } ) == 3 );
        REQUIRE( position_to_id( paths_set, { 0, 51 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 51 } ) == 24 );
        REQUIRE( position_to_id( paths_set, { 0, 52 } ) == 209 );
        REQUIRE( position_to_offset( paths_set, { 0, 52 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 53 } ) == 210 );
        REQUIRE( position_to_offset( paths_set, { 0, 53 } ) == 0 );
      }
    }
  }
}

SCENARIO( "Compress a PathSet", "[pathset]" )
{
  GIVEN( "A Vargraph and a PathSet containing sparse patches" )
  {
    typedef seqan::IndexEsa<> TIndexSpec;

    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    Dna5QPathSet< VarGraph, TIndexSpec > paths_set;
    Path< VarGraph > path( &vargraph, { 1, 2, 4, 6, 112, 123, 135, 200 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 1, 2, 4 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 123, 135, 200 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 3, 7, 12, 39 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 38, 45, 47, 87, 99 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 100, 190, 200, 205, 210 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 29, 100, 120, 130, 140 } );
    paths_set.add_path( std::move( path ) );
    path = Path< VarGraph >( &vargraph, { 150, 160 } );
    paths_set.add_path( std::move( path ) );

    std::vector< Path< VarGraph > > truth;
    truth.push_back( Path< VarGraph >( &vargraph, { 1, 2, 4, 6, 112, 123, 135, 200 } ) );
    truth.push_back( Path< VarGraph >( &vargraph, { 1, 2, 4, 123, 135, 200 } ) );
    truth.push_back( Path< VarGraph >( &vargraph, { 3, 7, 12, 39 } ) );
    truth.push_back( Path< VarGraph >( &vargraph, { 38, 45, 47, 87, 99, 100, 190, 200, 205, 210 } ) );
    truth.push_back( Path< VarGraph >( &vargraph, { 29, 100, 120, 130, 140, 150, 160 } ) );

    WHEN( "It is compressed" )
    {
      std::vector< Path< VarGraph > > compressed;
      compress( paths_set, compressed );

      THEN( "The non-overlapping patches with non-decreasing node IDs should be combined" )
      {
        REQUIRE( compressed.size() == 5 );
        for ( unsigned int i = 0; i < truth.size(); ++i ) {
          for ( unsigned int j = 0; j < truth[i].get_nodes().size(); ++j ) {
            REQUIRE( compressed[i].get_nodes()[j] == truth[i].get_nodes()[j] );
          }
        }
      }
    }
  }
}
