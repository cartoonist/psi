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

SCENARIO( "String set of PathSet with non-zero context", "[pathset]" )
{
  GIVEN( "A VarGraph" )
  {
    typedef seqan::IndexEsa<> TIndexSpec;

    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath.c_str() );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    WHEN( "A set of paths are added to a PathSet with non-zero context in lazy mode" )
    {
      uint64_t context = 10;
      Dna5QPathSet< VarGraph, TIndexSpec, Forward > paths_set( context, true );
      Path< VarGraph > path( &vargraph, { 205, 207, 209, 210 } );
      paths_set.add_path( ( path ) );
      path = Path< VarGraph >( &vargraph, { 187, 189, 191, 193, 194, 195, 197 });
      paths_set.add_path( ( path ) );
      path = Path< VarGraph >( &vargraph, { 167, 168, 171, 172, 174 });
      paths_set.add_path( ( path ) );

      paths_set.create_index();  /**< @brief Paths are added with this call. */

      THEN( "The paths sequence set should be trimmed according to context" )
      {
        REQUIRE( paths_set.string_set[0] == "GTTTCCTGTACTAAGGACAAAGGTGCGGGGAGATAA" );
        REQUIRE( paths_set.string_set[1] == "CAAGGGCTTTTAA" );
        REQUIRE( paths_set.string_set[2] == "CATTTGTCTTATTGTCCAGGA" );
      }
    }

    WHEN( "A set of paths are added to a PathSet with non-zero context in Forward direction" )
    {
      uint64_t context = 10;
      Dna5QPathSet< VarGraph, TIndexSpec, Forward > paths_set( context );
      Path< VarGraph > path( &vargraph, { 205, 207, 209, 210 } );
      paths_set.add_path( ( path ) );
      path = Path< VarGraph >( &vargraph, { 187, 189, 191, 193, 194, 195, 197 });
      paths_set.add_path( ( path ) );
      path = Path< VarGraph >( &vargraph, { 167, 168, 171, 172, 174 });
      paths_set.add_path( ( path ) );

      THEN( "The paths sequence set should be trimmed according to context" )
      {
        REQUIRE( paths_set.string_set[0] == "GTTTCCTGTACTAAGGACAAAGGTGCGGGGAGATAA" );
        REQUIRE( paths_set.string_set[1] == "CAAGGGCTTTTAA" );
        REQUIRE( paths_set.string_set[2] == "CATTTGTCTTATTGTCCAGGA" );
      }

      THEN( "The node ID/offset by position in the PathSet should be shifted according to the context" )
      {
        REQUIRE( position_to_id( paths_set, { 0, 0 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 0 } ) == 18 );
        REQUIRE( position_to_id( paths_set, { 0, 8 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 8 } ) == 26 );
        REQUIRE( position_to_id( paths_set, { 0, 9 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 9 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 12 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 12 } ) == 3 );
        REQUIRE( position_to_id( paths_set, { 0, 33 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 33 } ) == 24 );
        REQUIRE( position_to_id( paths_set, { 0, 34 } ) == 209 );
        REQUIRE( position_to_offset( paths_set, { 0, 34 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 35 } ) == 210 );
        REQUIRE( position_to_offset( paths_set, { 0, 35 } ) == 0 );

        REQUIRE( position_to_id( paths_set, { 1, 0 } ) == 187 );
        REQUIRE( position_to_offset( paths_set, { 1, 0 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 1 } ) == 189 );
        REQUIRE( position_to_offset( paths_set, { 1, 1 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 2 } ) == 191 );
        REQUIRE( position_to_offset( paths_set, { 1, 2 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 5 } ) == 191 );
        REQUIRE( position_to_offset( paths_set, { 1, 5 } ) == 3 );
        REQUIRE( position_to_id( paths_set, { 1, 6 } ) == 191 );
        REQUIRE( position_to_offset( paths_set, { 1, 6 } ) == 4 );
        REQUIRE( position_to_id( paths_set, { 1, 7 } ) == 193 );
        REQUIRE( position_to_offset( paths_set, { 1, 7 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 8 } ) == 194 );
        REQUIRE( position_to_offset( paths_set, { 1, 8 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 9 } ) == 195 );
        REQUIRE( position_to_offset( paths_set, { 1, 9 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 10 } ) == 197 );
        REQUIRE( position_to_offset( paths_set, { 1, 10 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 1, 11 } ) == 197 );
        REQUIRE( position_to_offset( paths_set, { 1, 11 } ) == 1 );
        REQUIRE( position_to_id( paths_set, { 1, 12 } ) == 197 );
        REQUIRE( position_to_offset( paths_set, { 1, 12 } ) == 2 );

        REQUIRE( position_to_id( paths_set, { 2, 0 } ) == 167 );
        REQUIRE( position_to_offset( paths_set, { 2, 0 } ) == 21 );
        REQUIRE( position_to_id( paths_set, { 2, 8 } ) == 167 );
        REQUIRE( position_to_offset( paths_set, { 2, 8 } ) == 29 );
        REQUIRE( position_to_id( paths_set, { 2, 9 } ) == 168 );
        REQUIRE( position_to_offset( paths_set, { 2, 9 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 2, 10 } ) == 171 );
        REQUIRE( position_to_offset( paths_set, { 2, 10 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 2, 11 } ) == 172 );
        REQUIRE( position_to_offset( paths_set, { 2, 11 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 2, 12 } ) == 174 );
        REQUIRE( position_to_offset( paths_set, { 2, 12 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 2, 20 } ) == 174 );
        REQUIRE( position_to_offset( paths_set, { 2, 20 } ) == 8 );
      }
    }

    WHEN( "A set of paths are added to a PathSet with non-zero context in Reversed direction" )
    {
      uint64_t context = 10;
      Dna5QPathSet< VarGraph, TIndexSpec, Reversed > paths_set( context );
      Path< VarGraph > path( &vargraph, { 205, 207, 209, 210 } );
      paths_set.add_path( ( path ) );
      path = Path< VarGraph >( &vargraph, { 187, 189, 191, 193, 194, 195, 197 });
      paths_set.add_path( ( path ) );
      path = Path< VarGraph >( &vargraph, { 167, 168, 171, 172, 174 });
      paths_set.add_path( ( path ) );

      THEN( "The paths sequence set should be trimmed according to context" )
      {
        REQUIRE( paths_set.string_set[0] == "AATAGAGGGGCGTGGAAACAGGAATCATGTCCTTTG" );
        REQUIRE( paths_set.string_set[1] == "AATTTTCGGGAAC" );
        REQUIRE( paths_set.string_set[2] == "AGGACCTGTTATTCTGTTTAC" );
      }

      THEN( "The node ID/offset by position in the PathSet should be shifted according to the context" )
      {
        REQUIRE( position_to_id( paths_set, { 0, 0 } ) == 210 );
        REQUIRE( position_to_offset( paths_set, { 0, 0 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 1 } ) == 209 );
        REQUIRE( position_to_offset( paths_set, { 0, 1 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 2 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 2 } ) == 24 );
        REQUIRE( position_to_id( paths_set, { 0, 20 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 20 } ) == 6 );
        REQUIRE( position_to_id( paths_set, { 0, 26 } ) == 207 );
        REQUIRE( position_to_offset( paths_set, { 0, 26 } ) == 0 );
        REQUIRE( position_to_id( paths_set, { 0, 27 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 27 } ) == 26 );
        REQUIRE( position_to_id( paths_set, { 0, 29 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 29 } ) == 24 );
        REQUIRE( position_to_id( paths_set, { 0, 35 } ) == 205 );
        REQUIRE( position_to_offset( paths_set, { 0, 35 } ) == 18 );

        // :TODO:Thu Jan 25 01:14:\@cartoonist: Add test cases for paths 1 and 2.
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
