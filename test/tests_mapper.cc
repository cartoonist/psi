/**
 *    @file  tests_mapper.cc
 *   @brief  Mapper class test cases.
 *
 *  This test contains test scenarios for Mapper class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sat Sep 23, 2017  22:04
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <string>
#include <utility>

#include "tests_base.h"
#include "mapper.h"
#include "traverser.h"
#include "logger.h"


using namespace grem;

SCENARIO ( "Pick genome-wide paths", "[mapper]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    std::string vgpath = _testdir + "/data/tiny/tiny.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::ifstream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    typedef seqan::IndexEsa<> TIndexSpec;
    typedef typename Traverser< TIndexSpec, BFS, ExactMatching >::Type TTraverser;
    Mapper< TTraverser > mapper( &vargraph, 30 );

    unsigned int nof_paths = 4;
    PathSet< TIndexSpec > paths;
    WHEN( "Some paths " + std::to_string( nof_paths ) + " are picked using a Mapper" )
    {
      mapper.pick_paths( paths, nof_paths );
      REQUIRE( length( paths.string_set ) == nof_paths );
      THEN( "The paths should be unique" )
      {
        REQUIRE( paths.string_set[0] ==
            "CAAATAAGATTTGAAAATTTTCTGGAGTTCTATAATATACCAACTCTCTG" );
        REQUIRE( paths.string_set[1] ==
            "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG" );
        REQUIRE( paths.string_set[2] != paths.string_set[0] );
        REQUIRE( paths.string_set[2] != paths.string_set[1] );
        REQUIRE( paths.string_set[3] != paths.string_set[0] );
        REQUIRE( paths.string_set[3] != paths.string_set[1] );
        REQUIRE( paths.string_set[3] != paths.string_set[2] );
      }
    }
  }
}


SCENARIO ( "Add starting points when using paths index", "[mapper]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    std::string vgpath = _testdir + "/data/tiny/tiny.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::ifstream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    typedef seqan::IndexEsa<> TIndexSpec;
    typedef typename Traverser< TIndexSpec, BFS, ExactMatching >::Type TTraverser;

    unsigned char k = 12;
    unsigned char nof_paths = 4;
    std::vector< std::pair< VarGraph::nodeid_type, VarGraph::offset_type > > truth;
    truth.push_back( std::make_pair( 1, 2 ) );
    truth.push_back( std::make_pair( 1, 3 ) );
    truth.push_back( std::make_pair( 1, 4 ) );
    truth.push_back( std::make_pair( 1, 5 ) );
    truth.push_back( std::make_pair( 1, 6 ) );
    truth.push_back( std::make_pair( 1, 7 ) );
    truth.push_back( std::make_pair( 2, 0 ) );
    truth.push_back( std::make_pair( 3, 0 ) );
    auto truth_itr = truth.begin();

    Mapper< TTraverser > mapper( &vargraph, k );
    PathSet< TIndexSpec > paths;

    WHEN( "Find starting points using " + std::to_string( nof_paths ) + " paths" )
    {
      mapper.pick_paths( paths, nof_paths );
      mapper.add_all_loci( paths, k );
      for ( const auto& locus : mapper.get_starting_loci() ) {
        REQUIRE( locus.node_id() == (*truth_itr).first );
        REQUIRE( locus.offset() == (*truth_itr).second );
        ++truth_itr;
      }
    }
  }
}