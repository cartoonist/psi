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
#include "utils.h"
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
    typedef seqan::Index< Dna5QStringSet<>, TIndexSpec > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;
    Mapper< TTraverser > mapper( &vargraph, 30 );

    unsigned int nof_paths = 4;
    Dna5QPathIndex< VarGraph, TIndexSpec > pindex;
    WHEN( "Some paths " + std::to_string( nof_paths ) + " are picked using a Mapper" )
    {
      mapper.pick_paths( pindex, nof_paths );
      REQUIRE( length( indexText( pindex.index ) ) == nof_paths );
      THEN( "The paths should be unique" )
      {
        REQUIRE( indexText( pindex.index )[0] ==
            "CAAATAAGATTTGAAAATTTTCTGGAGTTCTATAATATACCAACTCTCTG" );
        REQUIRE( indexText( pindex.index )[1] ==
            "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG" );
        REQUIRE( indexText( pindex.index )[2] != indexText( pindex.index )[0] );
        REQUIRE( indexText( pindex.index )[2] != indexText( pindex.index )[1] );
        REQUIRE( indexText( pindex.index )[3] != indexText( pindex.index )[0] );
        REQUIRE( indexText( pindex.index )[3] != indexText( pindex.index )[1] );
        REQUIRE( indexText( pindex.index )[3] != indexText( pindex.index )[2] );
      }
    }
  }
}

SCENARIO ( "Add starting loci when using paths index", "[mapper]" )
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
    typedef seqan::Index< Dna5QStringSet<>, TIndexSpec > TIndex;
    typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;

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
    Dna5QPathIndex< VarGraph, TIndexSpec > pindex;

    WHEN( "Find starting loci using " + std::to_string( nof_paths ) + " paths" )
    {
      mapper.pick_paths( pindex, nof_paths );
      mapper.add_all_loci( pindex.get_paths_set(), k );
      for ( const auto& locus : mapper.get_starting_loci() ) {
        REQUIRE( locus.node_id() == (*truth_itr).first );
        REQUIRE( locus.offset() == (*truth_itr).second );
        ++truth_itr;
      }
    }
  }
}

SCENARIO( "Load and save starting loci", "[mapper]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    std::string vgpath = _testdir + "/data/tiny/tiny.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::ifstream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    GIVEN( "A mapper on this graph with known starting loci" )
    {
      typedef seqan::IndexEsa<> TIndexSpec;
      typedef seqan::Index< Dna5QStringSet<>, TIndexSpec > TIndex;
      typedef typename Traverser< TIndex, BFS, ExactMatching >::Type TTraverser;

      int k = 12;
      int e = 10;
      Mapper< TTraverser > mapper( &vargraph, k );

      for ( int i = 325; i > 0; i -= 4 ) {
        mapper.add_start( i, i % 17 );
      }

      WHEN( "It is saved to the file" )
      {
        std::string prefix = get_tmpfile();
        mapper.save_starts( prefix, k, e );
        mapper.set_starting_loci( {} );

        THEN( "It should be loaded as it was" )
        {
          mapper.open_starts( prefix, k, e );
          REQUIRE( mapper.get_starting_loci().size() == 82 );

          int i = 325;
          for ( const auto& l : mapper.get_starting_loci() ) {
            REQUIRE( l.node_id() == i );
            REQUIRE( l.offset() == i % 17 );
            i -= 4;
          }
        }
      }
    }
  }
}
