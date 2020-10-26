/**
 *    @file  test_seedfinder.cpp
 *   @brief  SeedFinder class test cases.
 *
 *  This test contains test scenarios for SeedFinder class.
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

#include <gum/io_utils.hpp>
#include <psi/seed_finder.hpp>
#include <psi/traverser.hpp>
#include <psi/utils.hpp>

#include "test_base.hpp"


using namespace psi;

SCENARIO ( "Pick genome-wide paths", "[seedfinder]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    typedef gum::SeqGraph< gum::Dynamic > graph_type;
    typedef SeedFinderTraits< gum::Dynamic, Dna5QStringSet<>, seqan::IndexEsa<>, InMemory > finder_traits_type;

    std::string vgpath = test_data_dir + "/tiny/tiny.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath );

    SeedFinder< NoStats, finder_traits_type > finder( graph, 30 );
    finder.unset_as_finaliser();

    unsigned int nof_paths = 4;
    WHEN( "Some paths " + std::to_string( nof_paths ) + " are picked using a SeedFinder" )
    {
      finder.pick_paths( nof_paths, false );
      auto const& paths = finder.get_pindex().get_paths_set();
      REQUIRE( paths.size() == nof_paths );
      THEN( "The paths should be unique" )
      {
        REQUIRE( sequence( paths[0] ) ==
            "CAAATAAGATTTGAAAATTTTCTGGAGTTCTATAATATACCAACTCTCTG" );
        REQUIRE( sequence( paths[1] ) ==
            "CAAATAAGGCTTGGAAATTTTCTGGAGTTCTATTATATTCCAACTCTCTG" );
        REQUIRE( sequence( paths[2] ) != sequence( paths[0] ) );
        REQUIRE( sequence( paths[2] ) != sequence( paths[1] ) );
        REQUIRE( sequence( paths[3] ) != sequence( paths[0] ) );
        REQUIRE( sequence( paths[3] ) != sequence( paths[1] ) );
        REQUIRE( sequence( paths[3] ) != sequence( paths[2] ) );
      }
    }
  }
}

SCENARIO ( "Add starting loci when using paths index", "[seedfinder]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    typedef gum::SeqGraph< gum::Dynamic > graph_type;
    typedef graph_type::id_type id_type;
    typedef graph_type::offset_type offset_type;
    typedef SeedFinderTraits< gum::Dynamic, Dna5QStringSet<>, seqan::IndexEsa<> > finder_traits_type;

    std::string vgpath = test_data_dir + "/tiny/tiny.gfa";
    graph_type graph;
    gum::util::extend( graph, vgpath );

    unsigned char k = 12;
    unsigned char nof_paths = 4;
    WHEN( "Using " + std::to_string( nof_paths ) + " number of paths" )
    {
      std::vector< std::pair< id_type, offset_type > > truth;
      truth.push_back( std::make_pair( 1, 2 ) );
      truth.push_back( std::make_pair( 1, 3 ) );
      truth.push_back( std::make_pair( 1, 4 ) );
      truth.push_back( std::make_pair( 1, 5 ) );
      truth.push_back( std::make_pair( 1, 6 ) );
      truth.push_back( std::make_pair( 1, 7 ) );
      truth.push_back( std::make_pair( 2, 0 ) );
      truth.push_back( std::make_pair( 3, 0 ) );
      auto truth_itr = truth.begin();

      SeedFinder< WithStats, finder_traits_type > finder( graph, k );
      finder.unset_as_finaliser();
      finder.pick_paths( nof_paths, true, k );
      finder.index_paths();

      THEN( "Starting loci must have at least one uncovered " + std::to_string( k ) + "-path" )
      {
        finder.add_uncovered_loci( );
        for ( const auto& locus : finder.get_starting_loci() ) {
          REQUIRE( locus.node_id() == (*truth_itr).first );
          REQUIRE( locus.offset() == (*truth_itr).second );
          ++truth_itr;
        }
      }
    }

    nof_paths = 8;
    WHEN( "Using " + std::to_string( nof_paths ) + " number of paths" )
    {
      SeedFinder< NoStats, finder_traits_type > finder( graph, k );
      finder.unset_as_finaliser();
      finder.pick_paths( nof_paths, true, k );
      finder.index_paths();

      THEN( "All loci should be covered by path index" )
      {
        finder.add_uncovered_loci( );
        size_t nof_loci = finder.get_starting_loci().size();
        REQUIRE( nof_loci == 0 );
      }
    }

    k = 45;
    nof_paths = 32;
    WHEN( "Using " + std::to_string( nof_paths ) + " number of paths" )
    {
      SeedFinder< NoStats, finder_traits_type > finder( graph, k );
      finder.unset_as_finaliser();
      finder.pick_paths( nof_paths, false );
      finder.index_paths();

      THEN( "All loci should be covered by path index" )
      {
        finder.add_uncovered_loci( );
        size_t nof_loci = finder.get_starting_loci().size();
        REQUIRE( nof_loci == 0 );
      }
    }
  }
}

SCENARIO( "Load and save starting loci", "[seedfinder]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    typedef gum::SeqGraph< gum::Dynamic > graph_type;
    typedef SeedFinderTraits< gum::Dynamic, Dna5QStringSet<>, seqan::IndexEsa<> > finder_traits_type;

    std::string vgpath = test_data_dir + "/tiny/tiny.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath );

    GIVEN( "A SeedFinder on this graph with known starting loci" )
    {
      int k = 12;
      int e = 10;
      SeedFinder< WithStats, finder_traits_type > finder( graph, k );
      finder.unset_as_finaliser();

      for ( int i = 325; i > 0; i -= 4 ) {
        finder.add_start( i, i % 17 );
      }

      WHEN( "It is saved to the file" )
      {
        std::string prefix = get_tmpfile();
        finder.save_starts( prefix, k, e );
        finder.set_starting_loci( {} );

        THEN( "It should be loaded as it was" )
        {
          finder.open_starts( prefix, k, e );
          REQUIRE( finder.get_starting_loci().size() == 82 );

          int i = 325;
          for ( const auto& l : finder.get_starting_loci() ) {
            REQUIRE( l.node_id() == i );
            REQUIRE( l.offset() == i % 17 );
            i -= 4;
          }
        }
      }
    }
  }
}

SCENARIO( "Distance constraints verification", "[seedfinder]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    typedef gum::SeqGraph< gum::Succinct > graph_type;
    typedef typename graph_type::id_type id_type;
    typedef typename graph_type::offset_type offset_type;
    typedef SeedFinderTraits< gum::Succinct, Dna5QStringSet<>, seqan::IndexEsa<>, InMemory > finder_traits_type;
    typedef SeedFinder< NoStats, finder_traits_type > finder_type;
    typedef std::tuple< id_type, offset_type, id_type, offset_type > ends_type;

    std::string vgpath = test_data_dir + "/tiny/tiny.vg";
    graph_type graph;
    gum::util::load( graph, vgpath );

    unsigned int dlen = 10;
    unsigned int drad = 2;
    unsigned int seedlen = 30;
    finder_type finder( graph, seedlen );
    finder.unset_as_finaliser();

    auto ibyc =
        [&graph]( id_type cid ) {
          return graph.id_by_coordinate( cid );
        };

    std::vector< ends_type > distant =
        { { ibyc( 1 ), 0, ibyc( 1 ), 0 },
          { ibyc( 1 ), 0, ibyc( 1 ), 1 },
          { ibyc( 1 ), 0, ibyc( 1 ), 3 },
          { ibyc( 1 ), 0, ibyc( 1 ), 6 },
          { ibyc( 1 ), 0, ibyc( 1 ), 7 },
          { ibyc( 1 ), 0, ibyc( 7 ), 0 },
          { ibyc( 2 ), 0, ibyc( 9 ), 10 },
          { ibyc( 9 ), 1, ibyc( 9 ), 14 },
          { ibyc( 9 ), 5, ibyc( 9 ), 18 },
          { ibyc( 9 ), 18, ibyc( 11 ), 0 },
          { ibyc( 9 ), 18, ibyc( 11 ), 3 },
          { ibyc( 9 ), 18, ibyc( 15 ), 0 },
          { ibyc( 9 ), 18, ibyc( 15 ), 6 },
        };

    std::vector< ends_type > closed =
        { { ibyc( 1 ), 0, ibyc( 2 ), 0 },
          { ibyc( 1 ), 0, ibyc( 6 ), 0 },
          { ibyc( 1 ), 0, ibyc( 6 ), 2 },
          { ibyc( 9 ), 0, ibyc( 9 ), 8 },
          { ibyc( 9 ), 1, ibyc( 9 ), 13 },
          { ibyc( 9 ), 10, ibyc( 9 ), 18 },
          { ibyc( 9 ), 6, ibyc( 9 ), 18 },
          { ibyc( 9 ), 18, ibyc( 15 ), 1 },
          { ibyc( 9 ), 18, ibyc( 15 ), 5 },
        };

    WHEN( "Creating distance index" )
    {
      finder.create_distance_index( dlen, drad );

      THEN( "It should rejects nodes not complying with distance constraints" )
      {
        for ( auto ends : distant ) {
          REQUIRE( !finder.verify_distance( std::get<0>( ends ), std::get<1>( ends ),
                                            std::get<2>( ends ), std::get<3>( ends ) ) );
        }
      }

      THEN( "It should accept nodes complying with distance constraints" )
      {
        for ( auto ends : closed ) {
          REQUIRE( finder.verify_distance( std::get<0>( ends ), std::get<1>( ends ),
                                           std::get<2>( ends ), std::get<3>( ends ) ) );
        }
      }

      AND_WHEN( "It the index is loaded from disk" )
      {
        std::string prefix = get_tmpfile();
        finder.save_distance_index( prefix, dlen, drad );
        finder_type finder2( graph, seedlen );
        finder2.unset_as_finaliser();
        finder2.open_distance_index( prefix, dlen, drad );

        THEN( "It should rejects nodes not complying with distance constraints" )
        {
          for ( auto ends : distant ) {
            REQUIRE( !finder2.verify_distance( std::get<0>( ends ), std::get<1>( ends ),
                                               std::get<2>( ends ), std::get<3>( ends ) ) );
          }
        }

        THEN( "It should accept nodes complying with distance constraints" )
        {
          for ( auto ends : closed ) {
            REQUIRE( finder2.verify_distance( std::get<0>( ends ), std::get<1>( ends ),
                                              std::get<2>( ends ), std::get<3>( ends ) ) );
          }
        }
      }
    }
  }
}

/* NOTE: Put all test scenarios before this one! */
SCENARIO( "Finalise SeedFinder" )
{
  SeedFinder<>::finalise();
}
