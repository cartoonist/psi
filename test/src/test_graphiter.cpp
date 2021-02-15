/**
 *    @file  test_graphiter.cpp
 *   @brief  Graph iterator test cases.
 *
 *  Contains test cases for graph iterator header file.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Feb 13, 2017  22:24
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <algorithm>
#include <iterator>
#include <random>
#include <chrono>
#include <libgen.h>
#include <string>
#include <vector>

#include <psi/graph.hpp>
#include <psi/graph_iter.hpp>
#include <psi/pathindex.hpp>
#include <gum/seqgraph.hpp>
#include <gum/io_utils.hpp>

#include "test_base.hpp"


using namespace psi;

SCENARIO( "Get unique full haplotype using Haplotyper graph iterator", "[graph][iterator][haplotyper]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;

  GIVEN( "A tiny variation graph" )
  {
    std::string vgpath = test_data_dir + "/tiny/tiny.gfa";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );

    WHEN( "the eigth haplotypes are generated using Haplotyper" )
    {
      auto hap_itr = begin( graph, Haplotyper<>() );
      auto hap_end = end( graph, Haplotyper<>() );

      Path< graph_type > haplotype1( &graph );
      Path< graph_type > haplotype2( &graph );
      Path< graph_type > haplotype3( &graph );
      Path< graph_type > haplotype4( &graph );
      Path< graph_type > haplotype5( &graph );
      Path< graph_type > haplotype6( &graph );
      Path< graph_type > haplotype7( &graph );
      Path< graph_type > haplotype8( &graph );
      get_uniq_full_haplotype( haplotype1, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype2, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype3, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype4, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype5, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype6, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype7, hap_itr, hap_end, 1 );
      get_uniq_full_haplotype( haplotype8, hap_itr, hap_end, 1 );

      THEN( "they should be unique" )
      {
        std::string hapstr1 = sequence( haplotype1 );
        std::string hapstr2 = sequence( haplotype2 );
        std::string hapstr3 = sequence( haplotype3 );
        std::string hapstr4 = sequence( haplotype4 );
        std::string hapstr5 = sequence( haplotype5 );
        std::string hapstr6 = sequence( haplotype6 );
        std::string hapstr7 = sequence( haplotype7 );
        std::string hapstr8 = sequence( haplotype8 );
        unsigned int matched = 0;
        if ( hapstr1 == hapstr2 ) ++matched;
        if ( hapstr1 == hapstr3 ) ++matched;
        if ( hapstr2 == hapstr3 ) ++matched;
        if ( hapstr1 == hapstr4 ) ++matched;
        if ( hapstr2 == hapstr4 ) ++matched;
        if ( hapstr3 == hapstr4 ) ++matched;
        if ( hapstr1 == hapstr5 ) ++matched;
        if ( hapstr2 == hapstr5 ) ++matched;
        if ( hapstr3 == hapstr5 ) ++matched;
        if ( hapstr4 == hapstr5 ) ++matched;
        if ( hapstr1 == hapstr6 ) ++matched;
        if ( hapstr2 == hapstr6 ) ++matched;
        if ( hapstr3 == hapstr6 ) ++matched;
        if ( hapstr4 == hapstr6 ) ++matched;
        if ( hapstr5 == hapstr6 ) ++matched;
        if ( hapstr1 == hapstr7 ) ++matched;
        if ( hapstr2 == hapstr7 ) ++matched;
        if ( hapstr3 == hapstr7 ) ++matched;
        if ( hapstr4 == hapstr7 ) ++matched;
        if ( hapstr5 == hapstr7 ) ++matched;
        if ( hapstr6 == hapstr7 ) ++matched;
        if ( hapstr1 == hapstr8 ) ++matched;
        if ( hapstr2 == hapstr8 ) ++matched;
        if ( hapstr3 == hapstr8 ) ++matched;
        if ( hapstr4 == hapstr8 ) ++matched;
        if ( hapstr5 == hapstr8 ) ++matched;
        if ( hapstr6 == hapstr8 ) ++matched;
        if ( hapstr7 == hapstr8 ) ++matched;

        REQUIRE( matched == 0 );
      }
      AND_THEN( "they should have the correct length" )
      {
        REQUIRE( length( haplotype1 ) == 10 );
        REQUIRE( length( haplotype2 ) == 10 );
        REQUIRE( length( haplotype3 ) == 10 );
        REQUIRE( length( haplotype4 ) == 10 );
        REQUIRE( length( haplotype5 ) == 10 );
        REQUIRE( length( haplotype6 ) == 10 );
        REQUIRE( length( haplotype7 ) == 10 );
        REQUIRE( length( haplotype8 ) == 10 );
      }
      AND_THEN( "level of iterator should be the number of haplotypes" )
      {
        REQUIRE( hap_itr.level() == 8 );
      }
    }
  }

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = test_data_dir + "/small/x.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );

    WHEN( "the three haplotypes are generated using Haplotyper" )
    {
      auto hap_itr = begin( graph, Haplotyper<>() );
      auto hap_end = end( graph, Haplotyper<>() );

      Path< graph_type > haplotype1( &graph );
      Path< graph_type > haplotype2( &graph );
      Path< graph_type > haplotype3( &graph );
      get_uniq_full_haplotype( haplotype1, hap_itr, hap_end );
      get_uniq_full_haplotype( haplotype2, hap_itr, hap_end );
      get_uniq_full_haplotype( haplotype3, hap_itr, hap_end );

      THEN( "they should be unique" )
      {
        std::string hapstr1 = sequence( haplotype1 );
        std::string hapstr2 = sequence( haplotype2 );
        std::string hapstr3 = sequence( haplotype3 );
        REQUIRE( hapstr1 != hapstr2 );
        REQUIRE( hapstr2 != hapstr3 );
        REQUIRE( hapstr1 != hapstr3 );
      }
      AND_THEN( "they should have the correct length" )
      {
        REQUIRE( length( haplotype1 ) == 147 );
        REQUIRE( length( haplotype2 ) > 130 );  // randomised path.
        REQUIRE( length( haplotype3 ) > 130 );  // randomised path.
      }
      AND_THEN( "they all should cover 'merge' nodes" )
      {
        std::vector< Path< graph_type > > paths_set;
        paths_set.push_back( haplotype1 );
        paths_set.push_back( haplotype2 );
        paths_set.push_back( haplotype3 );

        REQUIRE( path_coverage( haplotype1.get_nodes().front(), paths_set ) == 3 );
        REQUIRE( path_coverage( 2, paths_set ) != 3 );
        REQUIRE( path_coverage( 6, paths_set ) == 3 );
        REQUIRE( path_coverage( 9, paths_set ) == 3 );
        REQUIRE( path_coverage( 18, paths_set ) == 3 );
        REQUIRE( path_coverage( 20, paths_set ) == 3 );
        REQUIRE( path_coverage( 210, paths_set ) == 3 );
        REQUIRE( path_coverage( 207, paths_set ) == 3 );
        REQUIRE( path_coverage( 205, paths_set ) == 3 );
        REQUIRE( path_coverage( 202, paths_set ) == 3 );
        REQUIRE( path_coverage( 200, paths_set ) == 3 );
        REQUIRE( path_coverage( 96, paths_set ) == 3 );
        REQUIRE( path_coverage( 99, paths_set ) == 3 );
        REQUIRE( path_coverage( 101, paths_set ) == 3 );
        REQUIRE( path_coverage( 104, paths_set ) == 3 );
        REQUIRE( path_coverage( haplotype1.get_nodes().back(), paths_set ) == 3 );
      }
      AND_THEN( "level of iterator should be the number of haplotypes" )
      {
        REQUIRE( hap_itr.level() == 3 );
      }
    }
  }
}

SCENARIO( "A Haplotyper graph iterator raise on end", "[graph][iterator][haplotyper]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;

  GIVEN( "A small variation graph and a Haplotyper iterator with `raise_on_end` enabled" )
  {
    std::string vgpath = test_data_dir + "/small/x.gfa";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );
    auto hap_itr = begin( graph, Haplotyper<>() );
    auto hap_end = end( graph, Haplotyper<>() );
    hap_itr.set_raise_on_end( true );

    WHEN( "A Haplotyper iterator reaches at end" )
    {
      auto loop =
          [&hap_itr, &hap_end]() {
            while ( hap_itr != hap_end ) ++hap_itr;
          };

      THEN( "It raise an exception if it is incremented" )
      {
        REQUIRE_THROWS( loop() );
      }
    }
  }
}

SCENARIO( "Extend a path to length k using Haplotyper graph iterator", "[graph][iterator][haplotyper]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;

  GIVEN( "A small variation graph and a Haplotyper graph iterator" )
  {
    std::string vgpath = test_data_dir + "/small/x.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );
    auto hap_itr = begin( graph, Haplotyper<>() );
    auto hap_end = end( graph, Haplotyper<>() );
    unsigned int k = 5;

    WHEN( "A path is extend to length " + std::to_string( k ) )
    {
      Path< graph_type > path( &graph );
      util::extend_to_k( path, hap_itr, hap_end, k );
      initialize(path);

      THEN( "Its length should be extended" )
      {
        REQUIRE( path.get_sequence_len() == 8 );
        REQUIRE( position_to_id( path, 7 ) == 1 );
        REQUIRE( position_to_offset( path, 7 ) == 7 );
      }
    }

    k = 14;

    WHEN( "A path is extend to length " + std::to_string( k ) )
    {
      Path< graph_type > path( &graph );
      util::extend_to_k( path, hap_itr, hap_end, k );
      initialize(path);

      THEN( "Its length should be extended" )
      {
        REQUIRE( path.get_sequence_len() == 14 );
        REQUIRE( position_to_id( path, 13 ) == 7 );
        REQUIRE( position_to_offset( path, 13 ) == 0 );
      }
    }
  }
}

SCENARIO( "Get unique patched haplotypes using Haplotyper graph iterator", "[graph][iterator][haplotyper]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = test_data_dir + "/small/x.gfa";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );
    unsigned int context_len = 10;

    WHEN( "Generate 32x patched haplotypes are generated using a Haplotyper iterator" )
    {
      auto hap_itr = begin( graph, Haplotyper<>() );
      auto hap_end = end( graph, Haplotyper<>() );
      std::vector< Path< graph_type > > pathset;

      for ( auto i = 0; i < 32; ++i )
        get_uniq_patched_haplotype( pathset, hap_itr, hap_end, context_len );

      for ( std::size_t i = 0; i < pathset.size(); ++i )
        initialize( pathset.at( i ) );

      THEN( "The number of patches should be in correct range" )
      {
        REQUIRE( pathset.size() >= 69 );
        REQUIRE( pathset.size() <= 76 );
      }
    }
  }
}

SCENARIO( "Traverse a sequence graph using backtracking algorithm", "[graph][iterator][backtracker]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;
  typedef graph_type::id_type id_type;

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = test_data_dir + "/small/x.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );

    unsigned int kmer_len = 20;

    WHEN( "enumerating all kmers of length " + std::to_string( kmer_len ) )
    {
      std::string truth_dir = test_data_dir + "/small/";
      std::string truth_filepath = truth_dir + std::to_string( kmer_len ) + "-mers";
      std::ifstream truth_stream( truth_filepath, std::ifstream::in );
      std::string true_kmer;
      id_type true_snode_id;
      unsigned int true_offset;

      auto bt_itr = begin( graph, Backtracker() );
      auto bt_end = end( graph, Backtracker() );
      std::vector< id_type > trav_path;
      std::string trav_seq = "";
      // :TODO:Mon May 22 11:16:\@cartoonist: add REQUIREs and assertions.
      // :TODO:Mon May 22 11:15:\@cartoonist: add a method to simulate a kmer.
      for ( unsigned int n_idx = 1; n_idx < graph.get_node_count(); ++n_idx ) {
        id_type start_node_id = graph.rank_to_id( n_idx );
        unsigned int label_len = graph.node_length( start_node_id );

        for ( unsigned int offset = 0; offset < label_len; ++offset ) {
          bt_itr.reset( start_node_id );

          while ( bt_itr != bt_end ) {
            while ( bt_itr != bt_end ) {
              trav_path.push_back( *bt_itr );
              if ( *bt_itr != start_node_id ) {
                trav_seq += graph.node_sequence( *bt_itr );
              }
              else {
                trav_seq = graph.node_sequence( *bt_itr ).substr( offset );
              }

              if ( trav_seq.length() < kmer_len ) ++bt_itr;
              else break;
            }

            const std::string &kmer = trav_seq.substr( 0, kmer_len );
            if ( kmer.length() == kmer_len ) {
              truth_stream >> true_kmer;
              truth_stream >> true_snode_id;
              truth_stream >> true_offset;
              REQUIRE( kmer == true_kmer );
              REQUIRE( start_node_id == true_snode_id );
              REQUIRE( offset == true_offset );
            }

            --bt_itr;

            unsigned int trav_len = trav_seq.length();
            id_type poped_id = 0;
            while ( !trav_path.empty() && poped_id != *bt_itr ) {
              poped_id = trav_path.back();
              trav_len -= graph.node_length( poped_id );
              trav_path.pop_back();
            }
            trav_seq = trav_seq.substr( 0, trav_len );
          }

          trav_seq.clear();
        }
      }
    }
  }
}

SCENARIO( "Sequence graph breadth-first traverse (BFS)", "[graph][iterator][bfs]" )
{
  typedef gum::SeqGraph< gum::Dynamic > graph_type;
  typedef graph_type::id_type id_type;

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = test_data_dir + "/small/x.gfa";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );

    WHEN( "traverse the graph using BFS graph iterator" )
    {
      auto bfs_itr = begin( graph, BFS() );
      auto bfs_end = end( graph, BFS() );

      THEN( "nodes should be traversed in BFS order" )
      {
        id_type truth = 1;
        while ( bfs_itr != bfs_end ) {
          REQUIRE( *bfs_itr == truth );  // The graph is such that its BFS is in order.
          ++truth;
          ++bfs_itr;
        }
        REQUIRE( truth == 211 );
      }
    }
  }

  GIVEN( "A variation graph with more than one connected component" )
  {
    std::string vgpath = test_data_dir + "/multi/multi.vg";
    graph_type graph;
    gum::util::extend( graph, vgpath, true );

    WHEN( "traverse the graph using BFS graph iterator" )
    {
      auto bfs_itr = begin( graph, BFS() );
      auto bfs_end = end( graph, BFS() );

      THEN( "nodes should be traversed in BFS order" )
      {
        id_type truth = 1;
        while ( bfs_itr != bfs_end ) {
          REQUIRE( *bfs_itr == truth );  // The graph is such that its BFS is in order.
          ++truth;
          ++bfs_itr;
        }
        REQUIRE( truth == 226 );
      }
    }
  }
}
