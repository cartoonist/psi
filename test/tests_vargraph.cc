/**
 *    @file  tests_vargraph.cc
 *   @brief  VarGraph module test cases.
 *
 *  Contains test cases for VarGraph header file.
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

#include "tests_base.h"
#include "vargraph.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;


SCENARIO( "Loading variation graph from a vg file", "[graph][input]" )
{
  GIVEN( "A small graph" )
  {
    std::string vgpath = _testdir + "/data/small/x";
    std::function<void(VarGraph&)> x_basic_test = []( VarGraph& xgraph )
    {
      REQUIRE( xgraph.node_count == 210 );
      REQUIRE( xgraph.edge_count == 291 );
    };

    WHEN( "The format is vg")
    {
      std::ifstream ifs( vgpath + ".vg", std::ifstream::in | std::ifstream::binary );
      VarGraph vargraph( ifs, false );

      THEN( "Should pass the basic test" )
      {
        x_basic_test( vargraph );
      }
    }
    WHEN( "The format is xg")
    {
      std::ifstream ifs( vgpath + ".xg", std::ifstream::in | std::ifstream::binary );
      VarGraph vargraph( ifs );
      THEN( "Should pass the basic test" )
      {
        x_basic_test( vargraph );
      }
    }
  }
}

// VarGraph graph iterators test scenarios.
SCENARIO ( "Get unique haplotype using Haplotyper graph iterator", "[graph][iterator][haplotyper]" )
{
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/tiny/tiny.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );

    WHEN ( "the eigth haplotypes are generated using Haplotyper" )
    {
      seqan::Iterator < VarGraph, Haplotyper >::Type hap_itr (vargraph);

      std::vector < VarGraph::nodeid_type > haplotype1;
      std::vector < VarGraph::nodeid_type > haplotype2;
      std::vector < VarGraph::nodeid_type > haplotype3;
      std::vector < VarGraph::nodeid_type > haplotype4;
      std::vector < VarGraph::nodeid_type > haplotype5;
      std::vector < VarGraph::nodeid_type > haplotype6;
      std::vector < VarGraph::nodeid_type > haplotype7;
      std::vector < VarGraph::nodeid_type > haplotype8;
      get_uniq_haplotype ( haplotype1, hap_itr, 1 );
      get_uniq_haplotype ( haplotype2, hap_itr, 1 );
      get_uniq_haplotype ( haplotype3, hap_itr, 1 );
      get_uniq_haplotype ( haplotype4, hap_itr, 1 );
      get_uniq_haplotype ( haplotype5, hap_itr, 1 );
      get_uniq_haplotype ( haplotype6, hap_itr, 1 );
      get_uniq_haplotype ( haplotype7, hap_itr, 1 );
      get_uniq_haplotype ( haplotype8, hap_itr, 1 );

      THEN ( "they should be unique" )
      {
        std::string hapstr1 = vargraph.get_string ( haplotype1 );
        std::string hapstr2 = vargraph.get_string ( haplotype2 );
        std::string hapstr3 = vargraph.get_string ( haplotype3 );
        std::string hapstr4 = vargraph.get_string ( haplotype4 );
        std::string hapstr5 = vargraph.get_string ( haplotype5 );
        std::string hapstr6 = vargraph.get_string ( haplotype6 );
        std::string hapstr7 = vargraph.get_string ( haplotype7 );
        std::string hapstr8 = vargraph.get_string ( haplotype8 );
        unsigned int matched = 0;
        if ( hapstr1 == hapstr2 ) matched += pow ( 2, 0 );
        if ( hapstr1 == hapstr3 ) matched += pow ( 2, 1 );
        if ( hapstr2 == hapstr3 ) matched += pow ( 2, 2 );
        if ( hapstr1 == hapstr4 ) matched += pow ( 2, 3 );
        if ( hapstr2 == hapstr4 ) matched += pow ( 2, 4 );
        if ( hapstr3 == hapstr4 ) matched += pow ( 2, 5 );
        if ( hapstr1 == hapstr5 ) matched += pow ( 2, 6 );
        if ( hapstr2 == hapstr5 ) matched += pow ( 2, 7 );
        if ( hapstr3 == hapstr5 ) matched += pow ( 2, 8 );
        if ( hapstr4 == hapstr5 ) matched += pow ( 2, 9 );
        if ( hapstr1 == hapstr6 ) matched += pow ( 2, 10 );
        if ( hapstr2 == hapstr6 ) matched += pow ( 2, 11 );
        if ( hapstr3 == hapstr6 ) matched += pow ( 2, 12 );
        if ( hapstr4 == hapstr6 ) matched += pow ( 2, 13 );
        if ( hapstr5 == hapstr6 ) matched += pow ( 2, 14 );
        if ( hapstr1 == hapstr7 ) matched += pow ( 2, 15 );
        if ( hapstr2 == hapstr7 ) matched += pow ( 2, 16 );
        if ( hapstr3 == hapstr7 ) matched += pow ( 2, 17 );
        if ( hapstr4 == hapstr7 ) matched += pow ( 2, 18 );
        if ( hapstr5 == hapstr7 ) matched += pow ( 2, 19 );
        if ( hapstr6 == hapstr7 ) matched += pow ( 2, 20 );
        if ( hapstr1 == hapstr8 ) matched += pow ( 2, 21 );
        if ( hapstr2 == hapstr8 ) matched += pow ( 2, 22 );
        if ( hapstr3 == hapstr8 ) matched += pow ( 2, 23 );
        if ( hapstr4 == hapstr8 ) matched += pow ( 2, 24 );
        if ( hapstr5 == hapstr8 ) matched += pow ( 2, 25 );
        if ( hapstr6 == hapstr8 ) matched += pow ( 2, 26 );
        if ( hapstr7 == hapstr8 ) matched += pow ( 2, 27 );

        if ( matched != 0 ) {
          std::cerr << "[WARNING] paths are not unique: " << matched << std::endl;
          REQUIRE ( matched > 2097151 );  // at least first seven paths should be unique.
        }
        else {
          REQUIRE (  matched == 0 );
        }
      }
      AND_THEN ( "they should have the correct length" )
      {
        REQUIRE ( haplotype1.size() == 10 );
        REQUIRE ( haplotype2.size() == 10 );
        REQUIRE ( haplotype3.size() == 10 );
        REQUIRE ( haplotype4.size() == 10 );
        REQUIRE ( haplotype5.size() == 10 );
        REQUIRE ( haplotype6.size() == 10 );
        REQUIRE ( haplotype7.size() == 10 );
        REQUIRE ( haplotype8.size() == 10 );
      }
    }
  }

  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );

    WHEN ( "the three haplotypes are generated using Haplotyper" )
    {
      seqan::Iterator < VarGraph, Haplotyper >::Type hap_itr (vargraph);

      std::vector < VarGraph::nodeid_type > haplotype1;
      std::vector < VarGraph::nodeid_type > haplotype2;
      std::vector < VarGraph::nodeid_type > haplotype3;
      get_uniq_haplotype ( haplotype1, hap_itr );
      get_uniq_haplotype ( haplotype2, hap_itr );
      get_uniq_haplotype ( haplotype3, hap_itr );

      THEN ( "they should be unique" )
      {
        std::string hapstr1 = vargraph.get_string ( haplotype1 );
        std::string hapstr2 = vargraph.get_string ( haplotype2 );
        std::string hapstr3 = vargraph.get_string ( haplotype3 );
        REQUIRE ( hapstr1 != hapstr2 );
        REQUIRE ( hapstr2 != hapstr3 );
        REQUIRE ( hapstr1 != hapstr3 );
      }
      AND_THEN ( "they should have the correct length" )
      {
        REQUIRE ( haplotype1.size() == 147 );
        REQUIRE ( haplotype2.size() > 130 );  // randomised path.
        REQUIRE ( haplotype3.size() > 130 );  // randomised path.
      }
      AND_THEN ( "they all should cover 'merge' nodes" )
      {
        std::vector< VarGraph::NodeCoverage > paths_coverage;
        VarGraph::NodeCoverage cpath;
        std::copy ( haplotype1.begin(), haplotype1.end(),
            std::inserter ( cpath, cpath.end() ) );
        paths_coverage.push_back ( cpath );
        cpath.clear();
        std::copy ( haplotype2.begin(), haplotype2.end(),
            std::inserter ( cpath, cpath.end() ) );
        paths_coverage.push_back ( cpath );
        cpath.clear();
        std::copy ( haplotype3.begin(), haplotype3.end(),
            std::inserter ( cpath, cpath.end() ) );
        paths_coverage.push_back ( cpath );

        REQUIRE ( get_path_coverage ( haplotype1.front(), paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 2, paths_coverage ) != 3 );
        REQUIRE ( get_path_coverage ( 6, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 9, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 18, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 20, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 210, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 207, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 205, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 202, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 200, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 96, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 99, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 101, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( 104, paths_coverage ) == 3 );
        REQUIRE ( get_path_coverage ( haplotype1.back(), paths_coverage ) == 3 );
      }
    }
  }
}

SCENARIO ( "Traverse a variation graph using backtracking algorithm", "[graph][iterator][backtracker]" )
{
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );

    unsigned int kmer_len = 20;

    WHEN ( "enumerating all kmers of length " + std::to_string ( kmer_len ) )
    {
      std::string truth_dir = _testdir + "/data/small/";
      std::string truth_filepath = truth_dir + std::to_string ( kmer_len ) + "-mers";
      std::ifstream truth_stream ( truth_filepath, std::ifstream::in );
      std::string true_kmer;
      VarGraph::nodeid_type true_snode_id;
      unsigned int true_offset;

      seqan::Iterator< VarGraph, Backtracker >::Type bt_itr ( vargraph );
      std::vector< VarGraph::nodeid_type > trav_path;
      std::string trav_seq = "";
      // :TODO:Mon May 22 11:16:\@cartoonist: add REQUIREs and assertions.
      // :TODO:Mon May 22 11:15:\@cartoonist: add a method to simulate a kmer.
      for ( unsigned int n_idx = 1; n_idx < vargraph.max_node_rank(); ++n_idx ) {
        VarGraph::nodeid_type start_node_id = vargraph.rank_to_id( n_idx );
        unsigned int label_len = vargraph.node_sequence( start_node_id ).length();

        for ( unsigned int offset = 0; offset < label_len; ++offset ) {
          go_begin ( bt_itr, start_node_id );

          while ( !at_end( bt_itr ) ) {
            while ( !at_end( bt_itr ) ) {
              trav_path.push_back ( *bt_itr );
              if ( *bt_itr != start_node_id ) {
                trav_seq += vargraph.node_sequence( *bt_itr );
              }
              else {
                trav_seq = vargraph.node_sequence( *bt_itr ).substr( offset );
              }

              if ( trav_seq.length() < kmer_len ) ++bt_itr;
              else break;
            }

            const std::string &kmer = trav_seq.substr ( 0, kmer_len );
            if ( kmer.length() == kmer_len ) {
              truth_stream >> true_kmer;
              truth_stream >> true_snode_id;
              truth_stream >> true_offset;
              REQUIRE ( kmer == true_kmer );
              REQUIRE ( start_node_id == true_snode_id );
              REQUIRE ( offset == true_offset );
            }

            --bt_itr;

            unsigned int trav_len = trav_seq.length();
            VarGraph::nodeid_type poped_id = 0;
            while ( !trav_path.empty() && poped_id != *bt_itr ) {
              poped_id = trav_path.back();
              trav_len -= vargraph.node_sequence( poped_id ).length();
              trav_path.pop_back();
            }
            trav_seq = trav_seq.substr ( 0, trav_len );
          }

          trav_seq.clear();
        }
      }
    }
  }
}

SCENARIO ( "Variation graph breadth-first traverse (BFS)", "[graph][iterator][bfs]" )
{
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );

    WHEN ( "traverse the graph using BFS graph iterator" )
    {
      seqan::Iterator < VarGraph, BFS >::Type bfs_itr (vargraph);

      THEN ( "nodes should be traversed in BFS order" )
      {
        VarGraph::nodeid_type truth = 1;
        while ( !at_end ( bfs_itr ) ) {
          REQUIRE ( *bfs_itr == truth );  // The graph is such that its BFS is in order.
          ++truth;
          ++bfs_itr;
        }
      }
    }
  }
}
