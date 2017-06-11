/*
 * =====================================================================================
 *
 * Filename: tests-vargraph.cc
 *
 * Created: Mon Feb 13, 2017  22:24
 *
 * Description: Test VarGraph class.
 *
 * Copyright (c) 2017, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari (cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

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


// Test if the parameter is a constant reference.
template<typename T>
bool is_const_ref(const T&) { return true; }

// Test if the parameter is a non-const pointer.
template<typename T>
bool is_mutable_ptr(T*) { return true; }


SCENARIO( "Loading variation graph from a vg file", "[input]" )
{
  GIVEN( "A small graph without graph name" )
  {
    std::string vgpath = _testdir + "/data/small/x.vg";
    VarGraph vargraph(vgpath.c_str());

    REQUIRE( vargraph.nodes_size() == 210 );

    grem::VarGraph::NodeID node_id = 0;
    REQUIRE( !vargraph.has_node(node_id) );

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<unsigned int> distribution(0, vargraph.nodes_size()-1);
    unsigned int random_idx = distribution(generator);

    auto& a_node = vargraph.node_at(random_idx);
    REQUIRE( vargraph.has_node(a_node.id()) );
    REQUIRE( vargraph.has_node(&a_node) );
    REQUIRE( vargraph.node_at(random_idx).id() == random_idx+1 );
    REQUIRE( vargraph.node_by(a_node.id()).sequence() == a_node.sequence() );

    REQUIRE( is_const_ref(vargraph.node_at(random_idx)) );
    REQUIRE( is_mutable_ptr(vargraph.mutable_node_at(random_idx)) );
    REQUIRE( is_const_ref(vargraph.node_by(a_node.id())) );
    REQUIRE( is_mutable_ptr(vargraph.mutable_node_by(a_node.id())) );

    REQUIRE( vargraph.edges_size() == 291 );

    random_idx = distribution(generator);

    auto& an_edge = vargraph.edge_at(random_idx);
  }
}

// VarGraph graph iterators test scenarios.
SCENARIO ( "Get unique haplotype using Haplotyper graph iterator", "[graph][iterator]" )
{
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.vg";
    VarGraph vargraph(vgpath.c_str());

    WHEN ( "the three haplotypes are generated using Haplotyper" )
    {
      seqan::Iterator < VarGraph, Haplotyper<> >::Type hap_itr (vargraph);

      std::vector < VarGraph::NodeID > haplotype1;
      std::vector < VarGraph::NodeID > haplotype2;
      std::vector < VarGraph::NodeID > haplotype3;
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
        REQUIRE ( haplotype2.size() == 134 );
        REQUIRE ( haplotype3.size() > 2 );      // randomized path
      }
    }
  }
}

SCENARIO ( "Traverse a variation graph using backtracking algorithm", "[graph][iterator]" )
{
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.vg";
    VarGraph vargraph(vgpath.c_str());

    unsigned int kmer_len = 20;

    WHEN ( "enumerating all kmers of length " + std::to_string ( kmer_len ) )
    {
      std::string truth_dir = _testdir + "/data/small/";
      std::string truth_filepath = truth_dir + std::to_string ( kmer_len ) + "-mers";
      std::ifstream truth_stream ( truth_filepath, std::ifstream::in );
      std::string true_kmer;
      VarGraph::NodeID true_snode_id;
      unsigned int true_offset;

      seqan::Iterator< VarGraph, Backtracker<> >::Type bt_itr ( vargraph );
      std::vector< VarGraph::NodeID > trav_path;
      std::string trav_seq = "";
      // :TODO:Mon May 22 11:16:\@cartoonist: add REQUIREs and assertions.
      // :TODO:Mon May 22 11:15:\@cartoonist: add a method to simulate a kmer.
      for ( unsigned int n_idx = 0; n_idx < vargraph.nodes_size(); ++n_idx ) {
        const VarGraph::Node & start_node = vargraph.node_at ( n_idx );
        VarGraph::NodeID start_node_id = start_node.id();
        unsigned int label_len = start_node.sequence().length();

        for ( unsigned int offset = 0; offset < label_len; ++offset ) {
          go_begin ( bt_itr, start_node_id );

          while ( !at_end( bt_itr ) ) {
            while ( !at_end( bt_itr ) ) {
              trav_path.push_back ( *bt_itr );
              if ( *bt_itr != start_node_id ) {
                trav_seq += vargraph.node_by ( *bt_itr ).sequence();
              }
              else {
                trav_seq = vargraph.node_by ( *bt_itr ).sequence().substr( offset );
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
            VarGraph::NodeID poped_id = 0;
            while ( !trav_path.empty() && poped_id != *bt_itr ) {
              poped_id = trav_path.back();
              trav_len -= vargraph.node_by ( poped_id ).sequence().length();
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
