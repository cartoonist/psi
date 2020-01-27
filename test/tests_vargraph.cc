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

SCENARIO( "Basic test for a path in a variation graph", "[graph][path]" )
{
  std::vector< VarGraph::nodeid_type > nodes
    = { 20, 21, 23, 25, 26, 28, 29, 30, 32, 34, 35, 37 };
  std::vector< VarGraph::nodeid_type > other_nodes
    = { 56, 123, 9, 10, 27, 9, 10 };
  auto path_basic_test = [&nodes, &other_nodes]( const auto& path ){
    REQUIRE( length( path ) == nodes.size() );
    for ( const auto& n : nodes ) {
      REQUIRE( contains( path, n ) );
    }
    for ( const auto& on : other_nodes ) {
      REQUIRE( !contains( path, on ) );
    }
    REQUIRE( contains( path, nodes.begin(), nodes.end() ));
    REQUIRE( !contains( path, other_nodes.begin(), other_nodes.end() ) );
  };

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::istream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    WHEN( "A path in the graph constructed at once" )
    {
      Path< VarGraph > path( &vargraph );
      path.set_nodes( nodes );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A path in the graph constructed incrementally" )
    {
      Path< VarGraph > path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( path, n );
      }

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Dynamic path in the graph constructed incrementally" )
    {
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Dynamic path constructed by a Normal path using assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );

      path.set_nodes( nodes );
      dyn_path = path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Normal path constructed by a Dynamic path using assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }

      path = dyn_path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Dynamic path constructed by a Dynamic path using assignment" )
    {
      Path< VarGraph, Dynamic > dyn_path2( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }

      dyn_path2 = dyn_path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path2 );
      }
    }

    WHEN( "A path is extended by another path" )
    {
      Path< VarGraph > path1( &vargraph );
      Path< VarGraph > path2( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path1, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( path2, nodes[ i ] );
      }
      path1 += path2;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Dynamic path is extended by another path" )
    {
      Path< VarGraph, Dynamic > path1( &vargraph );
      Path< VarGraph > path2( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path1, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( path2, nodes[ i ] );
      }
      path1 += path2;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A path is extended by a Dynamic Path" )
    {
      Path< VarGraph > path1( &vargraph );
      Path< VarGraph, Dynamic > path2( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path1, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( path2, nodes[ i ] );
      }
      path1 += path2;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Dynamic path is extended by another Dynamic path" )
    {
      Path< VarGraph, Dynamic > path1( &vargraph );
      Path< VarGraph, Dynamic > path2( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path1, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( path2, nodes[ i ] );
      }
      path1 += path2;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }
  }
}

SCENARIO( "Trim a path in a variation graph", "[graph][path]" )
{
  std::vector< VarGraph::nodeid_type > nodes
    = { 20, 21, 23, 25, 26, 28, 29, 30, 32, 34, 35, 37 };
  std::string init_sequence = "TGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTCCTTTGATTTA"
    "TTTGAAGTAACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAGTGGAAAGAATACAAGAT";

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::istream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    GIVEN( "A path in the graph" )
    {
      Path< VarGraph > path( &vargraph );
      //path.set_nodes( nodes );
      for ( const auto& n : nodes ) {
        add_node( path, n );
      }

      REQUIRE( path.get_sequence() == init_sequence );

      WHEN( "The forward sequence with non-zero context is computed" )
      {
        std::string s = sequence( path, 11 );

        THEN( "The sequence of the first and last nodes should be trimmed according to context" )
        {
          REQUIRE( s == init_sequence.substr( 31 ) );
        }
      }

      WHEN( "The forward sequence with non-zero context is computed" )
      {
        auto other = path;
        trim_back( other, 37 );
        std::string s = sequence( other, 11 );

        THEN( "The sequence of the first and last nodes should be trimmed according to context" )
        {
          REQUIRE( s == init_sequence.substr( 31, 81 ) );
        }
      }

      WHEN( "The reversed sequence is computed" )
      {
        std::string s = sequence( path, Reversed() );

        THEN( "The reversed sequence should be returned" )
        {
          std::string rev_s;
          std::copy( init_sequence.rbegin(),
              init_sequence.rend(),
              std::back_inserter( rev_s ) );
          REQUIRE( s == rev_s );
        }
      }

      WHEN( "The reversed sequence with non-zero context is computed" )
      {
        std::string s = sequence( path, Reversed(), 11 );

        THEN( "The sequence of the first and last nodes should be trimmed according to context" )
        {
          std::string truth = init_sequence.substr( 31 );
          std::string rev_s;
          std::copy( truth.rbegin(), truth.rend(), std::back_inserter( rev_s ) );
          REQUIRE( s == rev_s );
        }
      }

      WHEN( "The last node is trimmed" )
      {
        VarGraph::offset_type trimmed_len
          = path.get_sequence().length() - vargraph.node_length( path.get_nodes().back() );
        trim_back( path, 37 );

        THEN( "Its length and sequence should be decreased accordingly" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( 0, trimmed_len ) );
        }
      }

      WHEN( "Trim further" )
      {
        std::size_t trim_len = 0;
        auto it = path.get_nodes().end();
        for ( it -= 6; it != path.get_nodes().end(); ++it ) {
          trim_len += vargraph.node_length( *it );
        }
        VarGraph::offset_type trimmed_len = path.get_sequence().length() - trim_len;
        trim_back( path, 29 );
        THEN( "Its length and sequence should be decreased accordingly" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( 0, trimmed_len ) );
        }
      }

      WHEN( "Trim by providing zero as node ID" )
      {
        VarGraph::offset_type trimmed_len
          = path.get_sequence().length() - vargraph.node_length( path.get_nodes().back() );
        trim_back( path, 0 );
        THEN( "The last node should be trimmed" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( 0, trimmed_len ) );
        }
      }

      WHEN( "Trim by providing no parameter" )
      {
        VarGraph::offset_type trimmed_len
          = path.get_sequence().length() - vargraph.node_length( path.get_nodes().back() );
        trim_back( path );
        THEN( "The last node should be trimmed" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( 0, trimmed_len ) );
        }
      }

      WHEN( "Trim by providing unavailable node ID" )
      {
        trim_back( path, 70 );
        THEN( "It should be empty" )
        {
          REQUIRE( path.get_sequence().length() == 0 );
        }
      }
    }

    GIVEN( "A Dynamic Path in the graph" )
    {
      Path< VarGraph, Dynamic > path( &vargraph );
      for ( const auto& n : nodes ) {
        add_node( path, n );
      }

      REQUIRE( path.get_sequence() == init_sequence );

      WHEN( "The first node is trimmed" )
      {
        VarGraph::offset_type trim_len =
          vargraph.node_length( path.get_nodes().front() );
        VarGraph::offset_type trimmed_len = path.get_sequence().length() - trim_len;
        trim_front( path, 20 );

        THEN( "Its length and sequence should be decreased accordingly" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( trim_len ) );
        }
      }

      WHEN( "Trim further" )
      {
        std::size_t trim_len = 0;
        for ( auto it = path.get_nodes().begin(); it != path.get_nodes().end() - 8; ++it ) {
          trim_len += vargraph.node_length( *it );
        }
        VarGraph::offset_type trimmed_len = path.get_sequence().length() - trim_len;
        trim_front( path, 25 );
        THEN( "Its length and sequence should be decreased accordingly" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( trim_len ) );
        }
      }

      WHEN( "Trim by providing zero as node ID" )
      {
        VarGraph::offset_type trim_len =
          vargraph.node_length( path.get_nodes().front() );
        VarGraph::offset_type trimmed_len = path.get_sequence().length() - trim_len;
        trim_front( path, 0 );
        THEN( "The last node should be trimmed" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( trim_len ) );
        }
      }

      WHEN( "Trim by providing no parameter" )
      {
        VarGraph::offset_type trim_len =
          vargraph.node_length( path.get_nodes().front() );
        VarGraph::offset_type trimmed_len = path.get_sequence().length() - trim_len;
        trim_front( path );
        THEN( "The last node should be trimmed" )
        {
          REQUIRE( path.get_sequence().length() == trimmed_len );
          REQUIRE( path.get_sequence() == init_sequence.substr( trim_len ) );
        }
      }

      WHEN( "Trim by providing unavailable node ID" )
      {
        trim_front( path, 70 );
        THEN( "It should be empty" )
        {
          REQUIRE( path.get_sequence().length() == 0 );
        }
      }
    }
  }
}

SCENARIO( "Trim a path to the length of k", "[graph][path]" )
{
  GIVEN ( "A small variation graph and two paths: one Default and one Dynamic" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );
    Path< VarGraph > path( &vargraph );
    Path< VarGraph, Dynamic > dyn_path( &vargraph );
    path.set_nodes( { 2, 5, 6, 7, 9, 11, 12 } );
    dyn_path = path;
    unsigned int k = 5;

    WHEN( "Trim-back a path to the length of " + std::to_string( k ) )
    {
      trim_back_by_len( path, k );
      initialize( path );

      THEN( "It should be trimmed to the length of " + std::to_string( k ) )
      {
        REQUIRE( path.get_sequence_len() == 5 );
        REQUIRE( position_to_id( path, 0 ) == 2 );
        REQUIRE( position_to_offset( path, 0 ) == 0 );
      }
    }

    WHEN( "Trim-front a path to the length of " + std::to_string( k ) )
    {
      trim_front_by_len( dyn_path, k );
      initialize( dyn_path );

      THEN( "It should be moved one node forward preserving the length of at least " + std::to_string( k ) )
      {
        REQUIRE( dyn_path.get_sequence_len() == 5 );
        REQUIRE( position_to_id( dyn_path, 4 ) == 12 );
        REQUIRE( position_to_offset( dyn_path, 4 ) == 3 );
      }
    }
  }
}

SCENARIO( "Query node coordinates by position in the path", "[graph][path]" )
{
  std::vector< VarGraph::nodeid_type > nodes
    = { 20, 21, 23, 25, 26, 28, 29, 30, 32, 34, 35, 37 };
  std::string init_sequence = "TGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTCCTTTGATTTA"
    "TTTGAAGTAACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAGTGGAAAGAATACAAGAT";

  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::istream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    WHEN( "A path in the graph" )
    {
      Path< VarGraph > path( &vargraph );
      path.set_nodes( nodes );
      initialize( path );

      THEN( "Mapping positions to node coordinates should be correct" )
      {
        REQUIRE( position_to_id( path, 0 ) == 20 );
        REQUIRE( position_to_offset( path, 0 ) == 0 );
        REQUIRE( position_to_id( path, 18 ) == 20 );
        REQUIRE( position_to_offset( path, 18 ) == 18 );
        REQUIRE( position_to_id( path, 40 ) == 20 );
        REQUIRE( position_to_offset( path, 40 ) == 40 );
        REQUIRE( position_to_id( path, 41 ) == 21 );
        REQUIRE( position_to_offset( path, 41 ) == 0 );
        REQUIRE( position_to_id( path, 42 ) == 23 );
        REQUIRE( position_to_offset( path, 42 ) == 0 );
        REQUIRE( position_to_id( path, 43 ) == 23 );
        REQUIRE( position_to_offset( path, 43 ) == 1 );
        REQUIRE( position_to_id( path, 44 ) == 25 );
        REQUIRE( position_to_offset( path, 44 ) == 0 );
        REQUIRE( position_to_id( path, 100 ) == 32 );
        REQUIRE( position_to_offset( path, 100 ) == 16 );
        REQUIRE( position_to_id( path, 113 ) == 35 );
        REQUIRE( position_to_offset( path, 113 ) == 11 );
        REQUIRE( position_to_id( path, 116 ) == 37 );
        REQUIRE( position_to_offset( path, 116 ) == 2 );
        REQUIRE( position_to_id( path, 118 ) == 37 );
        REQUIRE( position_to_offset( path, 118 ) == 4 );
        REQUIRE_THROWS( position_to_id( path, 119 ) );
      }
    }
  }
}

// VarGraph graph iterators test scenarios.
SCENARIO ( "Get unique full haplotype using Haplotyper graph iterator", "[graph][iterator][haplotyper]" )
{
  GIVEN ( "A tiny variation graph" )
  {
    std::string vgpath = _testdir + "/data/tiny/tiny.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );

    WHEN ( "the eigth haplotypes are generated using Haplotyper" )
    {
      seqan::Iterator < VarGraph, Haplotyper >::Type hap_itr (vargraph);

      Path< VarGraph > haplotype1( &vargraph );
      Path< VarGraph > haplotype2( &vargraph );
      Path< VarGraph > haplotype3( &vargraph );
      Path< VarGraph > haplotype4( &vargraph );
      Path< VarGraph > haplotype5( &vargraph );
      Path< VarGraph > haplotype6( &vargraph );
      Path< VarGraph > haplotype7( &vargraph );
      Path< VarGraph > haplotype8( &vargraph );
      get_uniq_full_haplotype( haplotype1, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype2, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype3, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype4, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype5, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype6, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype7, hap_itr, 1 );
      get_uniq_full_haplotype( haplotype8, hap_itr, 1 );

      THEN ( "they should be unique" )
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
        REQUIRE( level( hap_itr ) == 8 );
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

      Path< VarGraph > haplotype1( &vargraph );
      Path< VarGraph > haplotype2( &vargraph );
      Path< VarGraph > haplotype3( &vargraph );
      get_uniq_full_haplotype( haplotype1, hap_itr );
      get_uniq_full_haplotype( haplotype2, hap_itr );
      get_uniq_full_haplotype( haplotype3, hap_itr );

      THEN ( "they should be unique" )
      {
        std::string hapstr1 = sequence( haplotype1 );
        std::string hapstr2 = sequence( haplotype2 );
        std::string hapstr3 = sequence( haplotype3 );
        REQUIRE( hapstr1 != hapstr2 );
        REQUIRE( hapstr2 != hapstr3 );
        REQUIRE( hapstr1 != hapstr3 );
      }
      AND_THEN ( "they should have the correct length" )
      {
        REQUIRE( length( haplotype1 ) == 147 );
        REQUIRE( length( haplotype2 ) > 130 );  // randomised path.
        REQUIRE( length( haplotype3 ) > 130 );  // randomised path.
      }
      AND_THEN ( "they all should cover 'merge' nodes" )
      {
        std::vector< Path< VarGraph > > paths_set;
        paths_set.push_back( haplotype1 );
        paths_set.push_back( haplotype2 );
        paths_set.push_back( haplotype3 );

        REQUIRE( get_path_coverage( haplotype1.get_nodes().front(), paths_set ) == 3 );
        REQUIRE( get_path_coverage( 2, paths_set ) != 3 );
        REQUIRE( get_path_coverage( 6, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 9, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 18, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 20, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 210, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 207, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 205, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 202, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 200, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 96, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 99, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 101, paths_set ) == 3 );
        REQUIRE( get_path_coverage( 104, paths_set ) == 3 );
        REQUIRE( get_path_coverage( haplotype1.get_nodes().back(), paths_set ) == 3 );
      }
      AND_THEN( "level of iterator should be the number of haplotypes" )
      {
        REQUIRE( level( hap_itr ) == 3 );
      }
    }
  }
}

SCENARIO( "A Haplotyper graph iterator raise on end", "[graph][iterator][haplotyper]" )
{
  GIVEN( "A small variation graph and a Haplotyper iterator with `raise_on_end` enabled" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );
    seqan::Iterator< VarGraph, Haplotyper >::Type hap_itr( vargraph );
    hap_itr.raise_on_end = true;

    WHEN( "A Haplotyper iterator reaches at end" )
    {
      while ( !at_end( hap_itr ) ) ++hap_itr;

      THEN( "It raise an exception if it is incremented" )
      {
        REQUIRE_THROWS( *hap_itr );
        REQUIRE_THROWS( ++hap_itr );
      }
    }
  }
}

SCENARIO( "Extend a path to length k using Haplotyper graph iterator", "[graph][iterator][haplotyper]" )
{
  GIVEN ( "A small variation graph and a Haplotyper graph iterator" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );
    seqan::Iterator< VarGraph, Haplotyper >::Type hap_itr( vargraph );
    unsigned int k = 5;

    WHEN( "A path is extend to length " + std::to_string( k ) )
    {
      Path< VarGraph > path( &vargraph );
      extend_to_k( path, hap_itr, k );
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
      Path< VarGraph > path( &vargraph );
      extend_to_k( path, hap_itr, k );
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
  GIVEN ( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );
    unsigned int context_len = 10;

    WHEN( "Generate 16x patched haplotypes are generated using a Haplotyper iterator" )
    {
      seqan::Iterator< VarGraph, Haplotyper >::Type hap_itr( vargraph );
      std::vector< Path< VarGraph > > pathset;

      for ( auto i = 0; i < 32; ++i )
        get_uniq_patched_haplotype( pathset, hap_itr, context_len );

      for ( std::size_t i = 0; i < pathset.size(); ++i )
        initialize( pathset.at( i ) );

      THEN( "The number of patches should be in correct range" )
      {
        REQUIRE( pathset.size() >= 49 );
        REQUIRE( pathset.size() <= 60 );
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
      seqan::Iterator< VarGraph, BFS >::Type bfs_itr (vargraph);

      THEN ( "nodes should be traversed in BFS order" )
      {
        VarGraph::nodeid_type truth = 1;
        while ( !at_end( bfs_itr ) ) {
          REQUIRE ( *bfs_itr == truth );  // The graph is such that its BFS is in order.
          ++truth;
          ++bfs_itr;
        }
        REQUIRE( truth == 211 );
      }
    }
  }
  GIVEN ( "A variation graph with more than one connected component" )
  {
    std::string vgpath = _testdir + "/data/multi/multi.xg";
    std::ifstream ifs( vgpath, std::ifstream::in | std::ifstream::binary );
    VarGraph vargraph( ifs );

    WHEN ( "traverse the graph using BFS graph iterator" )
    {
      seqan::Iterator< VarGraph, BFS >::Type bfs_itr (vargraph);

      THEN ( "nodes should be traversed in BFS order" )
      {
        VarGraph::nodeid_type truth = 1;
        while ( !at_end( bfs_itr ) ) {
          REQUIRE ( *bfs_itr == truth );  // The graph is such that its BFS is in order.
          ++truth;
          ++bfs_itr;
        }
        REQUIRE( truth == 226 );
      }
    }
  }
}
