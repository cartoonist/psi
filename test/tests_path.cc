/**
 *    @file  tests_path.cc
 *   @brief  Test cases for Path submodule of module VarGraph.
 *
 *  Contains test cases for Path template class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@gmail.com>
 *
 *  @internal
 *       Created:  Tue Feb 27, 2018  18:20
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

#include "tests_base.h"
#include "vargraph.h"
#include "logger.h"


using namespace grem;

SCENARIO( "Basic test for a simple path in a variation graph", "[graph][path]" )
{
  std::vector< VarGraph::nodeid_type > nodes
    = { 20, 21, 23, 25, 26, 28, 29, 30, 32, 34, 35, 37 };
  std::vector< VarGraph::nodeid_type > nodes_shuff
    = { 29, 32, 34, 28, 21, 23, 26, 25, 37, 35, 30, 20 };
  std::vector< VarGraph::nodeid_type > other_nodes
    = { 56, 123, 9, 10, 27, 9, 10 };
  std::vector< VarGraph::nodeid_type > other_nodes_sorted
    = { 9, 10, 27, 56, 123 };
  std::vector< VarGraph::nodeid_type > invld_nodes = { 0 };
  std::vector< VarGraph::nodeid_type > empty;
  std::string nodes_str = "TGCTATGTGTAACTAGTAATGGTAATGGATATGTTGGGCTTTTTCCTTTGATTTATTTGA"
    "AGTAACGTTTGACAATCTATCACTAGGGGTAATGTGGGGAAGTGGAAAGAATACAAGAT";
  auto common_path_basic_test =
    [&nodes, &other_nodes, &invld_nodes, &empty, &nodes_shuff]
    ( const auto& path ) {
      typedef typename std::remove_reference< decltype( path ) >::type TPath;
      REQUIRE( length( path ) == nodes.size() );
      for ( const auto& n : nodes ) {
        REQUIRE( contains( path, n ) );
      }
      for ( const auto& on : other_nodes ) {
        REQUIRE( !contains( path, on ) );
      }
      REQUIRE( contains( path, nodes.begin(), nodes.end() ) );
      REQUIRE( !contains( path, other_nodes.begin(), other_nodes.end() ) );
      REQUIRE( !contains( path, empty.begin(), empty.end() ) );
      REQUIRE( !contains( path, invld_nodes.begin(), invld_nodes.end() ) );
      if ( std::is_same< Micro, typename TPath::spec_type >::value ) {
        REQUIRE( contains( path, nodes_shuff.begin(), nodes_shuff.end() ) );
      }
      else {
        REQUIRE( !contains( path, nodes_shuff.begin(), nodes_shuff.end() ) );
        REQUIRE( rcontains( path, nodes.rbegin(), nodes.rend() ) );
        REQUIRE( !rcontains( path, other_nodes.rbegin(), other_nodes.rend() ) );
        REQUIRE( !rcontains( path, empty.rbegin(), empty.rend() ) );
        REQUIRE( !rcontains( path, invld_nodes.rbegin(), invld_nodes.rend() ) );
      }
    };

  auto path_basic_test =
    [&common_path_basic_test, &nodes_str]
    ( const auto& path ) {
      common_path_basic_test( path );
      REQUIRE( sequence( path ) == nodes_str );
      REQUIRE( path.get_sequence_len() == 119 );
      REQUIRE( position_to_id( path, 0 ) == 20 );
      REQUIRE( position_to_offset( path, 0 ) == 0 );
      REQUIRE( position_to_id( path, 118 ) == 37 );
      REQUIRE( position_to_offset( path, 118 ) == 4 );
    };

  GIVEN( "A small variation graph" )
  {
    std::string vgpath = _testdir + "/data/small/x.xg";
    std::ifstream gifs( vgpath, std::ifstream::in | std::istream::binary );
    if ( !gifs ) {
      throw std::runtime_error( "cannot open file " + vgpath );
    }
    VarGraph vargraph( gifs );

    WHEN( "An empty path is initialized" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      initialize( path );
      initialize( dyn_path );
      initialize( cmp_path );
      initialize( hap_path );

      THEN( "It has no effect on its state" )
      {
        REQUIRE( !path.is_initialized() );
        REQUIRE( !dyn_path.is_initialized() );
        REQUIRE( !cmp_path.is_initialized() );
        REQUIRE( hap_path.is_initialized() );
      }
    }

    WHEN( "A Default path in the graph constructed at once" )
    {
      Path< VarGraph > path( &vargraph );
      path.set_nodes( nodes );
      initialize( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }

      WHEN( "It is saved to file and loaded again" )
      {
        std::string tmp_fpath = SEQAN_TEMP_FILENAME();
        save( path, tmp_fpath );
        clear( path );

        REQUIRE( length( path ) == 0 );
        REQUIRE( sequence( path ) == "" );
        REQUIRE( path.get_sequence_len() == 0 );
        REQUIRE( path.is_initialized() == false );

        open( path, tmp_fpath );

        THEN( "It should pass basic tests" )
        {
          path_basic_test( path );
        }
      }
    }

    WHEN( "A Default path in the graph constructed incrementally" )
    {
      Path< VarGraph > path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( path, n );
      }
      initialize( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Dynamic path in the graph constructed at once" )
    {
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      dyn_path.set_nodes( nodes.begin(), nodes.end() );
      initialize( dyn_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }

      WHEN( "It is saved to file and loaded again" )
      {
        std::string tmp_fpath = SEQAN_TEMP_FILENAME();
        save( dyn_path, tmp_fpath );
        clear( dyn_path );

        REQUIRE( length( dyn_path ) == 0 );
        REQUIRE( sequence( dyn_path ) == "" );
        REQUIRE( dyn_path.get_sequence_len() == 0 );
        REQUIRE( dyn_path.is_initialized() == false );

        open( dyn_path, tmp_fpath );

        THEN( "It should pass basic tests" )
        {
          path_basic_test( dyn_path );
        }
      }
    }

    WHEN( "A Dynamic path in the graph constructed incrementally" )
    {
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }
      initialize( dyn_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Compact path in the graph constructed at once" )
    {
      Path< VarGraph, Compact > cmp_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( cmp_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path );
      }

      WHEN( "It is saved to file and loaded again" )
      {
        std::string tmp_fpath = SEQAN_TEMP_FILENAME();
        save( cmp_path, tmp_fpath );
        clear( cmp_path );

        REQUIRE( length( cmp_path ) == 0 );
        REQUIRE( sequence( cmp_path ) == "" );
        REQUIRE( cmp_path.get_sequence_len() == 0 );
        REQUIRE( cmp_path.is_initialized() == false );

        open( cmp_path, tmp_fpath );

        THEN( "It should pass basic tests" )
        {
          path_basic_test( cmp_path );
        }
      }
    }

    WHEN( "A Micro path in the graph constructed at once" )
    {
      Path< VarGraph, Micro > mcr_path;
      mcr_path.set_nodes( nodes );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( mcr_path );
      }

      WHEN( "It is saved to file and loaded again" )
      {
        std::string tmp_fpath = SEQAN_TEMP_FILENAME();
        save( mcr_path, tmp_fpath );
        clear( mcr_path );

        REQUIRE( length( mcr_path ) == 0 );

        open( mcr_path, tmp_fpath );

        THEN( "It should pass basic tests" )
        {
          common_path_basic_test( mcr_path );
        }
      }
    }

    WHEN( "A Micro path in the graph constructed incrementally" )
    {
      Path< VarGraph, Micro > mcr_path;
      for ( const auto& n: nodes ) {
        add_node( mcr_path, n );
      }

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( mcr_path );
      }
    }

    WHEN( "A Haplotype path in the graph constructed at once" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      hap_path.set_nodes( nodes );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }

      WHEN( "It is saved to file and loaded again" )
      {
        std::string tmp_fpath = SEQAN_TEMP_FILENAME();
        save( hap_path, tmp_fpath );
        clear( hap_path );

        REQUIRE( length( hap_path ) == 0 );
        REQUIRE( hap_path.is_initialized() == true );

        open( hap_path, tmp_fpath );

        THEN( "It should pass basic tests" )
        {
          common_path_basic_test( hap_path );
        }
      }
    }

    WHEN( "A Haplotype path in the graph constructed incrementally" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( hap_path, n );
      }
      initialize( hap_path );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "An existing path in the graph reset by another set of nodes" )
    {
      Path< VarGraph > path( &vargraph );
      path.set_nodes( other_nodes );
      initialize( path );
      path.set_nodes( nodes );
      initialize( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "An existing Haplotype path in the graph reset by another set of nodes" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      hap_path.set_nodes( other_nodes_sorted );
      hap_path.set_nodes( nodes );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Dynamic path constructed by a Default path using assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      path.set_nodes( nodes );
      initialize( path );

      dyn_path = path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Dynamic path constructed by a Default path using move assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      path.set_nodes( nodes );
      initialize( path );

      dyn_path = std::move( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Dynamic path constructed by a Compact path using assignment" )
    {
      Path< VarGraph, Compact > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( path );

      dyn_path = path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Dynamic path constructed by a Compact path using move assignment" )
    {
      Path< VarGraph, Compact > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( path );

      dyn_path = std::move( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Default path constructed by a Dynamic path using assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }
      initialize( dyn_path );

      path = dyn_path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Default path constructed by a Dynamic path using move assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }
      initialize( dyn_path );

      path = std::move( dyn_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Default path constructed by a Compact path using assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( cmp_path );

      path = cmp_path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Default path constructed by a Compact path using move assignment" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( cmp_path );

      path = std::move( cmp_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
      }
    }

    WHEN( "A Compact path constructed by a Dynamic path using assignment" )
    {
      Path< VarGraph, Compact > cmp_path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }
      initialize( dyn_path );

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path = dyn_path;
#pragma GCC diagnostic pop

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path );
      }
    }

    WHEN( "A Compact path constructed by a Dynamic path using move assignment" )
    {
      Path< VarGraph, Compact > cmp_path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }
      initialize( dyn_path );

      cmp_path = std::move( dyn_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path );
      }
    }

    WHEN( "A Compact path constructed by a Default path using assignment" )
    {
      Path< VarGraph, Compact > cmp_path( &vargraph );
      Path< VarGraph > path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( path );

      cmp_path = path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path );
      }
    }

    WHEN( "A Compact path constructed by a Default path using move assignment" )
    {
      Path< VarGraph, Compact > cmp_path( &vargraph );
      Path< VarGraph > path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( path );

      cmp_path = std::move( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path );
      }
    }

    WHEN( "A Haplotype path constructed by a Default path using assignment" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph > path( &vargraph );
      path.set_nodes( nodes );
      initialize( path );

      hap_path = path;

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Haplotype path constructed by a Dynamic path using assignment" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( dyn_path, n );
      }
      initialize( dyn_path );

      hap_path = dyn_path;

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Haplotype path constructed by a Compact path using assignment" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( cmp_path );

      hap_path = cmp_path;

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Dynamic path constructed by a Dynamic path using assignment" )
    {
      Path< VarGraph, Dynamic > dyn_path2( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      dyn_path.set_nodes( nodes.begin(), nodes.end() );
      initialize( dyn_path );

      dyn_path2 = dyn_path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path2 );
      }
    }

    WHEN( "A Dynamic path constructed by a Dynamic path using move assignment" )
    {
      Path< VarGraph, Dynamic > dyn_path2( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      dyn_path.set_nodes( nodes.begin(), nodes.end() );
      initialize( dyn_path );

      dyn_path2 = std::move( dyn_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path2 );
      }
    }

    WHEN( "A Default path constructed by a Default path using assignment" )
    {
      Path< VarGraph > path2( &vargraph );
      Path< VarGraph > path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( path, n );
      }
      initialize( path );

      path2 = path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path2 );
      }
    }

    WHEN( "A Default path constructed by a Default path using move assignment" )
    {
      Path< VarGraph > path2( &vargraph );
      Path< VarGraph > path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( path, n );
      }
      initialize( path );

      path2 = std::move( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path2 );
      }
    }

    WHEN( "A Compact path constructed by a Compact path using assignment" )
    {
      Path< VarGraph, Compact > cmp_path2( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( cmp_path );

      cmp_path2 = cmp_path;

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path2 );
      }
    }

    WHEN( "A Compact path constructed by a Compact path using move assignment" )
    {
      Path< VarGraph, Compact > cmp_path2( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( nodes );
#pragma GCC diagnostic pop
      initialize( cmp_path );

      cmp_path2 = std::move( cmp_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( cmp_path2 );
      }
    }

    WHEN( "A Micro path constructed by a Micro path using assignment" )
    {
      Path< VarGraph, Micro > path2;
      Path< VarGraph, Micro > path;
      path.set_nodes( nodes );

      path2 = path;

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( path2 );
      }
    }

    WHEN( "A Micro path constructed by a Micro path using move assignment" )
    {
      Path< VarGraph, Micro > path2;
      Path< VarGraph, Micro > path;
      for ( const auto& n: nodes ) {
        add_node( path, n );
      }

      path2 = std::move( path );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( path2 );
      }
    }

    WHEN( "A Haplotype path constructed by a Haplotype path using assignment" )
    {
      Path< VarGraph, Haplotype > hap_path2( &vargraph );
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      for ( const auto& n: nodes ) {
        add_node( hap_path, n );
      }
      initialize( hap_path );

      hap_path2 = hap_path;

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path2 );
      }
    }

    WHEN( "A Default path is extended by another Default path" )
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
      initialize( path1 );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Haplotype path is extended by another Haplotype path" )
    {
      Path< VarGraph, Haplotype > hap_path1( &vargraph );
      Path< VarGraph, Haplotype > hap_path2( &vargraph );
      Path< VarGraph, Haplotype >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( hap_path1, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( hap_path2, nodes[ i ] );
      }
      initialize( hap_path2 );
      hap_path1 += hap_path2;
      initialize( hap_path1 );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path1 );
      }
    }

    WHEN( "A Haplotype path is extended by a Default path" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Haplotype >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( hap_path, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( path, nodes[ i ] );
      }
      hap_path += path;
      initialize( hap_path );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Haplotype path is extended by a Dynamic path" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      Path< VarGraph, Haplotype >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( hap_path, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( dyn_path, nodes[ i ] );
      }
      hap_path += dyn_path;
      initialize( hap_path );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Haplotype path is extended by a Compact path" )
    {
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph, Compact > cmp_path( &vargraph );
      Path< VarGraph, Haplotype >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( hap_path, nodes[ i ] );
      }
      std::vector< VarGraph::nodeid_type > subset;
      for ( ; i < nodes.size(); ++i ) {
        subset.push_back( nodes[i] );
      }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      cmp_path.set_nodes( subset );
#pragma GCC diagnostic pop
      hap_path += cmp_path;
      initialize( hap_path );

      THEN( "It should pass basic tests" )
      {
        common_path_basic_test( hap_path );
      }
    }

    WHEN( "A Dynamic path is extended by a Default path" )
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
      initialize( path1 );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Dynamic path is extended by a Haplotype path" )
    {
      Path< VarGraph, Dynamic > dyn_path( &vargraph );
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph, Dynamic >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( dyn_path, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( hap_path, nodes[ i ] );
      }
      initialize( hap_path );
      dyn_path += hap_path;
      initialize( dyn_path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( dyn_path );
      }
    }

    WHEN( "A Default path is extended by a Dynamic Path" )
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
      initialize( path1 );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Default path is extended by a Haplotype path" )
    {
      Path< VarGraph > path( &vargraph );
      Path< VarGraph, Haplotype > hap_path( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path, nodes[ i ] );
      }
      for ( ; i < nodes.size(); ++i ) {
        add_node( hap_path, nodes[ i ] );
      }
      initialize( hap_path );
      path += hap_path;
      initialize( path );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path );
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
      initialize( path1 );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Default path is extended by a Compact Path" )
    {
      Path< VarGraph > path1( &vargraph );
      Path< VarGraph, Compact > path2( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path1, nodes[ i ] );
      }
      std::vector< VarGraph::nodeid_type > subset;
      for ( ; i < nodes.size(); ++i ) {
        subset.push_back( nodes[i] );
      }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      path2.set_nodes( subset );
#pragma GCC diagnostic pop
      path1 += path2;
      initialize( path1 );

      THEN( "It should pass basic tests" )
      {
        path_basic_test( path1 );
      }
    }

    WHEN( "A Dynamic path is extended by a Compact Path" )
    {
      Path< VarGraph, Dynamic > path1( &vargraph );
      Path< VarGraph, Compact > path2( &vargraph );
      Path< VarGraph >::size_type i = 0;
      for ( ; i < nodes.size() - 3; ++i ) {
        add_node( path1, nodes[ i ] );
      }
      std::vector< VarGraph::nodeid_type > subset;
      for ( ; i < nodes.size(); ++i ) {
        subset.push_back( nodes[i] );
      }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
      path2.set_nodes( subset );
#pragma GCC diagnostic pop
      path1 += path2;
      initialize( path1 );

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
