/**
 *    @file  test_hbitvector.cpp
 *   @brief  Test scenarios for `HBitVector` template class.
 *
 *  This test module contains all test scenarios for `HBitVector` class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Jul 24, 2023  14:49
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

//#include<cxxabi.h>

#include <sstream>

#include <psi/hbitvector.hpp>
#include <psi/graph.hpp>
#include <psi/crs_matrix.hpp>
#include <gum/gfa_utils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include "test_base.hpp"


using namespace psi;

TEMPLATE_SCENARIO_SIG(
    "L1 begin position in bitvectors", "[hbitvector]",
    ( ( typename T, int L ), T, L ),
    ( HBitVector< 128, uint32_t, uint64_t >, 11431 ),
    ( HBitVector< 128, uint32_t, uint64_t >, 4096 ),
    ( HBitVector< 128, uint32_t, uint64_t >, 128 ),
    ( HBitVector< 128, uint32_t, uint64_t >, 64 ),
    ( HBitVector< 128, uint32_t, uint64_t >, 38 ),
    ( HBitVector< 256, uint32_t, uint64_t >, 11431 ),
    ( HBitVector< 256, uint32_t, uint64_t >, 4096 ),
    ( HBitVector< 256, uint32_t, uint64_t >, 256 ),
    ( HBitVector< 256, uint32_t, uint64_t >, 128 ),
    ( HBitVector< 256, uint32_t, uint64_t >, 73 ),
    ( HBitVector< 1024, uint32_t, uint64_t >, 11431 ),
    ( HBitVector< 1024, uint32_t, uint64_t >, 4096 ),
    ( HBitVector< 1024, uint32_t, uint64_t >, 1024 ),
    ( HBitVector< 1024, uint32_t, uint64_t >, 64 ),
    ( HBitVector< 2048, uint32_t, uint64_t >, 11431 ),
    ( HBitVector< 2048, uint32_t, uint64_t >, 4096 ),
    ( HBitVector< 2048, uint32_t, uint64_t >, 2048 ),
    ( HBitVector< 2048, uint32_t, uint64_t >, 256 ),
    ( HBitVector< 128, uint32_t, uint32_t >, 11431 ),
    ( HBitVector< 128, uint32_t, uint32_t >, 4096 ),
    ( HBitVector< 128, uint32_t, uint32_t >, 128 ),
    ( HBitVector< 128, uint32_t, uint32_t >, 64 ),
    ( HBitVector< 128, uint32_t, uint32_t >, 38 ),
    ( HBitVector< 256, uint32_t, uint32_t >, 11431 ),
    ( HBitVector< 256, uint32_t, uint32_t >, 4096 ),
    ( HBitVector< 256, uint32_t, uint32_t >, 256 ),
    ( HBitVector< 256, uint32_t, uint32_t >, 128 ),
    ( HBitVector< 256, uint32_t, uint32_t >, 73 ),
    ( HBitVector< 1024, uint32_t, uint32_t >, 11431 ),
    ( HBitVector< 1024, uint32_t, uint32_t >, 4096 ),
    ( HBitVector< 1024, uint32_t, uint32_t >, 1024 ),
    ( HBitVector< 1024, uint32_t, uint32_t >, 64 ),
    ( HBitVector< 2048, uint32_t, uint32_t >, 11431 ),
    ( HBitVector< 2048, uint32_t, uint32_t >, 4096 ),
    ( HBitVector< 2048, uint32_t, uint32_t >, 2048 ),
    ( HBitVector< 2048, uint32_t, uint32_t >, 256 ) )
{
  using hbv_type = T;
  using bitset_type = typename hbv_type::bitset_type;
  using execution_space_type = typename hbv_type::execution_space;
  using policy_type = Kokkos::TeamPolicy< execution_space_type >;
  using member_type = typename policy_type::member_type;

  std::size_t len = L;
  auto width = gum::widthof< bitset_type >::value;
  auto nof_bitsets = ( L/width + ((bool)(L%width)) );
  auto bitset_size = sizeof( bitset_type );

  GIVEN( "A Kokkos Team Execution Policy consisting of "
         + std::to_string( len ) + " teams" )
  {
    auto policy = policy_type( len, Kokkos::AUTO );

    WHEN( "Setting scratch size" ) {
      policy = hbv_type::set_scratch_size( policy, len );

      THEN( "L1 scratch size should be equal to the given template parameter" ) {
        REQUIRE( policy.scratch_size( 0 ) == hbv_type::L1_SIZE_BYTES );
        REQUIRE( policy.scratch_size( 0 ) == hbv_type::l1_scratch_size() );
      }

      THEN( "L2 scratch size should be equal to the rest of bitvector" ) {
        std::ptrdiff_t l2_scratch_size =
            nof_bitsets * bitset_size - hbv_type::L1_SIZE_BYTES;
        std::size_t actual_l2_scratch_size = std::max( l2_scratch_size, 0L );
        REQUIRE( policy.scratch_size( 1 ) == actual_l2_scratch_size );
        REQUIRE( policy.scratch_size( 1 ) == hbv_type::l2_scratch_size( len ) );
        auto total_scratch_size = policy.scratch_size( 0 ) + policy.scratch_size( 1 );
        REQUIRE( total_scratch_size == hbv_type::capacity( len ) );
      }
    }
  }

  GIVEN( "A hierarchical bitvector of length " + std::to_string( len )
         + " and width " + std::to_string( width ) )
  {
    Kokkos::View< uint64_t* > true_begins( "t", len );
    auto h_tb = Kokkos::create_mirror_view( true_begins );

    {
      // enumerating ground truth values for begin positions
      std::size_t tbegin = 0;
      // When the l1-centre is in the first `buffer` bits, the begin position is 0.
      std::size_t buffer = hbv_type::L1_SIZE / 2;
      // When the l1-centre passes `r_centre`, the begin position does not follow it anymore.
      std::size_t r_centre = ( hbv_type::num_bitsets( len ) - hbv_type::L1_NUM_BITSETS ) * width + buffer;
      for ( std::size_t i = 0; i < len; ++i ) {
        // When on the left side of `r_centre`, `tbegin` follows the centre.
        if ( i < r_centre ) {
          if ( buffer == 0 ) {  // when initial/width `buffer` bits are consumed
            ++tbegin;           // increment `tbegin`
            buffer = width;     // set buffer to bitset width once the first buffer value exhausted
          }
          --buffer;             // consume buffer value when `i < r_centre`
        }
        h_tb( i ) = tbegin * width;
      }
    }

    Kokkos::deep_copy( true_begins, h_tb );

    WHEN( "Initialising the hierarchical bitvector inside a Kokkos kernel" )
    {
      auto policy = policy_type( len, Kokkos::AUTO );
      policy = hbv_type::set_scratch_size( policy, len );

      Kokkos::View< unsigned char* > flags( "flags", len );
      Kokkos::parallel_for(
          "psi::test_hbitvector::l1_begin", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            hbv_type h_bv( len, row, tm );
            Kokkos::single( Kokkos::PerTeam( tm ), [ & ]() {
              if ( h_bv.l1_begin == true_begins( row ) ) flags( row ) = 1;
              /*
              // NOTE: uncomment `#include <cxxabi.h>`
              else {  // For debugging on CPU
                std::stringstream ss;
                auto hbv_type_name = abi::__cxa_demangle(
                    typeid( hbv_type ).name(), NULL, NULL, NULL );
                ss << "=== " << hbv_type_name << " ===\n";
                std::free( hbv_type_name );
                ss << "Error: " << row << ": " << h_bv.l1_begin
                   << " != " << true_begins( row ) << "\n";
                std::cout << ss.str() << std::endl;
                assert( false );
              }
              */
            } );
          } );

      THEN( "The begin position of the bitvector for each team should be "
            "64-bit aligned and centred around the given row" )
      {
        std::size_t all_set = 0;
        Kokkos::parallel_reduce(
            "psi::test_hbit_vector::l1_begin_assess", len,
            KOKKOS_LAMBDA( const uint64_t i, std::size_t& all_set_local ) {
              if ( flags( i ) == 1 ) all_set_local += 1;
            },
            all_set );

        REQUIRE( all_set == len );
      }
    }
  }
}

SCENARIO( "Unorganised scenario", "[temp]" )
{
  typedef int scalar_t;
  typedef Kokkos::DefaultHostExecutionSpace host_space;
  typedef KokkosSparse::CrsMatrix< scalar_t, int32_t, host_space > xcrsmatrix_t;
  typedef psi::CRSMatrix< psi::crs_matrix::RangeDynamic, bool, uint32_t, uint64_t > range_crsmatrix_t;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  std::string graph_path = test_data_dir + "/small/x.gfa";
  graph_type graph;
  gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );

  auto a = psi::util::adjacency_matrix< xcrsmatrix_t >( graph );
  range_crsmatrix_t ra( a );

  auto h_a_entries = ra.entries_view();
  auto h_a_rowmap = ra.rowmap_view();

  auto a_entries = Kokkos::create_mirror_view_and_copy( Kokkos::DefaultExecutionSpace{}, h_a_entries );
  auto a_rowmap = Kokkos::create_mirror_view_and_copy( Kokkos::DefaultExecutionSpace{}, h_a_rowmap );

  auto ch_a_entries = Kokkos::create_mirror_view( a_entries );
  auto ch_a_rowmap = Kokkos::create_mirror_view( a_rowmap );

  for ( std::size_t i = 0; i < a_entries.extent( 0 ); ++i )
    REQUIRE( h_a_entries( i ) == ch_a_entries( i ) );
  for ( std::size_t i = 0; i < a_rowmap.extent( 0 ); ++i )
    REQUIRE( h_a_rowmap( i ) == ch_a_rowmap( i ) );
}
