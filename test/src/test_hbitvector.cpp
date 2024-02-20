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

#include <psi/hbitvector.hpp>
#include <psi/graph.hpp>
#include <psi/crs_matrix.hpp>
#include <psi/utils.hpp>
#include <gum/gfa_utils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

#include "test_base.hpp"


using namespace psi;

TEMPLATE_SCENARIO_SIG(
    "L1 begin position in bitvectors", "[hbitvector]",
    ( ( typename T, int L ), T, L ),
    ( HBitVector< 128 >, 11431 ),
    ( HBitVector< 128 >, 4096 ),
    ( HBitVector< 128 >, 128 ),
    ( HBitVector< 128 >, 64 ),
    ( HBitVector< 128 >, 38 ),
    ( HBitVector< 256 >, 11431 ),
    ( HBitVector< 256 >, 4096 ),
    ( HBitVector< 256 >, 256 ),
    ( HBitVector< 256 >, 128 ),
    ( HBitVector< 256 >, 73 ),
    ( HBitVector< 1024 >, 11431 ),
    ( HBitVector< 1024 >, 4096 ),
    ( HBitVector< 1024 >, 1024 ),
    ( HBitVector< 1024 >, 64 ),
    ( HBitVector< 2048 >, 11431 ),
    ( HBitVector< 2048 >, 4096 ),
    ( HBitVector< 2048 >, 2048 ),
    ( HBitVector< 2048 >, 256 ) )
{
  using hbv_type = T;
  using bitset_type = typename hbv_type::bitset_type;
  using policy_type = typename hbv_type::policy_type;
  using member_type = typename hbv_type::member_type;
  using size_type = typename hbv_type::size_type;

  std::size_t len = L;
  auto width = gum::widthof< bitset_type >::value;
  auto bitset_size = sizeof( bitset_type );

  size_type nof_bitsets = ( L / width + ( ( bool )( L % width ) ) );
  nof_bitsets = Kokkos::max( hbv_type::L1_NUM_BITSETS, nof_bitsets );

  GIVEN( "A Kokkos Team Execution Policy consisting of "
         + std::to_string( len ) + " teams" )
  {
    auto policy = policy_type( len, Kokkos::AUTO );

    WHEN( "Setting scratch size" ) {
      hbv_type::set_scratch_size( policy, len );

      THEN( "L1 scratch size should be equal to the given template parameter" ) {
        REQUIRE( policy.scratch_size( 0 ) == hbv_type::L1_SIZE_BYTES );
        REQUIRE( policy.scratch_size( 0 ) == hbv_type::l1_scratch_size() );
      }

      THEN( "L2 scratch size should be equal to the rest of bitvector" ) {
        std::ptrdiff_t l2_scratch_size =
            nof_bitsets * bitset_size - hbv_type::L1_SIZE_BYTES;
        std::size_t actual_l2_scratch_size = std::max( l2_scratch_size, 0L );
        REQUIRE( nof_bitsets == hbv_type::num_bitsets( len ) );
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
      hbv_type::set_scratch_size( policy, len );

      Kokkos::View< unsigned char* > flags( "flags", len );
      Kokkos::parallel_for(
          "psi::test_hbitvector::l1_begin", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            hbv_type hbv( tm, len, row );
            Kokkos::single( Kokkos::PerTeam( tm ), [=]() {
              if ( hbv.l1_begin == true_begins( row ) ) flags( row ) = 1;
              /*
              // NOTE: uncomment `#include <cxxabi.h>`
              else {  // For debugging on CPU
                std::stringstream ss;
                auto hbv_type_name = abi::__cxa_demangle(
                    typeid( hbv_type ).name(), NULL, NULL, NULL );
                ss << "=== " << hbv_type_name << " ===\n";
                std::free( hbv_type_name );
                ss << "Error: " << row << ": " << hbv.l1_begin
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
            KOKKOS_LAMBDA ( const uint64_t i, std::size_t& all_set_local ) {
              if ( flags( i ) == 1 ) all_set_local += 1;
            },
            all_set );

        REQUIRE( all_set == len );
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG(
    "Range set operation in hbitvector", "[hbitvector]",
    ( ( typename T, int L ), T, L ),
    ( HBitVector< 1024 >, 11431 ),
    ( HBitVector< 2048 >, 11431 ),
    ( HBitVector< 4096 >, 11431 ) )
{
  using hbv_type = T;
  using bitset_type = typename hbv_type::bitset_type;
  using policy_type = typename hbv_type::policy_type;
  using member_type = typename hbv_type::member_type;
  using size_type = typename hbv_type::size_type;

  std::size_t len = L;
  auto width = gum::widthof< bitset_type >::value;

  std::size_t nrows = 5;
  std::size_t nnz = 34;
  GIVEN( "A simple matrix with " + std::to_string( nrows ) + " rows and "
         + std::to_string( nnz ) + " non-zero values in rCRS format" )
  {
    Kokkos::View< uint32_t* > e( "entries", nnz );
    Kokkos::View< uint32_t* > row_map( "row_map", nrows + 1 );

    auto h_e = Kokkos::create_mirror_view( e );
    auto h_row_map = Kokkos::create_mirror_view( row_map );

    h_row_map( 0 ) = 0;
    h_e( 0 ) = 0;
    h_e( 1 ) = 63;
    h_e( 2 ) =  500;
    h_e( 3 ) =  511;
    h_e( 4 ) =  512;
    h_e( 5 ) =  639;
    h_e( 6 ) = 1472;
    h_e( 7 ) = 1535;
    h_e( 8 ) = 4091;
    h_e( 9 ) = 4200;
    h_row_map( 1 ) = 10;
    h_e( 10 ) = 1;
    h_e( 11 ) = 64;
    h_e( 12 ) = 500;
    h_e( 13 ) = 639;
    h_e( 14 ) = 1471;
    h_e( 15 ) = 1555;
    h_e( 16 ) = 11300;
    h_e( 17 ) = 11430;
    h_row_map( 2 ) = 18;
    h_row_map( 3 ) = 18;
    h_e( 18 ) = 0;
    h_e( 19 ) = 11430;
    h_row_map( 4 ) = 20;
    h_e( 20 ) =   32;
    h_e( 21 ) =   32;
    h_e( 22 ) =   65;
    h_e( 23 ) =  130;
    h_e( 24 ) =  140;
    h_e( 25 ) =  514;
    h_e( 26 ) =  543;
    h_e( 27 ) = 1034;
    h_e( 28 ) = 1036;
    h_e( 29 ) = 1036;
    h_e( 30 ) = 1038;
    h_e( 31 ) = 2080;
    h_e( 32 ) = 8911;
    h_e( 33 ) = 8911;
    h_row_map( 5 ) = 34;

    size_type true_crs_nnz = 0;
    for ( auto i = 0u; i < h_e.extent( 0 ) - 1; i += 2 )
      true_crs_nnz += h_e( i + 1 ) - h_e( i ) + 1;

    Kokkos::deep_copy( e, h_e );
    Kokkos::deep_copy( row_map, h_row_map );

    WHEN( "Accumulate entries across using a hierarchical bitvector of length "
          + std::to_string( len ) + " and width " + std::to_string( width )
          + " (Team-Sequential Partitioning)" )
    {
      Kokkos::View< unsigned char* > flags( "flags", nnz / 2 );
      Kokkos::View< size_type > crs_nnz( "crs_nnz" );

      // Zero initialise `flags`
      Kokkos::parallel_for(
          "psi::test_hbitvector::initialise_flags", nnz / 2,
          KOKKOS_LAMBDA ( const uint64_t i ) {
            flags( i ) = 0;
            if ( i == 0 ) crs_nnz() = 0;
          } );

      auto policy = policy_type( nrows, Kokkos::AUTO );
      hbv_type::set_scratch_size( policy, len );
      Kokkos::parallel_for(
          "psi::test_hbitvector::set_range", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            hbv_type hbv( tm, len, ( row + 1 ) * 1000 );
            auto e_idx = row_map( row );
            auto e_end = row_map( row + 1 );

            hbv.clear_l1( tm );
            hbv.clear_l2( tm );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, e_idx / 2, e_end / 2 ),
                [&]( const uint64_t jj ) {
                  auto j = jj * 2;
                  auto s = e( j );
                  auto f = e( j + 1 );
                  hbv.set( tm, s, f, TeamSequentialPartition{} );
                } );

            tm.team_barrier();

            size_type row_nnz = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange( tm, e_idx / 2, e_end / 2 ),
                [=]( const uint64_t jj, size_type& lrow_nnz ) {
                  auto j = jj * 2;
                  auto s = e( j );
                  auto f = e( j + 1 );
                  flags( jj ) = 1;
                  size_type rng_nnz = 0;
                  for ( auto i = s; i <= f; ++i )
                    if ( !hbv[ i ] ) {
                      flags( jj ) = 0;
                      break;
                    }
                    else ++rng_nnz;
                  lrow_nnz += rng_nnz;
                }, row_nnz );

            Kokkos::single( Kokkos::PerTeam( tm ), [=](){
              Kokkos::atomic_add( &crs_nnz(), row_nnz );
            } );
          } );

      THEN( "All bits within the non-zero ranges should be set" )
      {
        std::size_t all_set = 0;
        Kokkos::parallel_reduce(
            "psi::test_hbit_vector::set_range_assess",
            flags.extent( 0 ),
            KOKKOS_LAMBDA ( const uint64_t i, std::size_t& all_set_local ) {
              if ( flags( i ) == 1 )
                all_set_local += 1;
              else
                printf( "i: %" PRIu64 "\n", i );
            },
            all_set );

        auto h_crs_nnz = Kokkos::create_mirror_view( crs_nnz );
        Kokkos::deep_copy( h_crs_nnz, crs_nnz );

        REQUIRE( all_set == flags.extent( 0 ) );
        REQUIRE( h_crs_nnz() == true_crs_nnz );
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG(
    "Bitwise operations 'cnt', 'msb', and 'lsb' on bitsets", "[hbitvector]",
    ( ( typename T, int L ), T, L ),
    ( HBitVector< 1024 >, 11431 ),
    ( HBitVector< 2048 >, 11431 ),
    ( HBitVector< 4096 >, 11431 ) )
{
  using hbv_type = T;
  using bitset_type = typename hbv_type::bitset_type;
  using policy_type = typename hbv_type::policy_type;
  using member_type = typename hbv_type::member_type;
  using size_type = typename hbv_type::size_type;

  std::size_t len = L;
  auto width = gum::widthof< bitset_type >::value;

  std::size_t nrows = 5;
  std::size_t nnz = 34;
  GIVEN( "A simple matrix with " + std::to_string( nrows ) + " rows and "
         + std::to_string( nnz ) + " non-zero values in rCRS format" )
  {
    Kokkos::View< uint32_t* > e( "entries", nnz );
    Kokkos::View< uint32_t* > row_map( "row_map", nrows + 1 );

    auto h_e = Kokkos::create_mirror_view( e );
    auto h_row_map = Kokkos::create_mirror_view( row_map );

    h_row_map( 0 ) = 0;
    h_e( 0 ) = 0;
    h_e( 1 ) = 63;
    h_e( 2 ) =  500;
    h_e( 3 ) =  511;
    h_e( 4 ) =  512;
    h_e( 5 ) =  639;
    h_e( 6 ) = 1472;
    h_e( 7 ) = 1535;
    h_e( 8 ) = 4091;
    h_e( 9 ) = 4200;
    h_row_map( 1 ) = 10;
    h_e( 10 ) = 1;
    h_e( 11 ) = 64;
    h_e( 12 ) = 500;
    h_e( 13 ) = 639;
    h_e( 14 ) = 1471;
    h_e( 15 ) = 1555;
    h_e( 16 ) = 11300;
    h_e( 17 ) = 11430;
    h_row_map( 2 ) = 18;
    h_row_map( 3 ) = 18;
    h_e( 18 ) = 0;
    h_e( 19 ) = 11430;
    h_row_map( 4 ) = 20;
    h_e( 20 ) =   32;
    h_e( 21 ) =   32;
    h_e( 22 ) =   65;
    h_e( 23 ) =  130;
    h_e( 24 ) =  140;
    h_e( 25 ) =  514;
    h_e( 26 ) =  543;
    h_e( 27 ) = 1034;
    h_e( 28 ) = 1036;
    h_e( 29 ) = 1036;
    h_e( 30 ) = 1038;
    h_e( 31 ) = 2080;
    h_e( 32 ) = 8911;
    h_e( 33 ) = 8911;
    h_row_map( 5 ) = 34;

    Kokkos::deep_copy( e, h_e );
    Kokkos::deep_copy( row_map, h_row_map );

    WHEN( "Apply bitwise operations on bitsets in hbitvector of width "
          + std::to_string( width )
          + " populated by nnz values (Team-Sequential Partitioning)" )
    {
      Kokkos::View< unsigned char* > msb_flags( "flags", hbv_type::num_bitsets( len ) );
      Kokkos::View< unsigned char* > lsb_flags( "flags", hbv_type::num_bitsets( len ) );
      Kokkos::View< size_type > crs_nnz( "crs_nnz" );

      // Zero initialise `msb_flags` and `lsb_flags`
      Kokkos::parallel_for(
          "psi::test_hbitvector::initialise_flags",
          hbv_type::num_bitsets( len ), KOKKOS_LAMBDA ( const uint64_t i ) {
            msb_flags( i ) = 0;
            lsb_flags( i ) = 0;
            if ( i == 0 ) crs_nnz() = 0;
          } );

      auto policy = policy_type( nrows, Kokkos::AUTO );
      hbv_type::set_scratch_size( policy, len );
      Kokkos::parallel_for(
          "psi::test_hbitvector::set_range", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            hbv_type hbv( tm, len, ( row + 1 ) * 1000 );
            auto e_idx = row_map( row );
            auto e_end = row_map( row + 1 );

            hbv.clear_l1( tm );
            hbv.clear_l2( tm );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, e_idx / 2, e_end / 2 ),
                [&]( const uint64_t jj ) {
                  auto j = jj * 2;
                  auto s = e( j );
                  auto f = e( j + 1 );
                  hbv.set( tm, s, f, TeamSequentialPartition{} );
                } );

            tm.team_barrier();

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t j ) {
                  bitset_type mask = hbv_type::BITSET_ALL_SET
                                     - ( hbv_type::BITSET_ALL_SET >> 1 );
                  auto x = hbv( j );
                  auto truth = ( x & mask ) ? 1u : 0u;
                  msb_flags( j ) = ( hbv_type::msb( x ) == truth );

                  truth = ( x & hbv_type::BITSET_ONE );
                  lsb_flags( j ) = ( hbv_type::lsb( x ) == truth );
                } );

            size_type thread_nnz = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t j, size_type& l_nnz ) {
                  l_nnz += hbv_type::cnt( hbv( j ) );
                }, thread_nnz );

            Kokkos::single( Kokkos::PerTeam( tm ), [=]() {
              Kokkos::atomic_add( &crs_nnz(), thread_nnz );
            } );
          } );

      THEN( "Counting all set bits in hbitvector of all row should give the "
            "nnz of input matrix in CRS format" )
      {
        size_type true_crs_nnz = 0;
        Kokkos::parallel_reduce(
            "psi::test_hbit_vector::compute_crs_nnz",
            policy_type( nrows, Kokkos::AUTO ),
            KOKKOS_LAMBDA ( const member_type& tm, size_type& tcnnz ) {
              auto row = tm.league_rank();
              auto e_idx = row_map( row );
              auto e_end = row_map( row + 1 );

              size_type team_nnz = 0;
              Kokkos::parallel_reduce(
                  Kokkos::TeamThreadRange( tm, e_idx / 2, e_end / 2 ),
                  [&]( const uint64_t jj, size_type& row_nnz ) {
                    auto j = jj * 2;
                    auto s = e( j );
                    auto f = e( j + 1 );
                    row_nnz += f + 1 - s;
                  }, team_nnz );

              Kokkos::single( Kokkos::PerTeam( tm ),
                              [&]() { tcnnz += team_nnz; } );
            },
            true_crs_nnz );

        auto h_crs_nnz = Kokkos::create_mirror_view( crs_nnz );
        Kokkos::deep_copy( h_crs_nnz, crs_nnz );

        REQUIRE( h_crs_nnz() == true_crs_nnz );
      }

      THEN( "Calling 'msb' on a bitset should give its most significant bit" )
      {
        std::size_t all_set = 0;
        Kokkos::parallel_reduce(
            "psi::test_hbit_vector::set_range_assess",
            msb_flags.extent( 0 ),
            KOKKOS_LAMBDA ( const uint64_t i, std::size_t& all_set_local ) {
              if ( msb_flags( i ) == 1 )
                all_set_local += 1;
              else
                printf( "i: %" PRIu64 "\n", i );
            },
            all_set );

        REQUIRE( all_set == msb_flags.extent( 0 ) );
      }

      THEN( "Calling 'lsb' on a bitset should give its most significant bit" )
      {
        std::size_t all_set = 0;
        Kokkos::parallel_reduce(
            "psi::test_hbit_vector::set_range_assess",
            lsb_flags.extent( 0 ),
            KOKKOS_LAMBDA ( const uint64_t i, std::size_t& all_set_local ) {
              if ( lsb_flags( i ) == 1 )
                all_set_local += 1;
              else
                printf( "i: %" PRIu64 "\n", i );
            },
            all_set );

        REQUIRE( all_set == lsb_flags.extent( 0 ) );
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG(
    "Bitwise operations 'cnt01', 'cnt10', 'map01', and 'map10' on bitsets", "[hbitvector]",
    ( ( typename T, int L ), T, L ),
    ( HBitVector< 1024 >, 11431 ),
    ( HBitVector< 2048 >, 11431 ),
    ( HBitVector< 4096 >, 11431 ) )
{
  using hbv_type = T;
  using bitset_type = typename hbv_type::bitset_type;
  using policy_type = typename hbv_type::policy_type;
  using member_type = typename hbv_type::member_type;
  using size_type = typename hbv_type::size_type;

  std::size_t len = L;
  auto width = gum::widthof< bitset_type >::value;

  std::size_t nrows = 5;
  std::size_t nnz = 34;
  GIVEN( "A simple matrix with " + std::to_string( nrows ) + " rows and "
         + std::to_string( nnz ) + " non-zero values in rCRS format" )
  {
    Kokkos::View< uint32_t* > e( "entries", nnz );
    Kokkos::View< uint32_t* > row_map( "row_map", nrows + 1 );

    auto h_e = Kokkos::create_mirror_view( e );
    auto h_row_map = Kokkos::create_mirror_view( row_map );

    h_row_map( 0 ) = 0;
    h_e( 0 ) = 0;
    h_e( 1 ) = 63;
    h_e( 2 ) =  500;
    h_e( 3 ) =  511;
    h_e( 4 ) =  512;
    h_e( 5 ) =  639;
    h_e( 6 ) = 1472;
    h_e( 7 ) = 1535;
    h_e( 8 ) = 4091;
    h_e( 9 ) = 4200;
    h_row_map( 1 ) = 10;
    h_e( 10 ) = 1;
    h_e( 11 ) = 64;
    h_e( 12 ) = 500;
    h_e( 13 ) = 639;
    h_e( 14 ) = 1471;
    h_e( 15 ) = 1555;
    h_e( 16 ) = 11300;
    h_e( 17 ) = 11430;
    h_row_map( 2 ) = 18;
    h_row_map( 3 ) = 18;
    h_e( 18 ) = 0;
    h_e( 19 ) = 11430;
    h_row_map( 4 ) = 20;
    h_e( 20 ) =   32;
    h_e( 21 ) =   32;
    h_e( 22 ) =   65;
    h_e( 23 ) =  130;
    h_e( 24 ) =  140;
    h_e( 25 ) =  514;
    h_e( 26 ) =  543;
    h_e( 27 ) = 1034;
    h_e( 28 ) = 1036;
    h_e( 29 ) = 1036;
    h_e( 30 ) = 1038;
    h_e( 31 ) = 2080;
    h_e( 32 ) = 8911;
    h_e( 33 ) = 8911;
    h_row_map( 5 ) = 34;

    Kokkos::deep_copy( e, h_e );
    Kokkos::deep_copy( row_map, h_row_map );

    WHEN( "Counting '01's and '10's on a hbitvector of width "
          + std::to_string( width )
          + " populated by nnz values (Team-Sequential Partitioning)" )
    {
      // The answer
      size_type t_nnz = 32;
      Kokkos::View< uint32_t*, Kokkos::DefaultHostExecutionSpace > h_t_e(
          "true entries", t_nnz );
      Kokkos::View< uint32_t*, Kokkos::DefaultHostExecutionSpace > h_t_row_map(
          "true row_map", nrows + 1 );

      h_t_row_map( 0 ) = 0;
      h_t_e( 0 ) = 0;
      h_t_e( 1 ) = 63;
      h_t_e( 2 ) = 500;
      h_t_e( 3 ) = 639;
      h_t_e( 4 ) = 1472;
      h_t_e( 5 ) = 1535;
      h_t_e( 6 ) = 4091;
      h_t_e( 7 ) = 4200;
      h_t_row_map( 1 ) = 8;
      h_t_e( 8 ) = 1;
      h_t_e( 9 ) = 64;
      h_t_e( 10 ) = 500;
      h_t_e( 11 ) = 639;
      h_t_e( 12 ) = 1471;
      h_t_e( 13 ) = 1555;
      h_t_e( 14 ) = 11300;
      h_t_e( 15 ) = 11430;
      h_t_row_map( 2 ) = 16;
      h_t_row_map( 3 ) = 16;
      h_t_e( 16 ) = 0;
      h_t_e( 17 ) = 11430;
      h_t_row_map( 4 ) = 18;
      h_t_e( 18 ) = 32;
      h_t_e( 19 ) = 32;
      h_t_e( 20 ) = 65;
      h_t_e( 21 ) = 130;
      h_t_e( 22 ) = 140;
      h_t_e( 23 ) = 514;
      h_t_e( 24 ) = 543;
      h_t_e( 25 ) = 1034;
      h_t_e( 26 ) = 1036;
      h_t_e( 27 ) = 1036;
      h_t_e( 28 ) = 1038;
      h_t_e( 29 ) = 2080;
      h_t_e( 30 ) = 8911;
      h_t_e( 31 ) = 8911;
      h_t_row_map( 5 ) = 32;

      // Output views
      Kokkos::View< uint32_t* > c_rowmap( "acc_rowmap", nrows + 1 );
      auto h_c_rowmap = Kokkos::create_mirror_view( c_rowmap );

      // Allocating space required for hbitvector
      auto policy = policy_type( nrows, Kokkos::AUTO );
      hbv_type::set_scratch_size( policy, len );

      // Computing `c_rowmap`
      Kokkos::parallel_for(
          "psi::test_hbitvector::count_row_nnz", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            hbv_type hbv( tm, len, ( row + 1 ) * 1000 );
            auto e_idx = row_map( row );
            auto e_end = row_map( row + 1 );

            hbv.clear_l1( tm );
            hbv.clear_l2( tm );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, e_idx / 2, e_end / 2 ),
                [&]( const uint64_t jj ) {
                  auto j = jj * 2;
                  auto s = e( j );
                  auto f = e( j + 1 );
                  hbv.set( tm, s, f, TeamSequentialPartition{} );
                } );

            tm.team_barrier();

            size_type row_nnz = 0;
            Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t j, size_type& l_rnnz ) {
                  auto c = ( j != 0 ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                  auto x = hbv( j );
                  l_rnnz += hbv_type::cnt01( x, c ) + hbv_type::cnt10( x, c );
                },
                row_nnz );

            Kokkos::single( Kokkos::PerTeam( tm ), [=]() {
              c_rowmap( row + 1 ) = row_nnz;
              if ( row == 0 ) c_rowmap( 0 ) = 0;
            } );
          } );

      Kokkos::parallel_scan(
          "psi::test_hbitvector::compute_rowmap", nrows,
          KOKKOS_LAMBDA ( const uint64_t i, size_type& update, const bool final ) {
            const size_type val_ip1 = c_rowmap( i + 1 );
            update += val_ip1;
            if ( final )
              c_rowmap( i + 1 ) = update;
          } );

      Kokkos::deep_copy( h_c_rowmap, c_rowmap );

      THEN( "Total number of '01's and '10's should be equal to nnz" )
      {
        REQUIRE( h_c_rowmap( nrows ) == t_nnz );

        for ( auto i = 0u; i <= nrows; ++i ) {
          REQUIRE( h_c_rowmap( i ) == h_t_row_map( i ) );
        }
      }

      AND_WHEN( "'sel'ecting all set bits in the result of calling "
                "('map01' | 'map10') on bitsets" )
      {
        Kokkos::View< uint32_t* > c_e( "acc_entries", h_c_rowmap( nrows ) );
        auto h_c_e = Kokkos::create_mirror( c_e );

        // Calculating `c_e`
        Kokkos::parallel_for(
            "psi::test_hbitvector::accumulate_entries", policy,
            KOKKOS_LAMBDA ( const member_type& tm ) {
              auto row = tm.league_rank();
              hbv_type hbv( tm, len, ( row + 1 ) * 1000 );
              auto e_idx = row_map( row );
              auto e_end = row_map( row + 1 );

              hbv.clear_l1( tm );
              hbv.clear_l2( tm );

              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange( tm, e_idx / 2, e_end / 2 ),
                  [&]( const uint64_t jj ) {
                    auto j = jj * 2;
                    auto s = e( j );
                    auto f = e( j + 1 );
                    hbv.set( tm, s, f, TeamSequentialPartition{} );
                  } );

              tm.team_barrier();

              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                  [=]( const uint64_t j ) {
                    auto c = ( j != 0 ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                    auto x = hbv( j );
                    auto bounds
                        = hbv_type::map01( x, c ) | hbv_type::map10( x, c );

                    if ( bounds != 0 ) {
                      size_type c_idx = 0;
                      Kokkos::parallel_reduce(
                          Kokkos::ThreadVectorRange( tm, j ),
                          [=]( const uint64_t k, size_type& lc_idx ) {
                            auto c = ( k != 0 ) ? hbv_type::msb( hbv( k - 1 ) )
                                                : 0;
                            auto x = hbv( k );
                            lc_idx += hbv_type::cnt01( x, c )
                                      + hbv_type::cnt10( x, c );
                          },
                          c_idx );
                      c_idx += c_rowmap( row );

                      Kokkos::parallel_for(
                          Kokkos::ThreadVectorRange( tm,
                                                     hbv_type::cnt( bounds ) ),
                          [=]( const uint64_t k ) {
                            auto lidx = c_idx + k;
                            c_e( lidx ) = hbv_type::start_index( j )
                                          + hbv_type::sel( bounds, k + 1 )
                                          - lidx % 2;
                          } );
                    }
                  } );
            } );

        Kokkos::deep_copy( h_c_e, c_e );

        THEN( "Using 'sel', 'map01', and 'map10' on bitsets can reconstruct "
              "the entries" )
        {
          for ( auto i = 0u; i < t_nnz; ++i ) {
            REQUIRE( h_c_e( i ) == h_t_e( i ) );
          }
        }
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG(
    "Clearing L1 and L2 regions", "[hbitvector]",
    ( ( typename T, int L ), T, L ),
    ( HBitVector< 1024 >, 11431 ),
    ( HBitVector< 2048 >, 11431 ),
    ( HBitVector< 4096 >, 11431 ) )
{
  using hbv_type = T;
  using bitset_type = typename hbv_type::bitset_type;
  using policy_type = typename hbv_type::policy_type;
  using member_type = typename hbv_type::member_type;
  using size_type = typename hbv_type::size_type;

  std::size_t len = L;
  auto width = gum::widthof< bitset_type >::value;

  GIVEN( "A fully-set hierarchical bitvector of width "
         + std::to_string( width ) )
  {
    size_type nrows = 12;

    // Allocating space required for hbitvector
    auto policy = policy_type( nrows, Kokkos::AUTO );
    hbv_type::set_scratch_size( policy, len );

    WHEN( "L1 region is cleared completely" )
    {
      Kokkos::View< unsigned int > flag ( "" );
      auto h_flag = Kokkos::create_mirror_view( flag );

      h_flag() = 0;
      Kokkos::deep_copy( flag, h_flag );

      Kokkos::parallel_for(
          "psi::test_hbitvector::clear_l1", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            auto hbv = hbv_type( tm, len, row * 1000 );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [&hbv]( const uint64_t i ) {
                  hbv( i ) = hbv_type::BITSET_ALL_SET;
                } );

            hbv.clear_l1( tm );
            auto begin = hbv.l1_begin_bindex();
            auto end = begin + hbv.l1_num_bitsets();
            Kokkos::parallel_for( Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                                  [=]( const uint64_t i ) {
                                    if ( ( i < begin && hbv( i ) == 0 )
                                         || ( end <= i && hbv( i ) == 0 )
                                         || ( begin <= i && i < end && hbv( i ) != 0 ) ) {
                                      Kokkos::atomic_add( &flag(), 1 );
                                    }
                                  } );
          } );

      Kokkos::deep_copy( h_flag, flag );

      THEN( "All bitsets in L1 should be zero" )
      {
        REQUIRE( h_flag() == 0 );
      }
    }

    WHEN( "L2 region is cleared completely" )
    {
      Kokkos::View< unsigned int > flag ( "" );
      auto h_flag = Kokkos::create_mirror_view( flag );

      h_flag() = 0;
      Kokkos::deep_copy( flag, h_flag );

      Kokkos::parallel_for(
          "psi::test_hbitvector::clear_l2_all", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            auto hbv = hbv_type( tm, len, row * 1000 );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [&hbv]( const uint64_t i ) {
                  hbv( i ) = hbv_type::BITSET_ALL_SET;
                } );

            hbv.clear_l2( tm );
            auto begin = hbv.l1_begin_bindex();
            auto end = begin + hbv.l1_num_bitsets();
            Kokkos::parallel_for( Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                                  [=]( const uint64_t i ) {
                                    if ( ( i < begin && hbv( i ) != 0 )
                                         || ( end <= i && hbv( i ) != 0 )
                                         || ( begin <= i && i < end && hbv( i ) == 0 ) ) {
                                      Kokkos::atomic_add( &flag(), 1 );
                                    }
                                  } );
          } );

      Kokkos::deep_copy( h_flag, flag );

      THEN( "All bitsets in L2 should be zero" )
      {
        REQUIRE( h_flag() == 0 );
      }
    }

    WHEN( "A region of L2 indicated by by local bitset indices is cleared" )
    {
      Kokkos::View< unsigned int > flag ( "" );
      auto h_flag = Kokkos::create_mirror_view( flag );

      h_flag() = 0;
      Kokkos::deep_copy( flag, h_flag );

      size_type clen = 70;  // number of bitsets to clear in L2
      Kokkos::parallel_for(
          "psi::test_hbitvector::clear_l2_lbidx", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            auto hbv = hbv_type( tm, len, row * 1000 );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [&hbv]( const uint64_t i ) {
                  hbv( i ) = hbv_type::BITSET_ALL_SET;
                } );

            auto lb_bidx = ( hbv.l2_num_bitsets() - clen ) / 2;
            hbv.clear_l2( tm, lb_bidx, lb_bidx + clen );
            auto begin = lb_bidx + hbv.l1_begin_bindex() + hbv.l1_num_bitsets();
            if ( begin >= hbv.num_bitsets() ) begin -= hbv.num_bitsets();
            auto end = begin + clen;
            if ( end >= hbv.num_bitsets() ) end -= hbv.num_bitsets();

            typename hbv_type::bitset_type value = 0;
            if ( end < begin ) {
              auto tmp = end;
              end = begin;
              begin = tmp;
              value = hbv_type::BITSET_ALL_SET;
            }

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t i ) {
                  if ( ( i < begin && hbv( i ) == value )
                       || ( end <= i && hbv( i ) == value )
                       || ( begin <= i && i < end && hbv( i ) != value ) ) {
                    Kokkos::atomic_add( &flag(), 1 );
                  }
            } );
          } );

      Kokkos::deep_copy( h_flag, flag );

      THEN( "All bitsets in that L2 region should be zero" )
      {
        REQUIRE( h_flag() == 0 );
      }
    }

    WHEN( "A region of L2 indicated by by global bitset indices is cleared" )
    {
      Kokkos::View< unsigned int > flag ( "" );
      auto h_flag = Kokkos::create_mirror_view( flag );

      h_flag() = 0;
      Kokkos::deep_copy( flag, h_flag );

      size_type clen = 25;
      Kokkos::parallel_for(
          "psi::test_hbitvector::clear_l2_bidx", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            auto hbv = hbv_type( tm, len, row * 1000 );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [&hbv]( const uint64_t i ) {
                  hbv( i ) = hbv_type::BITSET_ALL_SET;
                } );

            auto begin = hbv.l1_begin_bindex() + hbv.l1_num_bitsets();
            // space_len: the length of bigger space between lo-L2 or hi-L2
            size_type space_len = hbv.num_bitsets() - begin;
            if ( hbv.l1_begin_bindex() >= hbv.num_bitsets() / 2 ) {
              space_len = hbv.l1_begin_bindex();
              begin = 0;
            }

            assert( space_len > clen );

            begin += ( space_len - clen ) / 2;
            auto end = begin + clen;
            hbv.clear_l2_by_bidx( tm, begin, end );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t i ) {
                  if ( ( i < begin && hbv( i ) == 0 )
                       || ( end <= i && hbv( i ) == 0 )
                       || ( begin <= i && i < end && hbv( i ) != 0 ) ) {
                    Kokkos::atomic_add( &flag(), 1 );
                  }
                } );
          } );

      Kokkos::deep_copy( h_flag, flag );

      THEN( "All bitsets in the region should be zero" )
      {
        REQUIRE( h_flag() == 0 );
      }
    }

    WHEN( "A region of L2 indicated is cleared by global bit indices (non-zero end offset)" )
    {
      Kokkos::View< unsigned int > flag ( "flag" );
      auto h_flag = Kokkos::create_mirror_view( flag );

      h_flag() = 0;
      Kokkos::deep_copy( flag, h_flag );

      size_type clen = 25;
      size_type s_offset = psi::random::random_index( width );
      size_type e_offset = psi::random::random_integer( 1, width - 1 );
      Kokkos::parallel_for(
          "psi::test_hbitvector::clear_l2_idx", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            auto hbv = hbv_type( tm, len, row * 1000 );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [&hbv]( const uint64_t i ) {
                  hbv( i ) = hbv_type::BITSET_ALL_SET;
                } );

            auto begin = hbv.l1_begin_bindex() + hbv.l1_num_bitsets();
            // space_len: the length of bigger space between lo-L2 or hi-L2
            size_type space_len = hbv.num_bitsets() - begin;
            if ( hbv.l1_begin_bindex() >= hbv.num_bitsets() / 2 ) {
              space_len = hbv.l1_begin_bindex();
              begin = 0;
            }

            assert( space_len > clen );

            begin += ( space_len - clen ) / 2;
            auto end = begin + clen;
            hbv.clear_l2_by_idx( tm, hbv_type::start_index( begin ) + s_offset,
                                 hbv_type::start_index( end ) + e_offset );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t i ) {
                  if ( ( i < begin && hbv( i ) == 0 )
                       || ( end < i && hbv( i ) == 0 )
                       || ( begin <= i && i <= end && hbv( i ) != 0 ) ) {  // end bitset should be cleared
                    Kokkos::atomic_add( &flag(), 1 );
                  }
                } );
          } );

      Kokkos::deep_copy( h_flag, flag );

      THEN( "All bitsets in the region should be zero" )
      {
        REQUIRE( h_flag() == 0 );
      }
    }

    WHEN( "A region of L2 indicated is cleared by global bit indices (zero end offset)" )
    {
      Kokkos::View< unsigned int > flag ( "flag" );
      auto h_flag = Kokkos::create_mirror_view( flag );

      h_flag() = 0;
      Kokkos::deep_copy( flag, h_flag );

      size_type clen = 25;
      size_type s_offset = psi::random::random_index( width );
      Kokkos::parallel_for(
          "psi::test_hbitvector::clear_l2_idx", policy,
          KOKKOS_LAMBDA ( const member_type& tm ) {
            auto row = tm.league_rank();
            auto hbv = hbv_type( tm, len, row * 1000 );

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [&hbv]( const uint64_t i ) {
                  hbv( i ) = hbv_type::BITSET_ALL_SET;
                } );

            auto begin = hbv.l1_begin_bindex() + hbv.l1_num_bitsets();
            // space_len: the length of bigger space between lo-L2 or hi-L2
            size_type space_len = hbv.num_bitsets() - begin;
            if ( hbv.l1_begin_bindex() >= hbv.num_bitsets() / 2 ) {
              space_len = hbv.l1_begin_bindex();
              begin = 0;
            }

            assert( space_len > clen );

            begin += ( space_len - clen ) / 2;
            auto end = begin + clen;
            hbv.clear_l2_by_idx( tm, hbv_type::start_index( begin ) + s_offset,
                                 hbv_type::start_index( end ) /* no offset */);

            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, hbv.num_bitsets() ),
                [=]( const uint64_t i ) {
                  if ( ( i < begin && hbv( i ) == 0 )
                       || ( end <= i && hbv( i ) == 0 )
                       || ( begin <= i && i < end && hbv( i ) != 0 ) ) {  // end bitset is not cleared
                    Kokkos::atomic_add( &flag(), 1 );
                  }
                } );
          } );

      Kokkos::deep_copy( h_flag, flag );

      THEN( "All bitsets in the region should be zero" )
      {
        REQUIRE( h_flag() == 0 );
      }
    }
  }
}
