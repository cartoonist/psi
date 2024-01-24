/**
 *    @file  test_range_sparse.cpp
 *   @brief  Test scenarios for range sparse matrix operations module
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu Sep 07, 2023  21:56
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spgemm.hpp>
#include <KokkosSparse_spadd.hpp>
#include <psi/range_sparse.hpp>
#include <psi/graph.hpp>
#include <gum/graph.hpp>
#include <gum/gfa_utils.hpp>
#include <string>

#include "test_base.hpp"


using namespace psi;

SCENARIO( "Sanity check of Kokkos views created by CRSMatrix", "[range_sparse]" )
{
  typedef int scalar_t;
  typedef Kokkos::DefaultHostExecutionSpace host_space;
  typedef KokkosSparse::CrsMatrix< scalar_t, int32_t, host_space > xcrsmatrix_t;
  typedef psi::CRSMatrix< psi::crs_matrix::RangeDynamic, bool, uint32_t, uint64_t > range_crsmatrix_t;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  GIVEN( "A sequence graph" )
  {
    std::string graph_path = test_data_dir + "/small/x.gfa";
    graph_type graph;
    gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );

    auto a = psi::util::adjacency_matrix< xcrsmatrix_t >( graph );
    range_crsmatrix_t ra( a );

    WHEN( "Create a Kokkos view on host for CRSMatrix" )
    {
      auto h_a_entries = ra.entries_view();
      auto h_a_row_map = ra.rowmap_view();

      auto a_entries = ra.entries_device_view( host_space{} );
      auto a_row_map = ra.rowmap_device_view( host_space{} );

      THEN( "The data should not be copied" )
      {
        REQUIRE( h_a_entries.data() == a_entries.data() );
        REQUIRE( h_a_row_map.data() == a_row_map.data() );
      }
    }

    WHEN( "Create a Device mirror view for CRSMatrix" )
    {
      auto h_a_entries = ra.entries_view();
      auto h_a_row_map = ra.rowmap_view();

      auto a_entries = ra.entries_device_view( Kokkos::DefaultExecutionSpace{} );
      auto a_row_map = ra.rowmap_device_view( Kokkos::DefaultExecutionSpace{} );

      auto ch_a_entries = Kokkos::create_mirror_view( a_entries );
      auto ch_a_row_map = Kokkos::create_mirror_view( a_row_map );

      Kokkos::deep_copy( ch_a_entries, a_entries );
      Kokkos::deep_copy( ch_a_row_map, a_row_map );

      THEN( "The data should be identical on both host and device" )
      {
        for ( std::size_t i = 0; i < a_entries.extent( 0 ); ++i )
          REQUIRE( h_a_entries( i ) == ch_a_entries( i ) );
        for ( std::size_t i = 0; i < a_row_map.extent( 0 ); ++i )
          REQUIRE( h_a_row_map( i ) == ch_a_row_map( i ) );
      }
    }
  }
}

template< typename TXCRSMatrix1, typename TXCRSMatrix2 >
TXCRSMatrix1
copy_xcrs( TXCRSMatrix2 mat )
{
  using dst_matrix_type  = TXCRSMatrix1;
  using src_matrix_type  = TXCRSMatrix2;

  using ordinal_t = typename src_matrix_type::ordinal_type;
  using size_type = typename src_matrix_type::size_type;
  using scalar_t = typename src_matrix_type::value_type;

  static_assert(
      std::is_same< ordinal_t, typename dst_matrix_type::ordinal_type >::value,
      "ordinal types of two matrices should be identical" );
  static_assert(
      std::is_same< size_type, typename dst_matrix_type::size_type >::value,
      "size types of two matrices should be identical" );
  static_assert(
      std::is_same< scalar_t, typename dst_matrix_type::value_type >::value,
      "value types of two matrices should be identical" );

  using dst_exec_space = typename dst_matrix_type::execution_space;
  using dst_graph_type   = typename dst_matrix_type::staticcrsgraph_type;

  auto row_map =
      Kokkos::create_mirror_view_and_copy( dst_exec_space{},
                                           mat.graph.row_map );
  auto entries =
      Kokkos::create_mirror_view_and_copy( dst_exec_space{},
                                           mat.graph.entries );
  auto values
      = Kokkos::create_mirror_view_and_copy( dst_exec_space{},
                                             mat.values );

  // using dst_row_map_type = typename dst_graph_type::row_map_type;
  // using dst_entries_type = typename dst_graph_type::entries_type;
  // using dst_values_type = typename dst_matrix_type::values_type;

  // typename dst_row_map_type::non_const_type row_map( "rowmap", mat.numRows() + 1 );
  // typename dst_entries_type::non_const_type entries( "entries", mat.nnz() );
  // typename dst_values_type::non_const_type values( "values", mat.nnz() );

  // Kokkos::deep_copy( row_map, mat.graph.row_map );
  // Kokkos::deep_copy( entries, mat.graph.entries );
  // Kokkos::deep_copy( values, mat.values );

  dst_graph_type crs_graph( entries, row_map );
  return dst_matrix_type( "moved", mat.numRows(), values, crs_graph );
}

template< typename TXCRSMatrix, typename TRCRSMatrix >
inline bool
is_same_host( TXCRSMatrix& x_mat, TRCRSMatrix& r_mat )
{
  typedef Kokkos::DefaultHostExecutionSpace host_space;
  typedef typename TXCRSMatrix::size_type size_type;
  typedef typename TXCRSMatrix::ordinal_type ordinal_type;
  typedef make_basic_t< TRCRSMatrix > TBCRSMatrix;

  TBCRSMatrix b_mat;
  b_mat.assign( r_mat );

  if ( x_mat.numRows() != b_mat.numRows() || x_mat.numCols() != b_mat.numCols()
       || x_mat.nnz() != b_mat.nnz() )
    return false;

  {
    ordinal_type num_matches = 0;
    Kokkos::parallel_reduce(
        "psi::test_range_sparse::compare_rowmap",
        Kokkos::RangePolicy< host_space >( 0, x_mat.numRows() + 1 ),
        KOKKOS_LAMBDA( const uint64_t i, ordinal_type& l_nm ) {
          if ( x_mat.graph.row_map( i ) == b_mat.rowMap( i ) ) ++l_nm;
        },
        num_matches );

    if ( num_matches != x_mat.numRows() + 1 ) return false;
  }

  {
    size_type num_matches = 0;
    Kokkos::parallel_reduce(
        "psi::test_range_sparse::compare_entries",
        Kokkos::RangePolicy< host_space >( 0, x_mat.nnz() ),
        KOKKOS_LAMBDA( const uint64_t i, size_type& l_nm ) {
          if ( x_mat.graph.entries( i ) == b_mat.entry( i ) ) ++l_nm;
        },
        num_matches );

    if ( num_matches != x_mat.nnz() ) return false;
  }

  return true;
}

template< typename TXCRSMatrix, typename TRCRSMatrix >
inline bool
is_same( TXCRSMatrix& x_mat, TRCRSMatrix& r_mat )
{
  typedef Kokkos::DefaultHostExecutionSpace host_space;
  typedef typename TXCRSMatrix::memory_space xcrs_memory_space;
  typedef typename TXCRSMatrix::HostMirror xcrs_host_mirror;

  if constexpr ( Kokkos::SpaceAccessibility<
                     host_space, xcrs_memory_space >::accessible ) {
    return is_same_host( x_mat, r_mat );
  }
  else {
    auto h_xmat = copy_xcrs< xcrs_host_mirror >( x_mat );
    return is_same_host( h_xmat, r_mat );
  }
}

template< typename TXCRSMatrix >
inline TXCRSMatrix
kokkos_kernels_spgemm( TXCRSMatrix const& a, TXCRSMatrix const& b )
{
  typedef typename TXCRSMatrix::ordinal_type ordinal_t;
  typedef typename TXCRSMatrix::size_type size_type;
  typedef typename TXCRSMatrix::value_type scalar_t;
  typedef typename TXCRSMatrix::execution_space execution_space;
  typedef typename TXCRSMatrix::memory_space memory_space;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<
    size_type, ordinal_t, scalar_t,
    execution_space, memory_space, memory_space > kernel_handle_t;

  kernel_handle_t handle;
  handle.set_team_work_size( 16 );
  handle.set_dynamic_scheduling( true );

  KokkosSparse::SPGEMMAlgorithm spgemm_algorithm =
    // spgemm algorithm, limited by configuration at compile-time and set via the handle
    // Other options: {SPGEMM_KK_SPEED, SPGEMM_KK_MEMSPEED, SPGEMM_MKL}
      KokkosSparse::SPGEMM_KK_MEMORY;
  handle.create_spgemm_handle( spgemm_algorithm );

  TXCRSMatrix c;

  {
#ifdef PSI_STATS
    Kokkos::Timer timer;
#endif

    KokkosSparse::spgemm_symbolic( handle, a, false, b, false, c );
    execution_space{}.fence();

#ifdef PSI_STATS
    auto duration = timer.seconds();
    std::cout << "Kokkos::SpGEMM_symbolic time: " << duration * 1000 << "ms"
              << std::endl;
#endif
  }

  {
#ifdef PSI_STATS
    Kokkos::Timer timer;
#endif

    KokkosSparse::spgemm_numeric( handle, a, false, b, false, c );
    execution_space{}.fence();

#ifdef PSI_STATS
    auto duration = timer.seconds();
    std::cout << "Kokkos::SpGEMM_numeric time: " << duration * 1000 << "ms"
              << std::endl;
#endif
  }

//  Kokkos::parallel_for(
//      "psi::test_range_sparse::set_values",
//      Kokkos::RangePolicy< execution_space >( 0, c.nnz() ),
//      KOKKOS_LAMBDA ( const uint64_t i ) {
//        c.values( i ) = 1;
//      } );

  return c;
}

template< typename TXCRSMatrix >
inline TXCRSMatrix
kokkos_kernels_power( TXCRSMatrix const& a, unsigned int n )
{
  assert( a.numRows() == a.numCols() );

#ifdef PSI_STATS
  Kokkos::Timer timer;
#endif
  auto c = create_identity_matrix< TXCRSMatrix >( a.numRows() );
  auto a_copy = a;

  while ( true ) {
    if ( n & 1 ) c = kokkos_kernels_spgemm( c, a_copy );
    n = n >> 1;
    if ( n <= 0 ) break;
    a_copy = kokkos_kernels_spgemm( a_copy, a_copy );
  }

  typename TXCRSMatrix::execution_space{}.fence();

#ifdef PSI_STATS
  auto duration = timer.seconds();
  std::cout << "KokkosKernels::power time: " << duration * 1000 << "ms"
            << std::endl;
#endif

  return c;
}

TEMPLATE_SCENARIO_SIG( "Validation and verification of range SpGEMM", "[range_sparse]",
    ( ( typename TSpec, typename TScalar, typename TOrdinal, typename TSize, int N, int NNZ ),
      TSpec, TScalar, TOrdinal, TSize, N, NNZ ),
    ( ( crs_matrix::RangeDynamic ), char, int32_t, uint64_t, 6521, 200000 ) )
{
  typedef Kokkos::DefaultExecutionSpace device_t;
  typedef CRSMatrix< TSpec, bool, TOrdinal, TSize > rcrsmatrix_t;
  typedef KokkosSparse::CrsMatrix< TScalar, TOrdinal, device_t > xcrsmatrix_t;
  typedef typename xcrsmatrix_t::HostMirror xcrs_host_mirror;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  GIVEN( "A random square matrix of order " + std::to_string( N ) + " with "
         + std::to_string( NNZ ) + " non-zero values" )
  {
    rcrsmatrix_t rrand_mat;
    xcrsmatrix_t xrand_mat = create_random_binary_matrix< xcrsmatrix_t >( N, NNZ, rrand_mat );

    REQUIRE( rrand_mat.nnz() == xrand_mat.nnz() );
    REQUIRE( rrand_mat.nnz() == NNZ );
    REQUIRE( is_same( xrand_mat, rrand_mat ) );

    WHEN( "It is multiplied to itself" )
    {
      auto xc = kokkos_kernels_spgemm( xrand_mat, xrand_mat );
      auto rc = range_spgemm( rrand_mat, rrand_mat );
      REQUIRE( is_same( xc, rc ) );
    }
  }

  GIVEN( "A sequence graph" )
  {
    std::string graph_path = test_data_dir + "/middle/m.gfa";
    graph_type graph;
    gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );
    auto h_a = psi::util::adjacency_matrix< xcrs_host_mirror >( graph );
    auto a = copy_xcrs< xcrsmatrix_t >( h_a );
    rcrsmatrix_t ra( h_a );

    REQUIRE( ra.nnz() == a.nnz() );
    REQUIRE( is_same( h_a, ra ) );

    WHEN( "It is multiplied to itself" )
    {
      auto xc = kokkos_kernels_spgemm( a, a );
      auto rc = range_spgemm( ra, ra );
      REQUIRE( is_same( xc, rc ) );
    }
  }
}

TEMPLATE_SCENARIO_SIG( "Validation and verification of range power", "[range_sparse]",
    ( ( typename TSpec, typename TScalar, typename TOrdinal, typename TSize, int N, int NNZ, int K ),
      TSpec, TScalar, TOrdinal, TSize, N, NNZ, K ),
    ( ( crs_matrix::RangeDynamic ), char, int32_t, uint64_t, 6521, 4000, 100 ) )
{
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef execution_space::device_type device_t;
  typedef CRSMatrix< TSpec, bool, TOrdinal, TSize > rcrsmatrix_t;
  typedef KokkosSparse::CrsMatrix< TScalar, TOrdinal, device_t > xcrsmatrix_t;
  typedef typename xcrsmatrix_t::HostMirror xcrs_host_mirror;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  GIVEN( "A random square matrix of order " + std::to_string( N ) + " with "
         + std::to_string( NNZ ) + " non-zero values" )
  {
    rcrsmatrix_t rrand_mat;
    xcrsmatrix_t xrand_mat = create_random_binary_matrix< xcrsmatrix_t >( N, NNZ, rrand_mat );

    REQUIRE( rrand_mat.nnz() == xrand_mat.nnz() );
    REQUIRE( rrand_mat.nnz() == NNZ );
    REQUIRE( is_same( xrand_mat, rrand_mat ) );

    WHEN( "It is raised to the power of " + std::to_string( K ) )
    {
      auto xc = kokkos_kernels_power( xrand_mat, K );
      auto rc = range_power( rrand_mat, K );
      REQUIRE( is_same( xc, rc ) );
    }
  }

  GIVEN( "A sequence graph" )
  {
    std::string graph_path = test_data_dir + "/middle/m.gfa";
    graph_type graph;
    gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );
    auto h_a = psi::util::adjacency_matrix< xcrs_host_mirror >( graph );
    auto a = copy_xcrs< xcrsmatrix_t >( h_a );
    rcrsmatrix_t ra( h_a );

    REQUIRE( ra.nnz() == a.nnz() );
    REQUIRE( is_same( h_a, ra ) );

    WHEN( "It is raised to the power of " + std::to_string( K ) )
    {
      auto xc = kokkos_kernels_power( a, K );
      auto rc = range_power( ra, K );
      REQUIRE( is_same( xc, rc ) );
    }
  }
}
