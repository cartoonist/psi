/**
 *    @file  rcrs_benchmark.cpp
 *   @brief  Benchmark Range CRS matrix operations
 *
 *  Performance benchmarking of Range CRS matrix operations compared to general
 *  matrix representation as well as generic compressed row storage format (CRS).
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu May 11, 2023  18:34
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spgemm.hpp>
#include <KokkosSparse_spadd.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_NestedSort.hpp>

#include <psi/basic_types.hpp>
#include <psi/utils.hpp>
#include <psi/hbitvector.hpp>
#include <psi/graph.hpp>
#include <psi/range_sparse.hpp>
#include <gum/graph.hpp>
#include <gum/io_utils.hpp>

#include "vg/vg.pb.h"
#include "vg/stream.hpp"

using namespace psi;


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

#ifdef PSI_STATS
    Kokkos::Timer timer;
#endif

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

#ifdef PSI_STATS
  auto duration = timer.seconds();
  std::cout << "copy time: " << duration * 1000 << "ms"
            << std::endl;
#endif

  return dst_matrix_type( "moved", mat.numRows(), values, crs_graph );
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

  handle.destroy_spgemm_handle();

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
kokkos_kernels_spadd( TXCRSMatrix const& a, TXCRSMatrix const& b )
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
  handle.create_spadd_handle( true /* sorted rows */ );

  TXCRSMatrix c;

  {
#ifdef PSI_STATS
    Kokkos::Timer timer;
#endif

    KokkosSparse::spadd_symbolic( &handle, a, b, c );
    execution_space{}.fence();

#ifdef PSI_STATS
    auto duration = timer.seconds();
    std::cout << "Kokkos::SpAdd_symbolic time: " << duration * 1000 << "ms"
              << std::endl;
#endif
  }

  {
#ifdef PSI_STATS
    Kokkos::Timer timer;
#endif

    KokkosSparse::spadd_numeric( &handle, 1, a, 1, b, c );
    execution_space{}.fence();

#ifdef PSI_STATS
    auto duration = timer.seconds();
    std::cout << "Kokkos::SpAdd_numeric time: " << duration * 1000 << "ms"
              << std::endl;
#endif
  }

  handle.destroy_spadd_handle();

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
    if ( n == 0 ) break;
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

template< typename TOrdinal=int32_t, typename TScalar=char >
void
benchmark_range_spgemm_graph( const std::string& graph_path, int d, bool verbose )
{
  typedef TScalar scalar_t;
  typedef TOrdinal ordinal_t;
  typedef Kokkos::DefaultExecutionSpace execution_space;
  typedef execution_space::device_type device_t;
  typedef KokkosSparse::CrsMatrix< scalar_t, ordinal_t, device_t > xcrsmatrix_t;
  typedef typename xcrsmatrix_t::HostMirror xcrs_host_mirror;
  typedef typename xcrsmatrix_t::size_type size_type;
  typedef psi::CRSMatrix< psi::crs_matrix::RangeDynamic, bool, ordinal_t, size_type > range_crsmatrix_t;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  // Loading input graph
  auto parse_vg = []( std::istream& in ) -> vg::Graph {
    vg::Graph merged;
    std::function< void( vg::Graph& ) > handle_chunks =
      [&]( vg::Graph& other ) {
        gum::util::merge_vg( merged, static_cast< vg::Graph const& >( other ) );
      };
    stream::for_each( in, handle_chunks );
    return merged;
  };
  graph_type graph;
  gum::ExternalLoader< vg::Graph > loader{ parse_vg };

  std::cout << "Loading input graph..." << std::endl;
  gum::util::load( graph, graph_path, loader, true );

  std::cout << "Creating adjacency matrix..." << std::endl;
  std::vector< graph_type::rank_type > comp_ranks;
  gum::util::for_each_start_node( graph, [&comp_ranks]( auto rank, auto ) {
    comp_ranks.push_back( rank );
    return true;
  } );
  auto h_a = psi::util::adjacency_matrix< xcrs_host_mirror >(
      graph, comp_ranks[ 0 ], comp_ranks[ 1 ] );

  std::cout << "Copying adjacency matrix to device..." << std::endl;
  auto a = copy_xcrs< xcrsmatrix_t >( h_a );

  auto I = create_identity_matrix< xcrsmatrix_t >( a.numRows() );
  auto avi = kokkos_kernels_spadd( a, I );

  range_crsmatrix_t ra( h_a );
  auto rI = create_range_identity_matrix< range_crsmatrix_t >( a.numRows() );
  auto ravi = range_spadd( ra, rI );

  auto c = kokkos_kernels_power( avi, d );
  execution_space{}.fence();

  //if ( verbose ) {
  //  psi::print( a );
  //  psi::print( c );
  //}

  auto rc = range_power( ravi, d );
  execution_space{}.fence();

  double comp_rate = static_cast< double >( rc.nnz() ) / rc.rowMap( rc.numRows() );
  std::cout << "distance matrix of rank " << rc.numRows() << "x"
            << rc.numCols() << " holds " << rc.nnz()
            << " non-zero elements with compression rates " << comp_rate
            << " (" << rc.rowMap( rc.numRows() ) << ")" << std::endl;

  //if ( verbose ) {
  //  psi::print( rc, std::string( "RC" ) );
  //}
}

template< typename TOrdinal=uint32_t, typename TScalar=char >
void
benchmark_range_spgemm_random( TOrdinal n, std::size_t nnz, bool verbose )
{
  //auto a = create_random_binary_matrix< xcrsmatrix_t >( n, nnz, ra );
  //execution_space.fence();
  //range_crsmatrix_t rc;
  //range_spgemm( ra, ra );
}

template< typename TOrdinal=int32_t, typename TSize=std::size_t >
struct Options {
  typedef TOrdinal ordinal_type;
  typedef TSize size_type;

  std::string graph_path;
  TOrdinal n;
  TSize nnz;
  int d;
  bool verbose;

  Options() : graph_path( "" ), n( 0 ), nnz( 0 ), verbose( false ) {}
};

template< typename TOrdinal=int32_t, typename TSize=std::size_t >
Options< TOrdinal, TSize >
parse_arguments( int argc, char* argv[] )
{
  typedef Options< TOrdinal, TSize > option_type;

  option_type opts;

  // Parse command line arguments
  for ( int i = 0; i < argc; i++ ) {
    if ( ( strcmp( argv[ i ], "-g" ) == 0 ) || ( strcmp( argv[ i ], "--graph" ) == 0 ) ) {
      opts.graph_path = argv[ ++i ];
      std::cout << "graph <- " << opts.graph_path << std::endl;
    }
    if ( ( strcmp( argv[ i ], "-n" ) == 0 ) || ( strcmp( argv[ i ], "--order" ) == 0 ) ) {
      opts.n = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "n <- " << opts.n << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-z" ) == 0 ) || ( strcmp( argv[ i ], "--nnz" ) == 0 ) ) {
      opts.nnz = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "nnz <- " << opts.nnz << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-d" ) == 0 ) || ( strcmp( argv[ i ], "--dist" ) == 0 ) ) {
      opts.d = atoi( argv[ ++i ] );
      std::cout << "dist <- " << opts.d << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-v" ) == 0 ) || ( strcmp( argv[ i ], "--verbose" ) == 0 ) ) {
      opts.verbose = true;
      std::cout << "verbose <- true" << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-h" ) == 0 ) || ( strcmp( argv[ i ], "--help" ) == 0 ) ) {
      std::cout << "  Options:\n"
                << "  * Adjacency matrix as input\n"
                << "    --graph (-g) <path>: file path, determines input graph file path (gfa, vg)\n"
                << "  * Random matrix as input\n"
                << "    --order (-n) <int>:  exponent num, determines number of rows 2^num (default: 2^12 = 4096)\n"
                << "    --nnz (-z) <int>:    exponent num, determines total matrix size 2^num (default: 2^22 = 4096*1024)\n"
                << "    --dist (-d) <int>:   outer distance, i.e. fragment size (default: 100)\n"
                << "  --verbose (-v):        enable verbose output"
                << "  --help (-h):           print this message"
                << std::endl;
      exit( 1 );
    }
  }

  return opts;
}

template< typename TOrdinal=int32_t, typename TSize=std::size_t >
void
check_options( Options< TOrdinal, TSize >& opts )
{
  if ( opts.nnz == 0 && opts.n == 0 ) opts.nnz = pow( 2, 22 );

  // If opts.n is provided and opts.nnz is not, set it 1000*n.
  if ( opts.nnz == 0 ) opts.nnz = opts.n * 1000;

  // If opts.n is not provided, set order to 2^(ceil(log2(opts.nnz))/2).
  if ( opts.n == 0 ) {
    int power = 1;
    auto temp = opts.nnz;
    while ( ( temp /= 2 ) ) ++power;
    opts.n = pow( 2, power/2 );
  }

  std::cout << "nnz = " << opts.nnz << ", n = " << opts.n << std::endl;

  // Check sizes
  if ( opts.n < 0 || opts.nnz < 0 ) {
    std::cout << "nnz must be greater than zero" << std::endl;
    exit( EXIT_FAILURE );
  }

  if ( opts.nnz / opts.n > static_cast< TSize >( opts.n ) ) {
    std::cout << "Non-zero values cannot be fit (nnz > n*n)" << std::endl;
    exit( EXIT_FAILURE );
  }

  if ( opts.d == 0 ) opts.d = 100;
}

int
main( int argc, char* argv[] )
{
  auto opts = parse_arguments( argc, argv );
  check_options( opts );

  Kokkos::initialize( argc, argv );

  if ( opts.graph_path.empty() )
    benchmark_range_spgemm_random( opts.n, opts.nnz, opts.verbose );
  else
    benchmark_range_spgemm_graph( opts.graph_path, opts.d, opts.verbose );

  Kokkos::finalize();
  return EXIT_SUCCESS;
}
