/**
 *    @file  range_sparse.hpp
 *   @brief  Range sparse matrix operations
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu Sep 07, 2023  20:20
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_RANGE_SPARSE_HPP_
#define PSI_RANGE_SPARSE_HPP_

#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <KokkosKernels_Utils.hpp>

#include "crs_matrix.hpp"
#include "range_sparse_base.hpp"
#include "basic_types.hpp"
#include "utils.hpp"


namespace psi {
  /**
   *  @brief  Print a matrix of type `KokkosSparse::CrsMatrix`.
   *
   *  NOTE: All Views associated with the input matrix are assumed to be on
   *  device memory space. If they are host accessible deep copy does not copy
   *  anything.
   */
  template< typename TXCRSMatrix >
  static void
  print( const TXCRSMatrix& m, bool verbose=true, bool print_all=false )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename xcrsmatrix_t::index_type::non_const_type::HostMirror const_host_entries_t;
    typedef typename xcrsmatrix_t::values_type::non_const_type::HostMirror const_host_values_t;
    typedef typename xcrsmatrix_t::row_map_type::non_const_type::HostMirror const_host_row_map_t;

    auto label = m.values.label();
    auto nrows = m.numRows();
    auto ncols = m.numCols();
    auto nnz = m.nnz();

    std::cout << "[INFO] Matrix '" << label << "'"
              << " (" << nrows << "x" << ncols << ") with " << nnz << " non-zero elements:\n";

    if ( verbose ) {
      const_host_entries_t entries = Kokkos::create_mirror_view( m.graph.entries );
      const_host_values_t values = Kokkos::create_mirror_view( m.values );
      const_host_row_map_t rowmap = Kokkos::create_mirror_view( m.graph.row_map );

      Kokkos::deep_copy( entries, m.graph.entries );
      Kokkos::deep_copy( values, m.values );
      Kokkos::deep_copy( rowmap, m.graph.row_map );

      std::cout << "   ... ┬─\n";
      std::cout << "   ... ├─ entries  (" << entries.extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( entries, print_all );
      std::cout << "   ... ├─ values   (" << values.extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( values, print_all );
      std::cout << "   ... ╰─ row map  (" << rowmap.extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( rowmap, print_all );
      std::cout << "   ... " << std::endl;

      auto width = std::log( nrows ) + 1;
      std::cout << "   ... " << std::setw( width - 2 ) << label << " = [" << std::endl;
      for ( std::size_t i = 0; i < nrows; ++i ) {
        std::cout << "   ... " << std::setw( width ) << i << ": ";
        auto end = rowmap( i + 1 );
        for ( auto j = rowmap( i ); j < end; ++j ) {
          std::cout << " " << entries( j );
        }
        std::cout << "\n";
      }
      std::cout << "   ... " << std::setw( width + 2 ) << "]" << std::endl;
    }
  }

  /**
   *  @brief  Print a matrix of type 'Range' spacialised `psi::CRSMatrix`.
   *
   *  NOTE: All Views associated with the matrix are assumed to be on host
   * memory space.
   */
  template< typename TSpec, typename TOrdinal, typename TSize >
  static /* void */ std::enable_if_t<
      std::is_same< typename psi::crs_matrix::Group< TSpec >::type,
                    psi::crs_matrix::RangeGroup >::value >
  print( psi::CRSMatrix< TSpec, bool, TOrdinal, TSize >& m,
         std::string label = {}, bool verbose = true, bool print_all = false )
  {
    if ( label.empty() ) label = "A";

    std::cout << "[INFO] Matrix '" << label << "'"
              << " (" << m.numRows() << "x" << m.numCols() << ") with "
              << m.nnz() << " non-zero elements:\n";

    if ( verbose ) {
      std::cout << "   ... ┬─\n";
      std::cout << "   ... ├─ entries  (" << m.entries_view().extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( m.entries_view(), print_all );
      std::cout << "   ... ├─ values   (" << m.entries_view().extent( 0 ) << "): 1 ... 1\n";
      std::cout << "   ... ╰─ row map  (" << m.rowmap_view().extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( m.rowmap_view(), print_all );
      std::cout << "   ... " << std::endl;

      auto width = std::log( m.numRows() ) + 1;
      std::cout << "   ... " << std::setw( width - 2 ) << label << " = [" << std::endl;
      for ( std::size_t i = 0; i < m.numRows(); ++i ) {
        std::cout << "   ... " << std::setw( width ) << i << ": ";
        auto end = m.rowMap( i + 1 );
        for ( auto j = m.rowMap( i ); j < end; j += 2 ) {
          std::cout << " (" << m.entry( j ) << ", " << m.entry( j + 1 ) << ")";
        }
        std::cout << "\n";
      }
      std::cout << "   ... " << std::setw( width + 2 ) << "]" << std::endl;
    }
  }

  /**
   *  @brief  Create the identity matrix of order `n` in CRS format.
   *
   *  NOTE: `TXCRSMatrix` should be a `KokkosSparse::CrsMatrix`-like type.
   */
  template< typename TXCRSMatrix >
  inline TXCRSMatrix
  create_identity_matrix( typename TXCRSMatrix::ordinal_type n )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename xcrsmatrix_t::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > range_policy_t;

    typename xcrsmatrix_t::values_type::non_const_type i_values(
        Kokkos::ViewAllocateWithoutInitializing( "I" ), n );
    typename xcrsmatrix_t::row_map_type::non_const_type i_row_map(
        Kokkos::ViewAllocateWithoutInitializing( "I_rowmap" ), n + 1 );
    typename xcrsmatrix_t::index_type::non_const_type i_entries(
        Kokkos::ViewAllocateWithoutInitializing( "I_entries" ), n );

    Kokkos::parallel_for(
        "psi::crs_matrix::create_identity_matrix", range_policy_t( 0, n ),
        KOKKOS_LAMBDA ( const uint64_t i ) {
          i_values( i ) = 1;
          i_row_map( i + 1 ) = i + 1;
          i_entries( i ) = i;
          if ( i == 0 ) i_row_map( 0 ) = 0;
        } );

    return xcrsmatrix_t( "Identity Matrix", n, n, n, i_values, i_row_map, i_entries );
  }

  template< typename TRowMapDeviceView, typename TEntriesDeviceView >
  inline void
  create_range_identity_matrix(
      TRowMapDeviceView& i_rowmap, TEntriesDeviceView& i_entries,
      typename TEntriesDeviceView::value_type n )
  {
    typedef typename TEntriesDeviceView::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > range_policy_t;

    i_entries = TEntriesDeviceView(
        Kokkos::ViewAllocateWithoutInitializing( "I" ), n * 2 );
    i_rowmap = TRowMapDeviceView(
        Kokkos::ViewAllocateWithoutInitializing( "I_rowmap" ), n + 1 );

    Kokkos::parallel_for(
        "psi::crs_matrix::create_range_identity_matrix",
        range_policy_t( 0, n ), KOKKOS_LAMBDA ( const uint64_t ii ) {
          i_entries( ii * 2 ) = ii;
          i_entries( ii * 2 + 1 ) = ii;
          i_rowmap( ii + 1 ) = ( ii + 1 ) * 2;
          if ( ii == 0 ) i_rowmap( 0 ) = 0;
        } );
  }

  /**
   *  @brief  Create the identity matrix of order `n` in RCRS format.
   *
   *  NOTE: `TRCRSMatrix` should be a `psi::CRSMatrix`-like type from Range group.
   */
  template< typename TRCRSMatrix, typename TExecSpace=Kokkos::DefaultExecutionSpace >
  inline TRCRSMatrix
  create_range_identity_matrix( typename TRCRSMatrix::ordinal_type n, TExecSpace space={} )
  {
    typedef TRCRSMatrix rcrsmatrix_t;
    typedef typename rcrsmatrix_t::spec_type rcrs_spec_type;
    typedef typename crs_matrix::Group< rcrs_spec_type >::type group_type;

    auto i_entries = rcrsmatrix_t::make_entries_device_view( space );
    auto i_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );

    create_range_identity_matrix( i_rowmap, i_entries, n );

    if constexpr ( std::is_same< group_type, crs_matrix::RangeGroup >::value )
      return rcrsmatrix_t( n, i_entries, i_rowmap, n /* nnz */ );
    else
      return rcrsmatrix_t( n, i_entries, i_rowmap );
  }

  /**
   *  @brief  Create a random square matrix (on-host) [SLOW].
   *
   *  @param  n     order of the output matrix
   *  @param  nnz   number of non-zero values
   *  @param  lower lower bound of value range (inclusive)
   *  @param  upper upper bound of value range (exclusive)
   *
   *  @return a `KokkosSparse::CrsMatrix`-like square matrix of order `n` with
   *  `nnz` number of non-zero random values in range [lower, upper).
   *
   *  NOTE: Use `create_random_matrix` instead.
   *  NOTE: The output matrix is constructed on **host** memory space, then it
   *  is transferred (deep-copied) to device memory space. It is suitable for
   *  small matrices.
   *
   *  NOTE: The final matrix is on device memory and host mirrors get
   *  deallocated on return.
   */
  template< typename TXCRSMatrix >
  inline TXCRSMatrix
  create_random_matrix_on_host( typename TXCRSMatrix::ordinal_type n,
                                typename TXCRSMatrix::size_type nnz,
                                typename TXCRSMatrix::value_type lower=std::numeric_limits< typename TXCRSMatrix::value_type >::min(),
                                typename TXCRSMatrix::value_type upper=std::numeric_limits< typename TXCRSMatrix::value_type >::max() )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename xcrsmatrix_t::value_type value_type;
    typedef typename xcrsmatrix_t::size_type size_type;
    typedef typename xcrsmatrix_t::values_type::non_const_type values_t;
    typedef typename xcrsmatrix_t::row_map_type::non_const_type row_map_t;
    typedef typename xcrsmatrix_t::index_type::non_const_type entries_t;

    assert( n > 1 && nnz > 0 && ( nnz / n ) <= static_cast< size_type >( n ) );

    values_t a_values( Kokkos::ViewAllocateWithoutInitializing( "R" ), nnz );
    row_map_t a_row_map( Kokkos::ViewAllocateWithoutInitializing( "rowmap" ),
                         n + 1 );
    entries_t a_entries( Kokkos::ViewAllocateWithoutInitializing( "entries" ),
                         nnz );

    auto h_a_entries = Kokkos::create_mirror_view( a_entries );
    auto h_a_values = Kokkos::create_mirror_view( a_values );
    auto h_a_row_map = Kokkos::create_mirror_view( a_row_map );

    // Zero initialisation: the rest will be initialized later
    h_a_row_map( 0 ) = 0;

    Kokkos::parallel_for(
        "psi::crs_matrix::::create_random_matrix_on_host::random_values",
        Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >( 0, nnz ),
        KOKKOS_LAMBDA ( const uint64_t i ) {
          value_type v = 0;
          while ( v == 0 ) v = psi::random::random_integer( lower, upper + 1 );
          h_a_values( i ) = v;
        } );

    {
      // Distributing nnz values into rows
      std::size_t i = 0;
      while ( i < nnz ) {
        auto idx = psi::random::random_index( n );
        do {
          if ( h_a_row_map( idx + 1 ) < n ) {
            ++h_a_row_map( idx + 1 );
            ++i;
            break;
          }
          idx = ( idx + 1 ) % n;
        } while( true );
      }

      //for ( i = 1; i < n; ++i ) h_a_row_map( i + 1 ) += h_a_row_map( i );
      Kokkos::parallel_scan(
          "psi::crs_matrix::::create_random_matrix_on_host::compute_row_map",
          Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >( 0, n ),
          KOKKOS_LAMBDA ( const int i, size_type& update, const bool final ) {
            // Load old value in case we update it before accumulating
            const size_type val_ip1 = h_a_row_map( i + 1 );
            update += val_ip1;
            if ( final )
              h_a_row_map( i + 1 )
                  = update;  // only update array on final pass
          } );
    }

    Kokkos::parallel_for(
        "psi::crs_matrix::::create_random_matrix_on_host::random_entries",
        Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >( 0, n ),
        KOKKOS_LAMBDA ( const uint64_t i ) {
          auto l = h_a_row_map( i );
          auto u = h_a_row_map( i + 1 );
          auto begin = h_a_entries.data() + l;
          auto end = h_a_entries.data() + u;
          std::sample( psi::RangeIterator< decltype( n ) >{ 0 },
                       psi::RangeIterator{ n }, begin, u - l,
                       psi::random::gen );
          std::sort( begin, end );
        } );

    Kokkos::deep_copy( a_entries, h_a_entries );
    Kokkos::deep_copy( a_values, h_a_values );
    Kokkos::deep_copy( a_row_map, h_a_row_map );

    assert( h_a_row_map( n ) == nnz );

    return xcrsmatrix_t( "Random Matrix", n, n, nnz, a_values, a_row_map, a_entries );
  }

  /**
   *  @brief  Create a random square matrix (on-device) [FAST].
   *
   *  @param  n     order of the output matrix
   *  @param  nnz   number of non-zero values
   *  @param  lower lower bound of value range (inclusive)
   *  @param  upper upper bound of value range (exclusive)
   *
   *  @return a `KokkosSparse::CrsMatrix`-like square matrix of order `n` with
   * `nnz` number of non-zero random values in range [lower, upper).
   *
   *  This function is usually much faster than `create_random_matrix_on_host`.
   *
   *  NOTE: The output matrix is constructed on device memory space. It does
   *    not performs any deep copy between device and host.
   */
  template< typename TXCRSMatrix,
            typename THostSpace=Kokkos::DefaultHostExecutionSpace >
  inline TXCRSMatrix
  create_random_matrix(
      typename TXCRSMatrix::ordinal_type n,
      typename TXCRSMatrix::size_type nnz,
      typename TXCRSMatrix::value_type lower
      = std::numeric_limits< typename TXCRSMatrix::value_type >::min(),
      typename TXCRSMatrix::value_type upper
      = std::numeric_limits< typename TXCRSMatrix::value_type >::max() )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename TXCRSMatrix::execution_space execution_space;
    typedef typename xcrsmatrix_t::value_type value_type;
    typedef typename xcrsmatrix_t::ordinal_type ordinal_type;
    typedef typename xcrsmatrix_t::size_type size_type;
    typedef typename xcrsmatrix_t::values_type::non_const_type values_t;
    typedef typename xcrsmatrix_t::row_map_type::non_const_type row_map_t;
    typedef typename xcrsmatrix_t::index_type::non_const_type entries_t;
    typedef Kokkos::RangePolicy< execution_space > range_policy_t;
    typedef typename Kokkos::Random_XorShift64_Pool< execution_space > random_pool_t;
    typedef typename random_pool_t::generator_type generator_t;

    assert( n > 1 && nnz > 0 && ( nnz / n ) <= static_cast< size_type >( n ) );

    // NOTE: (nnz/n) should fit in `ordinal_type` (asserted above) as
    //       `sizeof( size_type ) >= sizeof( ordinal_type )`
    // NOTE: `ordinal_type` to match with the signness of `n`
    ordinal_type nnz_per_row = nnz / n;

    values_t r_values( Kokkos::ViewAllocateWithoutInitializing( "R" ), nnz );
    row_map_t r_row_map( "rowmap", n + 1 );
    entries_t r_entries( Kokkos::ViewAllocateWithoutInitializing( "entries" ),
                         nnz );

    random_pool_t random_pool( psi::random::rd() );

    Kokkos::parallel_for( "psi::crs_matrix::create_random_matrix::random_values",
                          range_policy_t( 0, nnz ),
                          KOKKOS_LAMBDA ( const uint64_t i ) {
                            generator_t generator = random_pool.get_state();
                            value_type v = 0;
                            while ( v == 0 ) {
                              v = Kokkos::rand< generator_t,
                                                value_type >::draw( generator,
                                                                    lower,
                                                                    upper );
                            }
                            r_values( i ) = v;
                            random_pool.free_state( generator );
                          } );

    if ( nnz_per_row == n ) {  // if nnz = n*n (n*n may be really large to fit in an integer)
      Kokkos::parallel_for(
          "psi::crs_matrix::create_random_matrix::fill_nnz",
          range_policy_t( 0, n ),
          KOKKOS_LAMBDA ( const uint64_t i ) { r_row_map( i + 1 ) = n; } );
    }
    else {
      auto d_nnz = nnz;
      if ( nnz_per_row > n / 2 ) {  // in this case, distribute the zero values is cheaper.
        d_nnz = ( n - ( nnz_per_row ) ) * n - ( nnz % n );
      }

      Kokkos::parallel_for(
          "psi::crs_matrix::create_random_matrix::distribute_nnz",
          range_policy_t( 0, d_nnz ), KOKKOS_LAMBDA ( const uint64_t ) {
            generator_t generator = random_pool.get_state();

            ordinal_type idx =
              Kokkos::rand< generator_t, ordinal_type >::draw( generator, n );

            bool exchanged = false;
            do {
              auto ptr = &r_row_map( idx + 1 );
              auto value = Kokkos::atomic_load( ptr );
              while ( value < static_cast< size_type >( n ) ) {
                exchanged = Kokkos::atomic_compare_exchange_strong(
                    ptr, value, value + 1 );
                if ( exchanged ) break;
                value = Kokkos::atomic_load( ptr );
              }
              idx = ( idx + 1 ) % n;
            } while ( !exchanged );

            random_pool.free_state( generator );
          } );

      if ( d_nnz != nnz ) {
        Kokkos::parallel_for(
            "psi::crs_matrix::create_random_matrix::reverse_nnz_dist",
            range_policy_t( 0, n ), KOKKOS_LAMBDA ( const uint64_t i ) {
              r_row_map( i + 1 ) = n - r_row_map( i + 1 );
            } );
      }
    }

    Kokkos::parallel_scan(
        "psi::crs_matrix::create_random_matrix::compute_row_map",
        range_policy_t( 0, n ),
        KOKKOS_LAMBDA ( const int i, size_type& partial_sum,
                        const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type value = r_row_map( i + 1 );
          partial_sum += value;
          // only update array on final pass
          if ( final ) r_row_map( i + 1 ) = partial_sum;
          if ( i == 0 ) r_row_map( 0 ) = 0;
        } );

    Kokkos::parallel_for(
        "psi::crs_matrix::create_random_matrix::random_entries",
        range_policy_t( 0, n ), KOKKOS_LAMBDA ( const uint64_t i ) {
          auto l = r_row_map( i );
          auto u = r_row_map( i + 1 );
          auto begin = r_entries.data() + l;
          assert( l <= u );
          ordinal_type k = u - l;
          if ( k != 0 ) {
            generator_t generator = random_pool.get_state();

            // Reservoir sampling algorithm
            ordinal_type j = 0;
            for ( j = 0; j < k; ++j ) *( begin + j ) = j;
            for ( ; j < n; ++j ) {
              ordinal_type r = Kokkos::rand< generator_t, ordinal_type >::draw(
                  generator, j );
              if ( r < k ) *( begin + r ) = j;
            }

            random_pool.free_state( generator );
          }
        } );

    auto func = SortEntriesFunctor( r_row_map, r_entries );

    Kokkos::parallel_for(
        "psi::crs_matrix::create_random_matrix::sort_entries",
        func.policy( n ), func );

    return xcrsmatrix_t( "Random Matrix", n, n, nnz, r_values, r_row_map,
                         r_entries );
  }

  /**
   *  @brief  Create a random binary `KokkosSparse::CrsMatrix`.
   */
  template< typename TXCRSMatrix >
  inline TXCRSMatrix
  create_random_binary_matrix( typename TXCRSMatrix::ordinal_type n,
                               typename TXCRSMatrix::size_type nnz )
  {
    return create_random_matrix< TXCRSMatrix >( n, nnz, 1, 2 );
  }

  /**
   *  @brief  Create a random binary matrix and return in both KokkosKernels
   *          CRS and PSI Range CRS formats.
   */
  template< typename TXCRSMatrix, typename TRCRSMatrix >
  inline TXCRSMatrix
  create_random_binary_matrix( typename TXCRSMatrix::ordinal_type n,
                               typename TXCRSMatrix::size_type nnz,
                               TRCRSMatrix& range_crs )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename TXCRSMatrix::HostMirror xcrsmatrix_host_t;
    typedef TRCRSMatrix range_crsmatrix_t;

    typename xcrsmatrix_t::index_type::non_const_type::HostMirror h_r_entries;
    typename xcrsmatrix_t::values_type::non_const_type::HostMirror h_r_values;
    typename xcrsmatrix_t::row_map_type::non_const_type::HostMirror h_r_row_map;

    auto crs = create_random_binary_matrix< xcrsmatrix_t >( n, nnz );

    h_r_values = Kokkos::create_mirror_view( crs.values );
    h_r_row_map = Kokkos::create_mirror_view( crs.graph.row_map );
    h_r_entries = Kokkos::create_mirror_view( crs.graph.entries );

    Kokkos::deep_copy( h_r_values, crs.values );
    Kokkos::deep_copy( h_r_row_map, crs.graph.row_map );
    Kokkos::deep_copy( h_r_entries, crs.graph.entries );

    const xcrsmatrix_host_t h_crs( "R host copy", n, n, nnz, h_r_values, h_r_row_map, h_r_entries );
    range_crs = range_crsmatrix_t( h_crs );

    return crs;
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC >
  inline void
  range_spadd_symbolic( const THandle&,
                        TRowMapDeviceViewA a_rowmap,
                        TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap,
                        TEntriesDeviceViewB b_entries,
                        TRowMapDeviceViewC& c_rowmap )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::parallel_for(
        "psi::crs_matrix::range_spadd_symbolic::count_row_nnz",
        policy_type( 0, a_nrows ), KOKKOS_LAMBDA ( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          auto b_idx = b_rowmap( row );
          auto b_end = b_rowmap( row + 1 );

          assert( ( a_end - a_idx ) % 2 == 0 );
          assert( ( b_end - b_idx ) % 2 == 0 );

          size_type count = 0;
          if ( a_idx < a_end && b_idx < b_end ) {
            ordinal_type lo;
            ordinal_type hi;
            if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
              lo = a_entries( a_idx );
              hi = a_entries( a_idx + 1 );
              a_idx += 2;
            }
            else {
              lo = b_entries( b_idx );
              hi = b_entries( b_idx + 1 );
              b_idx += 2;
            }

            while ( a_idx < a_end && b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }
              else {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = PSI_MACRO_MIN( lo, rs );
                hi = PSI_MACRO_MAX( hi, re );
                continue;
              }
              lo = rs;
              hi = re;
              count += 2;
            }

            while ( a_idx < a_end || b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_idx == a_end ) {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }
              else {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = PSI_MACRO_MIN( lo, rs );
                hi = PSI_MACRO_MAX( hi, re );
                continue;
              }
              count += 2;
              break;
            }
            count += 2;  // last one from the previous loop
          }

          count += ( b_end - b_idx ) + ( a_end - a_idx );  // the rest

          c_rowmap( row + 1 ) = count;
          if ( row == 0 ) c_rowmap( 0 ) = 0;
        } );

    Kokkos::parallel_scan(
        "psi::crs_matrix::range_spadd_symbolic::computing_row_map_c",
        policy_type( 0, a_nrows ),
        KOKKOS_LAMBDA ( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC >
  inline void
  range_spadd_numeric( const THandle&,
                       TRowMapDeviceViewA a_rowmap,
                       TEntriesDeviceViewA a_entries,
                       TRowMapDeviceViewB b_rowmap,
                       TEntriesDeviceViewB b_entries,
                       TRowMapDeviceViewC c_rowmap,
                       TEntriesDeviceViewC& c_entries )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::parallel_for(
        "psi::crs_matrix::range_spadd_numeric::count_row_nnz",
        policy_type( 0, a_nrows ), KOKKOS_LAMBDA ( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          auto b_idx = b_rowmap( row );
          auto b_end = b_rowmap( row + 1 );

          assert( ( a_end - a_idx ) % 2 == 0 );
          assert( ( b_end - b_idx ) % 2 == 0 );

          size_type c_idx = c_rowmap( row );
          if ( a_idx < a_end && b_idx < b_end ) {
            ordinal_type lo;
            ordinal_type hi;
            if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
              lo = a_entries( a_idx );
              hi = a_entries( a_idx + 1 );
              a_idx += 2;
            }
            else {
              lo = b_entries( b_idx );
              hi = b_entries( b_idx + 1 );
              b_idx += 2;
            }

            while ( a_idx < a_end && b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }
              else {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = PSI_MACRO_MIN( lo, rs );
                hi = PSI_MACRO_MAX( hi, re );
                continue;
              }
              c_entries( c_idx++ ) = lo;
              c_entries( c_idx++ ) = hi;
              lo = rs;
              hi = re;
            }

            while ( a_idx < a_end || b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_idx == a_end ) {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }
              else {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = PSI_MACRO_MIN( lo, rs );
                hi = PSI_MACRO_MAX( hi, re );
                continue;
              }
              c_entries( c_idx++ ) = lo;
              c_entries( c_idx++ ) = hi;
              lo = rs;
              hi = re;
              break;
            }
            c_entries( c_idx++ ) = lo;
            c_entries( c_idx++ ) = hi;
          }
          for ( ; a_idx < a_end; ++a_idx )
            c_entries( c_idx++ ) = a_entries( a_idx );
          for ( ; b_idx < b_end; ++b_idx )
            c_entries( c_idx++ ) = b_entries( b_idx );
        } );
  }

  /**
   *  @brief Symbolic phase of computing matrix c as the product of a and b
   *  (ThreadRangePartition-BTreeAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` is allocated on device and is of
   *        size `a.numRows()+1`.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC >
  inline void
  _range_spgemm_symbolic( const THandle&,
                          TRowMapDeviceViewA a_rowmap,
                          TEntriesDeviceViewA a_entries,
                          TRowMapDeviceViewB b_rowmap,
                          TEntriesDeviceViewB b_entries,
                          TRowMapDeviceViewC& c_rowmap,
                          ThreadRangePartition, BTreeAccumulator )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;
    typedef std::map< ordinal_type, ordinal_type > btree_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::Experimental::UniqueToken< execution_space > token;
    std::vector< btree_type > maps( token.size() );
    auto maps_ptr = &maps;

    Kokkos::parallel_for(
        "psi::crs_matrix::range_spgemm_symbolic::count_row_nnz",
        policy_type( 0, a_nrows ), KOKKOS_LAMBDA ( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          int id = token.acquire();
          auto& acc = maps_ptr->operator[]( id );
          for ( ; a_idx != a_end; a_idx += 2 ) {
            auto b_idx = b_rowmap( a_entries( a_idx ) );
            auto b_end = b_rowmap( a_entries( a_idx + 1 ) + 1 );
            for ( ; b_idx != b_end; b_idx += 2 ) {
              auto old = acc[ b_entries( b_idx ) ];
              acc[ b_entries( b_idx ) ]
                  = PSI_MACRO_MAX( old, b_entries( b_idx + 1 ) );
            }
          }

          ordinal_type count = 0;
          if ( !acc.empty() ) {
            auto it = acc.begin();
            ordinal_type lo = it->first;
            ordinal_type hi = it->second;
            ++it;
            for ( ; it != acc.end(); ++it ) {
              if ( it->first <= hi + 1 ) {  // merge
                lo = PSI_MACRO_MIN( lo, it->first );
                hi = PSI_MACRO_MAX( hi, it->second );
                continue;
              }
              lo = it->first;
              hi = it->second;
              count += 2;
            }
            count += 2;
            acc.clear();
          }

          token.release(id);

          c_rowmap( row + 1 ) = count;
          if ( row == 0 ) c_rowmap( 0 ) = 0;
        } );

    Kokkos::parallel_scan(
        "psi::crs_matrix::range_spgemm_symbolic::computing_row_map_c",
        policy_type( 0, a_nrows ),
        KOKKOS_LAMBDA ( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );
  }

  /**
   *  @brief Numeric phase of computing matrix c as the product of a and b
   *         (ThreadRangePartition-BTreeAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` and `c_entries` are allocated on
   *        device with sufficient space.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC >
  inline void
  _range_spgemm_numeric( const THandle&,
                        TRowMapDeviceViewA a_rowmap,
                        TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap,
                        TEntriesDeviceViewB b_entries,
                        TRowMapDeviceViewC c_rowmap,
                        TEntriesDeviceViewC& c_entries,
                        ThreadRangePartition, BTreeAccumulator )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_entries_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;
    typedef std::map< ordinal_type, ordinal_type > btree_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::Experimental::UniqueToken< execution_space > token;
    std::vector< btree_type > maps( token.size() );
    auto maps_ptr = &maps;

    Kokkos::parallel_for(
        "psi::crs_matrix::range_spgemm_numeric::compute_numeric",
        policy_type( 0, a_nrows ), KOKKOS_LAMBDA ( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          int id = token.acquire();
          auto& acc = maps_ptr->operator[]( id );
          for ( ; a_idx < a_end; a_idx += 2 ) {
            auto b_idx = b_rowmap( a_entries( a_idx ) );
            auto b_end = b_rowmap( a_entries( a_idx + 1 ) + 1 );
            for ( ; b_idx != b_end; b_idx += 2 ) {
              auto old = acc[ b_entries( b_idx ) ];
              acc[ b_entries( b_idx ) ]
                  = PSI_MACRO_MAX( old, b_entries( b_idx + 1 ) );
            }
          }

          if ( !acc.empty() ) {
            size_type c_idx = c_rowmap( row );
            auto it = acc.begin();
            ordinal_type lo = it->first;
            ordinal_type hi = it->second;
            ++it;
            for ( ; it != acc.end(); ++it ) {
              if ( it->first <= hi + 1 ) {  // merge
                lo = PSI_MACRO_MIN( lo, it->first );
                hi = PSI_MACRO_MAX( hi, it->second );
                continue;
              }
              c_entries( c_idx++ ) = lo;
              c_entries( c_idx++ ) = hi;
              lo = it->first;
              hi = it->second;
            }
            c_entries( c_idx++ ) = lo;
            c_entries( c_idx++ ) = hi;
            acc.clear();
          }
          token.release( id );
        } );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline void
  range_spgemm_symbolic( const THandle& handle,
                        TRowMapDeviceViewA a_rowmap,
                        TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap,
                        TEntriesDeviceViewB b_entries,
                        TRowMapDeviceViewC& c_rowmap, TSparseConfig = {} )
  {
    typedef typename TEntriesDeviceViewA::non_const_value_type ordinal_type;
    typedef typename TSparseConfig::partition_type partition_type;
    typedef typename TSparseConfig::accumulator_type accumulator_type;

    assert( handle.a_ncols == ( b_rowmap.extent( 0 ) - 1 ) );
    ASSERT( handle.a_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );
    ASSERT( handle.b_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );

    _range_spgemm_symbolic( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                            c_rowmap, partition_type(), accumulator_type() );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline void
  range_spgemm_numeric( const THandle& handle,
                        TRowMapDeviceViewA a_rowmap, TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap, TEntriesDeviceViewB b_entries,
                        TRowMapDeviceViewC c_rowmap, TEntriesDeviceViewC& c_entries,
                        TSparseConfig={} )
  {
    typedef typename TEntriesDeviceViewC::value_type ordinal_type;
    typedef typename TSparseConfig::partition_type partition_type;
    typedef typename TSparseConfig::accumulator_type accumulator_type;

    assert( handle.a_ncols == ( b_rowmap.extent( 0 ) - 1 ) );
    ASSERT( handle.a_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );
    ASSERT( handle.b_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );

    _range_spgemm_numeric( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                           c_rowmap, c_entries, partition_type(),
                           accumulator_type() );
  }

  /**
   *  @brief Computing matrix c as the product of a and b.
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline void
  range_spgemm( const THandle& handle,
                TRowMapDeviceViewA a_rowmap, TEntriesDeviceViewA a_entries,
                TRowMapDeviceViewB b_rowmap, TEntriesDeviceViewB b_entries,
                TRowMapDeviceViewC& c_rowmap, TEntriesDeviceViewC& c_entries,
                TSparseConfig config={} )
  {
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;

    assert( handle.a_ncols == ( b_rowmap.extent( 0 ) - 1 ) );
    ASSERT( handle.a_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );
    ASSERT( handle.b_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );

    ordinal_type n = a_rowmap.extent( 0 ) - 1;
    c_rowmap = c_row_map_type(
        Kokkos::ViewAllocateWithoutInitializing( "c_rowmap" ),
        a_rowmap.extent( 0 ) );

#ifdef PSI_STATS
    Kokkos::Timer timer;
#endif

    range_spgemm_symbolic( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                           c_rowmap, config );

#ifdef PSI_STATS
    double d = timer.seconds();
    std::cout << "psi::range_spgemm_symbolic time: " << d * 1000 << "ms"
              << std::endl;
#endif

    size_type c_rnnz;
    Kokkos::deep_copy( c_rnnz, Kokkos::subview( c_rowmap, n ) );
    c_entries = c_entries_type( Kokkos::ViewAllocateWithoutInitializing( "C" ),
                                c_rnnz );

#ifdef PSI_STATS
    timer.reset();
#endif

    range_spgemm_numeric( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                          c_rowmap, c_entries, config );

#ifdef PSI_STATS
    d = timer.seconds();
    std::cout << "psi::range_spgemm_numeric time: " << d * 1000 << "ms"
              << std::endl;
#endif
  }

  template< typename TRCRSMatrix,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline TRCRSMatrix
  range_spgemm( TRCRSMatrix const& a, TRCRSMatrix const& b,
                TSparseConfig config={} )
  {
    typedef TRCRSMatrix range_crsmatrix_t;
    typedef TSparseConfig config_type;
    typedef typename range_crsmatrix_t::ordinal_type ordinal_type;
    typedef typename config_type::execution_space execution_space;

    assert( a.numCols() == b.numRows() );
    ASSERT( a.numCols() <= std::numeric_limits< ordinal_type >::max() - 1 );
    ASSERT( b.numCols() <= std::numeric_limits< ordinal_type >::max() - 1 );

    execution_space space{};

    auto a_entries = a.entries_device_view( space );
    auto a_rowmap = a.rowmap_device_view( space );
    auto b_entries = b.entries_device_view( space );
    auto b_rowmap = b.rowmap_device_view( space );

    auto c_entries = range_crsmatrix_t::make_entries_device_view( space );
    auto c_rowmap = range_crsmatrix_t::make_rowmap_device_view( space );

    SparseRangeHandle handle( a, b );

    range_spgemm( handle, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap,
                  c_entries, config );

    // FIXME: since entries and rowmap arrays of range CRS is not a view, there
    // would be an extra copy here and the `c_entries` and `c_rowmap` cannot be
    // moved when the views are on the same memory space (the RCRS ctor does call
    // `deep_copy`).
    return TRCRSMatrix( b.numCols(), c_entries, c_rowmap );
  }

  template< typename TRCRSMatrix,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline TRCRSMatrix
  range_power( TRCRSMatrix const& a, unsigned int k, TSparseConfig config={} )
  {
    typedef TRCRSMatrix rcrsmatrix_t;
    typedef TSparseConfig config_type;
    typedef typename config_type::execution_space execution_space;

    execution_space space{};

    auto c_entries = rcrsmatrix_t::make_entries_device_view( space );
    auto c_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );

    create_range_identity_matrix( c_rowmap, c_entries, a.numRows() );

    auto a2n_entries = a.entries_device_view( space );
    auto a2n_rowmap = a.rowmap_device_view( space );

    auto tmp_entries = rcrsmatrix_t::make_entries_device_view( space );
    auto tmp_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );

    SparseRangeHandle handle( a, a );

    while ( true ) {
      if ( k & 1 ) {
        range_spgemm( handle, c_rowmap, c_entries, a2n_rowmap, a2n_entries,
                      tmp_rowmap, tmp_entries, config );
        c_entries = tmp_entries;
        c_rowmap = tmp_rowmap;
      }

      k = k >> 1;

      if ( k != 0 ) {
        range_spgemm( handle, a2n_rowmap, a2n_entries, a2n_rowmap, a2n_entries,
                      tmp_rowmap, tmp_entries, config );
        a2n_entries = tmp_entries;
        a2n_rowmap = tmp_rowmap;
      }
      else break;
    }

    return TRCRSMatrix( a.numCols(), c_entries, c_rowmap );
  }
}  // namespace psi

#endif  // PSI_RANGE_SPARSE_HPP_
