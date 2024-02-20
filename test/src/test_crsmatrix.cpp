/**
 *    @file  test_crsmatrix.cpp
 *   @brief  CRS matrix module test scenarios.
 *
 *  Contains test cases for CRS matrix module.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Nov 12, 2020  23:20
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <numeric>
#include <algorithm>

#include <psi/crs_matrix.hpp>
#include <pairg/spgemm_utility.hpp>

#include "test_base.hpp"


using namespace psi;

namespace test_util {
  template< typename TCRSMatrix >
  inline auto
  numRows( TCRSMatrix const& matrix )
  {
    return matrix.numRows();
  }

  template< typename TValue, std::size_t NRows, std::size_t NCols >
  constexpr inline auto
  numRows( std::array< std::array< TValue, NCols >, NRows > const& matrix )
  {
    return NRows;
  }

  template< typename TCRSMatrix >
  inline auto
  numCols( TCRSMatrix const& matrix )
  {
    return matrix.numCols();
  }

  template< typename TValue, std::size_t NRows, std::size_t NCols >
  constexpr inline auto
  numCols( std::array< std::array< TValue, NCols >, NRows > const& matrix )
  {
    return NCols;
  }

  /* NOTE: The `matrix` cannot be const reference because of Buffered specialisations. */
  template< typename TCRSMatrix >
  inline auto
  access( TCRSMatrix& matrix, std::size_t i, std::size_t j )
  {
    return matrix( i, j );
  }

  template< typename TValue, std::size_t NRows, std::size_t NCols >
  inline auto
  access( std::array< std::array< TValue, NCols >, NRows >& matrix,
          std::size_t i, std::size_t j )
  {
    return matrix[i][j];
  }

  /* NOTE: The `matrix` cannot be const reference because of Buffered specialisations. */
  template< typename TValue, std::size_t NRows, std::size_t NCols >
  inline std::size_t
  get_nnz( std::array< std::array< TValue, NCols >, NRows >& matrix )
  {
    using host_space = Kokkos::DefaultHostExecutionSpace;
    using rank_type = Kokkos::Rank< 2 >;

    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );
    auto mat_ptr = &matrix;

    std::size_t nnz = 0;
    Kokkos::parallel_reduce(
        "psi::test_crsmatrix::compute_nnz",
        Kokkos::MDRangePolicy< host_space, rank_type >( { 0, 0 }, { nrows, ncols } ),
        KOKKOS_LAMBDA ( const uint64_t i, const uint64_t j, std::size_t& l_nnz ) {
          if ( access( *mat_ptr, i, j ) == 1 ) l_nnz += 1;
        }, nnz );

    return nnz;
  }

  template< typename TCRSMatrix >
  inline std::size_t
  get_nnz( TCRSMatrix& matrix )
  {
    return matrix.nnz();
  }

  template< typename TMatrix >
  inline void
  zero_matrix( TMatrix& matrix )
  {
    using host_space = Kokkos::DefaultHostExecutionSpace;
    using rank_type = Kokkos::Rank< 2 >;

    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );
    auto mat_ptr = &matrix;

    Kokkos::parallel_for(
        "psi::test_crsmatrix::zero_matrix",
        Kokkos::MDRangePolicy< host_space, rank_type >( { 0, 0 }, { nrows, ncols } ),
        KOKKOS_LAMBDA ( const uint64_t i, const uint64_t j ) {
          ( *mat_ptr )[i][j] = 0;
        } );
  }

  template< typename TMatrix >
  inline void
  random_matrix( TMatrix& matrix, std::size_t nnz )
  {
    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );

    assert( nnz <= nrows * ncols );

    // initialise the matrix by zero
    zero_matrix( matrix );

    std::size_t i = 0;
    while ( i < nnz ) {
      auto r = random::random_index( nrows, rnd::get_rgn() );
      auto c = random::random_index( ncols, rnd::get_rgn() );
      if ( matrix[r][c] == 0 ) {
        matrix[r][c] = 1;
        ++i;
      }
    }
  }

  template< typename TMatrix >
  inline void
  random_matrix_ranged( TMatrix& matrix, std::size_t nnz )
  {
    static const std::size_t MIN_NOF_FRAGS = 1;
    static const std::size_t MAX_NOF_FRAGS = 4;

    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );

    assert( nnz <= nrows * ncols );

    // initialise the matrix by zero
    zero_matrix( matrix );

    std::size_t base_zp = nnz / nrows;     // minimum number of non-zero value per row
    std::size_t remainders = nnz % nrows;  // remainder of the division scattered across rows
    std::size_t nfrags = std::max( MIN_NOF_FRAGS, std::min( MAX_NOF_FRAGS, base_zp ) );
    std::size_t filled = 0;
    std::vector< std::size_t > lens;
    std::vector< std::size_t > spaces;
    std::vector< std::size_t > therange( ncols );

    std::iota( therange.begin(), therange.end(), 1 );

    auto partition = [&]( std::vector< std::size_t >& out, std::size_t len,
                            std::size_t nfrag, bool can_empty = false ) {
      assert( 1 <= len && len <= ncols );
      if ( can_empty ) {
        std::generate_n( out.begin(), nfrag - 1, [&](){ return random::random_index( len ); } );
      }
      else {
        std::sample( therange.begin(), therange.begin() + len - 1, out.begin(), nfrag - 1, rnd::get_rgn() );
      }
      out.back() = len;
      std::sort( out.begin(), out.end() );
      std::adjacent_difference( out.begin(), out.end(), out.begin() );
    };

    for ( std::size_t i = 0; i < nrows; ++i ) {
      //if ( nnz < filled + base_zp ) base_zp = nnz - filled;
      assert( filled <= nnz );
      if ( filled == nnz ) break;
      auto c_nfrags = random::random_integer( 1ul, nfrags, rnd::get_rgn() );  // number of fragments in the current row
      auto nnzp = base_zp;  // number of non-zero values in the current row
      if ( remainders ) {
        ++nnzp;
        --remainders;
      }

      lens.resize( c_nfrags );
      partition( lens, nnzp, c_nfrags );
      spaces.resize( c_nfrags );
      partition( spaces, ncols - nnzp, c_nfrags + 1, true );

      std::size_t j = 0;
      auto cspace = spaces.begin();
      for ( auto it = lens.begin(); it != lens.end(); ++it ) {
        j += *cspace++;
        for ( std::size_t k = 0; k < *it; ++k ) {
          matrix[i][j] = 1;
          ++j;
        }
      }

      filled += nnzp;
    }
  }

  /* NOTE: The `matrix` cannot be const reference because of Buffered specialisations. */
  template< typename TMatrix,
            typename TExternalMatrix = KokkosSparse::CrsMatrix<
                char, int, Kokkos::DefaultHostExecutionSpace > >
  inline TExternalMatrix
  to_external_crs( TMatrix& matrix, std::size_t nnz=0 )
  {
    typedef Kokkos::DefaultHostExecutionSpace host_space;
    typedef Kokkos::TeamPolicy< host_space > policy_type;
    typedef typename policy_type::member_type member_type;

    typedef TExternalMatrix xcrsmat_type;
    typedef typename xcrsmat_type::staticcrsgraph_type staticcrsgraph_type;
    typedef typename staticcrsgraph_type::size_type size_type;
    typedef typename staticcrsgraph_type::entries_type::non_const_type entries_type;
    typedef typename staticcrsgraph_type::row_map_type::non_const_type row_map_type;
    typedef typename xcrsmat_type::values_type::non_const_type values_type;

    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );
    if ( nnz == 0 ) nnz = get_nnz( matrix );

    entries_type entries( "entries", nnz );
    values_type values( "values", nnz );
    row_map_type rowmap( "rowmap", nrows + 1 );

    Kokkos::parallel_for(
        "psi::test_crsmatrix::to_external_crs::set_values",
        Kokkos::RangePolicy< host_space >( 0, nnz ),
        KOKKOS_LAMBDA ( const uint64_t i ) {
          values( i ) = 1;  // boolean
        } );

    auto mat_ptr = &matrix;

    Kokkos::parallel_for(
        "psi::test_crsmatrix::to_external_crs::count_row_nnz",
        policy_type( nrows, Kokkos::AUTO ),
        KOKKOS_LAMBDA ( const member_type& tm ) {
          auto i = tm.league_rank();
          size_type row_nnz = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange( tm, ncols ),
              [=]( const uint64_t j, size_type& l_row_nnz ) {
                if ( ( *mat_ptr )[ i ][ j ] == 1 ) ++l_row_nnz;
              }, row_nnz );

          Kokkos::single( Kokkos::PerTeam( tm ),
                          [=]() { rowmap( i + 1 ) = row_nnz; } );
        } );

    rowmap( 0 ) = 0;
    Kokkos::parallel_scan(
        "psi::test_crsmatrix::to_external_crs::computing_rowmap",
        Kokkos::RangePolicy< host_space >( 0, nrows ),
        KOKKOS_LAMBDA ( const uint64_t i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            rowmap( i + 1 ) = update;  // only update array on final pass
        } );

    Kokkos::parallel_for(
        "psi::test_crsmatrix::to_external_crs::fill_entries",
        Kokkos::RangePolicy< host_space >( 0, nrows ),
        KOKKOS_LAMBDA ( const uint64_t i ) {
          auto eidx = rowmap( i );
          for ( size_type j = 0; j < ncols; ++j ) {
            if ( ( *mat_ptr )[ i ][ j ] == 1 ) entries( eidx++ ) = j;
          }
          assert( eidx == rowmap( i + 1 ) );
        } );

    assert( rowmap( nrows ) == nnz );

    return xcrsmat_type( "matrix", nrows, ncols, nnz, values, rowmap, entries );
  }

  template< typename TMatrix1, typename TMatrix2 >
  inline void
  fill_block( TMatrix1& matrix, TMatrix2 const& block,
              unsigned int& si, unsigned int& sj )
  {
    using host_space = Kokkos::DefaultHostExecutionSpace;
    using rank_type = Kokkos::Rank< 2 >;

    auto mat_ptr = &matrix;
    auto block_ptr = &block;

    Kokkos::parallel_for(
        "psi::test_crsmatrix::fill_block",
        Kokkos::MDRangePolicy< host_space, rank_type >(
            { 0, 0 }, { block.size(), block[ 0 ].size() } ),
        KOKKOS_LAMBDA ( const uint64_t i, const uint64_t j ) {
          ( *mat_ptr )[si + i][sj + j] = ( *block_ptr )[i][j];
        } );

    si += numRows( block );
    sj += numCols( block );
  }

  template< typename TMatrix1, typename TMatrix2 >
  inline void
  is_identical( TMatrix1& matrix1, TMatrix2& matrix2 )
  {
    using host_space = Kokkos::DefaultHostExecutionSpace;
    using rank_type = Kokkos::Rank< 2 >;

    REQUIRE( numRows( matrix1 ) != 0 );
    REQUIRE( numCols( matrix1 ) != 0 );
    REQUIRE( numRows( matrix1 ) == numRows( matrix2 ) );
    REQUIRE( numCols( matrix1 ) == numCols( matrix2 ) );
    REQUIRE( get_nnz( matrix1 ) == get_nnz( matrix2 ) );

    auto M = numRows( matrix1 );
    auto N = numCols( matrix1 );
    auto mat1_ptr = &matrix1;
    auto mat2_ptr = &matrix2;

    std::size_t num_matches = 0;

    Kokkos::parallel_reduce(
        "psi::test_crsmatrix::elem-wise-compare",
        Kokkos::MDRangePolicy< host_space, rank_type >( { 0, 0 }, { M, N } ),
        KOKKOS_LAMBDA ( const uint64_t i, const uint64_t j, std::size_t& l_nm ) {
          if ( access( *mat1_ptr, i, j ) == access( *mat2_ptr, i, j ) ) l_nm += 1;
        }, num_matches );

    REQUIRE( num_matches == M*N );
  }

  template< typename TCRSMatrix1, typename TCRSMatrix2 >
  inline void
  is_identical_crs( TCRSMatrix1& matrix1, TCRSMatrix2& matrix2 )
  {
    using host_space = Kokkos::DefaultHostExecutionSpace;

    REQUIRE( numRows( matrix1 ) != 0 );
    REQUIRE( numCols( matrix1 ) != 0 );
    REQUIRE( numRows( matrix1 ) == numRows( matrix2 ) );
    REQUIRE( numCols( matrix1 ) == numCols( matrix2 ) );
    REQUIRE( get_nnz( matrix1 ) == get_nnz( matrix2 ) );

    std::size_t nrows = numRows( matrix1 );
    std::size_t entries_size = matrix1.rowMap( nrows );

    auto mat1_ptr = &matrix1;
    auto mat2_ptr = &matrix2;

    std::size_t num_matches = 0;

    Kokkos::parallel_reduce(
        "psi::test_crsmatrix::crs_compare_entries",
        Kokkos::RangePolicy< host_space >( 0, entries_size ),
        KOKKOS_LAMBDA ( const uint64_t i, std::size_t& l_nm ) {
          if ( mat1_ptr->entry( i ) == mat2_ptr->entry( i ) ) l_nm += 1;
        }, num_matches );
    REQUIRE( num_matches == entries_size );

    num_matches = 0;
    Kokkos::parallel_reduce(
        "psi::test_crsmatrix::crs_compare_rowmap",
        Kokkos::RangePolicy< host_space >( 0, nrows + 1 ),
        KOKKOS_LAMBDA ( const uint64_t i, std::size_t& l_nm ) {
          if ( mat1_ptr->rowMap( i ) == mat2_ptr->rowMap( i ) ) l_nm += 1;
        }, num_matches );
    REQUIRE( num_matches == nrows + 1 );
  }

  template< typename TCRSMatrix >
  inline void
  print_matrix( TCRSMatrix& mat, std::string title="" )
  {
    if ( title.empty() ) title = "matrix:";
    auto nrows = mat.numRows();
    std::cout << title << std::endl;
    std::size_t j = 0;
    for ( std::size_t i = 0; i < nrows; ++i ) {
      auto end = mat.rowMap( i + 1 );
      for ( ; j < end; ++j ) {
        std::cout << mat.entry( j ) << " ";
      }
      std::cout << std::endl;
    }
  }
}

TEMPLATE_SCENARIO( "Generic functionality of Boolean CRSMatrices", "[crsmatrix][bool]",
                   crs_matrix::Dynamic,
                   crs_matrix::Compressed,
                   crs_matrix::Buffered,
                   crs_matrix::FullyBuffered,
                   crs_matrix::RangeDynamic,
                   crs_matrix::RangeBuffered,
                   crs_matrix::RangeFullyBuffered,
                   crs_matrix::RangeCompressed )
{
  typedef TestType spec_type;
  typedef CRSMatrix< spec_type, bool > crsmat_type;

  GIVEN( "A tiny matrix" )
  {
    constexpr const std::size_t nnz = 25;
    constexpr const std::size_t nrows = 10;
    constexpr const std::size_t ncols = 10;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::zero_matrix( simple );
    simple[ 0 ][ 3 ] = 1;
    simple[ 0 ][ 4 ] = 1;
    simple[ 0 ][ 5 ] = 1;
    simple[ 0 ][ 6 ] = 1;
    simple[ 1 ][ 7 ] = 1;
    simple[ 1 ][ 8 ] = 1;
    simple[ 1 ][ 9 ] = 1;
    simple[ 2 ][ 0 ] = 1;
    simple[ 2 ][ 1 ] = 1;
    simple[ 2 ][ 6 ] = 1;
    simple[ 4 ][ 1 ] = 1;
    simple[ 4 ][ 2 ] = 1;
    simple[ 5 ][ 0 ] = 1;
    simple[ 5 ][ 1 ] = 1;
    simple[ 5 ][ 2 ] = 1;
    simple[ 5 ][ 3 ] = 1;
    simple[ 5 ][ 4 ] = 1;
    simple[ 5 ][ 5 ] = 1;
    simple[ 5 ][ 6 ] = 1;
    simple[ 5 ][ 7 ] = 1;
    simple[ 5 ][ 8 ] = 1;
    simple[ 5 ][ 9 ] = 1;
    simple[ 8 ][ 9 ] = 1;
    simple[ 9 ][ 8 ] = 1;
    simple[ 9 ][ 9 ] = 1;

    REQUIRE( test_util::get_nnz( simple ) == nnz );

    WHEN( "A CRSMatrix instance is constructed by an external CRS matrix" )
    {
      crsmat_type matrix( test_util::to_external_crs( simple, nnz ) );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }
    }
  }

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2400;
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix( simple, nnz );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    REQUIRE( test_util::get_nnz( simple ) == nnz );

    WHEN( "A CRSMatrix instance is constructed by an external CRS matrix" )
    {
      crsmat_type matrix( test_util::to_external_crs( simple, nnz ) );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }

      AND_WHEN( "Another matrix is constructed using move constructor from the matrix" )
      {
        crsmat_type matrix2( std::move( matrix ) );

        THEN( "It should be identical to the original matrix" )
        {
          test_util::is_identical( matrix2, simple );
        }
      }
    }
  }

  GIVEN( "Two external matrices as blocks" )
  {
    constexpr const std::size_t nnz1 = 2400;
    constexpr const std::size_t nrows1 = 200;
    constexpr const std::size_t ncols1 = 200;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t nnz2 = 4000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1 + ncols2 >, nrows1 + nrows2 > appended;

    test_util::random_matrix_ranged( block1, nnz1 );
    test_util::random_matrix_ranged( block2, nnz2 );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );
    test_util::zero_matrix( appended );
    unsigned int i = 0;
    unsigned int j = 0;
    test_util::fill_block( appended, block1, i, j );
    test_util::fill_block( appended, block2, i, j );

    WHEN( "A CRSMatrix is constructed by two external CRS blocks" )
    {
      auto provider =
          [&block1, &block2, &nrows1, &ncols1, &nnz1, &nnz2]
          ( auto callback ) {
            callback( test_util::to_external_crs( block1, nnz1 ), 0, 0 );
            callback( test_util::to_external_crs( block2, nnz2 ), nrows1, ncols1 );
          };

      crsmat_type matrix( nrows1 + nrows2, ncols1 + ncols2, provider );

      THEN( "It should be identical to the original matrix" )
      {
        REQUIRE( matrix.nnz() == nnz1 + nnz2 );
        test_util::is_identical( matrix, appended );
      }

      AND_WHEN( "Another matrix is constructed using move constructor from the matrix" )
      {
        crsmat_type matrix2( std::move( matrix ) );

        THEN( "It should be identical to the original matrix" )
        {
          test_util::is_identical( matrix2, appended );
        }
      }
    }
  }

  GIVEN( "Two external matrices as non-consecutive blocks" )
  {
    constexpr const std::size_t nnz1 = 2400;
    constexpr const std::size_t nrows1 = 200;
    constexpr const std::size_t ncols1 = 200;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t znrows = 3;
    constexpr const std::size_t zncols = 6;

    constexpr const std::size_t nnz2 = 4000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1+zncols+ncols2 >, nrows1+znrows+nrows2 > appended;

    test_util::random_matrix_ranged( block1, nnz1 );
    test_util::random_matrix_ranged( block2, nnz2 );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );
    test_util::zero_matrix( appended );
    unsigned int i = 0;
    unsigned int j = 0;
    test_util::fill_block( appended, block1, i, j );
    i += znrows;
    j += zncols;
    test_util::fill_block( appended, block2, i, j );

    WHEN( "A CRSMatrix is constructed by two external CRS blocks" )
    {
      auto provider =
          [&]( auto callback ) {
            callback( test_util::to_external_crs( block1, nnz1 ), 0, 0 );
            callback( test_util::to_external_crs( block2, nnz2 ), nrows1 + znrows, ncols1 + zncols );
          };

      crsmat_type matrix( nrows1 + znrows + nrows2, ncols1 + zncols + ncols2, provider );

      THEN( "It should be identical to the original matrix" )
      {
        REQUIRE( matrix.nnz() == nnz1 + nnz2 );
        test_util::is_identical( matrix, appended );
      }

      AND_WHEN( "Another matrix is constructed using move constructor from the matrix" )
      {
        crsmat_type matrix2( std::move( matrix ) );

        THEN( "It should be identical to the original matrix" )
        {
          test_util::is_identical( matrix2, appended );
        }
      }
    }
  }
}

TEMPLATE_SCENARIO( "Specialised functionalities of non-Buffered Boolean CRSMatrix",
                   "[crsmatrix][bool]",
                   crs_matrix::Dynamic,
                   crs_matrix::Compressed,
                   crs_matrix::RangeDynamic,
                   crs_matrix::RangeCompressed )
{
  typedef TestType spec_type;
  typedef CRSMatrix< spec_type, bool > crsmat_type;

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2400;
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix_ranged( simple, nnz );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    REQUIRE( test_util::get_nnz( simple ) == nnz );

    WHEN( "A CRSMatrix instance is constructed by an external CRS matrix" )
    {
      crsmat_type matrix( test_util::to_external_crs( simple, nnz ) );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }

      AND_WHEN( "It is written to a file" )
      {
        std::string tmpfpath = get_tmpfile();
        save( matrix, tmpfpath );

        THEN( "It should be identical to the original matrix when loaded" )
        {
          crsmat_type matrix2;
          open( matrix2, tmpfpath );
          test_util::is_identical( matrix2, simple );
        }
      }

      AND_WHEN( "Another matrix is constructed using copy constructor from the matrix" )
      {
        crsmat_type matrix2( matrix );

        THEN( "It should be identical to the original matrix" )
        {
          test_util::is_identical( matrix2, simple );
        }
      }
    }
  }

  GIVEN( "Two external matrices as blocks" )
  {
    constexpr const std::size_t nnz1 = 2400;
    constexpr const std::size_t nrows1 = 200;
    constexpr const std::size_t ncols1 = 200;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t nnz2 = 4000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1 + ncols2 >, nrows1 + nrows2 > appended;

    test_util::random_matrix_ranged( block1, nnz1 );
    test_util::random_matrix_ranged( block2, nnz2 );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );
    test_util::zero_matrix( appended );
    unsigned int i = 0;
    unsigned int j = 0;
    test_util::fill_block( appended, block1, i, j );
    test_util::fill_block( appended, block2, i, j );

    WHEN( "A CRSMatrix is constructed by two external CRS blocks" )
    {
      auto provider =
          [&block1, &block2, &nrows1, &ncols1, &nnz1, &nnz2]
          ( auto callback ) {
            callback( test_util::to_external_crs( block1, nnz1 ), 0, 0 );
            callback( test_util::to_external_crs( block2, nnz2 ), nrows1, ncols1 );
          };

      crsmat_type matrix( nrows1 + nrows2, ncols1 + ncols2, provider );

      THEN( "It should be identical to the original matrix" )
      {
        REQUIRE( matrix.nnz() == nnz1 + nnz2 );
        test_util::is_identical( matrix, appended );
      }

      AND_WHEN( "It is written to a file" )
      {
        std::string tmpfpath = get_tmpfile();
        save( matrix, tmpfpath );

        THEN( "It should be identical to the original matrix when loaded" )
        {
          crsmat_type matrix2;
          open( matrix2, tmpfpath );
          test_util::is_identical( matrix2, appended );
        }
      }

      AND_WHEN( "Another matrix is constructed using copy constructor from the matrix" )
      {
        crsmat_type matrix2( matrix );

        THEN( "It should be identical to the original matrix" )
        {
          test_util::is_identical( matrix2, appended );
        }
      }
    }
  }

  GIVEN( "Two external matrices as non-consecutive blocks" )
  {
    constexpr const std::size_t nnz1 = 2400;
    constexpr const std::size_t nrows1 = 200;
    constexpr const std::size_t ncols1 = 200;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t znrows = 3;
    constexpr const std::size_t zncols = 6;

    constexpr const std::size_t nnz2 = 4000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1+zncols+ncols2 >, nrows1+znrows+nrows2 > appended;

    test_util::random_matrix_ranged( block1, nnz1 );
    test_util::random_matrix_ranged( block2, nnz2 );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );
    test_util::zero_matrix( appended );
    unsigned int i = 0;
    unsigned int j = 0;
    test_util::fill_block( appended, block1, i, j );
    i += znrows;
    j += zncols;
    test_util::fill_block( appended, block2, i, j );

    WHEN( "A CRSMatrix is constructed by two external CRS blocks" )
    {
      auto provider =
          [&]( auto callback ) {
            callback( test_util::to_external_crs( block1, nnz1 ), 0, 0 );
            callback( test_util::to_external_crs( block2, nnz2 ), nrows1 + znrows, ncols1 + zncols );
          };

      crsmat_type matrix( nrows1 + znrows + nrows2, ncols1 + zncols + ncols2, provider );

      THEN( "It should be identical to the original matrix" )
      {
        REQUIRE( matrix.nnz() == nnz1 + nnz2 );
        test_util::is_identical( matrix, appended );
      }

      AND_WHEN( "It is written to a file" )
      {
        std::string tmpfpath = get_tmpfile();
        save( matrix, tmpfpath );

        THEN( "It should be identical to the original matrix when loaded" )
        {
          crsmat_type matrix2;
          open( matrix2, tmpfpath );
          test_util::is_identical( matrix2, appended );
        }
      }

      AND_WHEN( "Another matrix is constructed using copy constructor from the matrix" )
      {
        crsmat_type matrix2( matrix );

        THEN( "It should be identical to the original matrix" )
        {
          test_util::is_identical( matrix2, appended );
        }
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG( "Specialised functionalities of Compressed Boolean CRS Matrix", "[crsmatrix][bool]",
                       ( ( typename T, typename U, typename V, int W /* dummy type unless compile error */ ), T, U, V, W ),
                       ( crs_matrix::Compressed, crs_matrix::Buffered, crs_matrix::FullyBuffered, 0 ),
                       ( crs_matrix::RangeCompressed, crs_matrix::RangeBuffered, crs_matrix::RangeFullyBuffered, 0 ) )
{
  typedef T crsmat_spec_type;
  typedef U crsmat_buffered_spec_type;
  typedef V crsmat_fully_buffered_spec_type;
  typedef CRSMatrix< crsmat_spec_type, bool > crsmat_type;

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2000;
    constexpr const std::size_t nrows = 1000;
    constexpr const std::size_t ncols = 500;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix( simple, nnz );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    REQUIRE( test_util::get_nnz( simple ) == nnz );

    WHEN( "The CRSMatrix is constructed by a Buffered CRS matrix")
    {
      crsmat_type matrix( test_util::to_external_crs( simple, nnz ),
                          crsmat_buffered_spec_type() );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }
    }

    WHEN( "The CRSMatrix is assigned by a Buffered CRS matrix")
    {
      crsmat_type matrix;
      make_buffered_t< crsmat_type > bmatrix( test_util::to_external_crs( simple, nnz ) );
      matrix.assign( bmatrix );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }
    }

    WHEN( "The CRSMatrix is constructed by a FullyBuffered CRS matrix")
    {
      crsmat_type matrix( test_util::to_external_crs( simple, nnz ),
                          crsmat_fully_buffered_spec_type() );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }
    }

    WHEN( "The CRSMatrix is assigned by a FullyBuffered CRS matrix")
    {
      crsmat_type matrix;
      make_fully_buffered_t< crsmat_type > bmatrix( test_util::to_external_crs( simple, nnz ) );
      matrix.assign( bmatrix );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG( "Specialised functionalities of Range Boolean CRS Matrix", "[crsmatrix][bool]",
                   ( ( typename T, typename U, int V /* dummy type unless compile error */ ), T, U, V ),
                   ( crs_matrix::RangeDynamic, crs_matrix::Dynamic, 0 ),
                   ( crs_matrix::RangeDynamic, crs_matrix::Compressed, 0 ),
                   ( crs_matrix::RangeDynamic, crs_matrix::Buffered, 0 ),
                   ( crs_matrix::RangeDynamic, crs_matrix::FullyBuffered, 0 ),
                   ( crs_matrix::RangeBuffered, crs_matrix::Dynamic, 0 ),
                   ( crs_matrix::RangeBuffered, crs_matrix::Compressed, 0 ),
                   ( crs_matrix::RangeBuffered, crs_matrix::Buffered, 0 ),
                   ( crs_matrix::RangeBuffered, crs_matrix::FullyBuffered, 0 ),
                   ( crs_matrix::RangeFullyBuffered, crs_matrix::Dynamic, 0 ),
                   ( crs_matrix::RangeFullyBuffered, crs_matrix::Compressed, 0 ),
                   ( crs_matrix::RangeFullyBuffered, crs_matrix::Buffered, 0 ),
                   ( crs_matrix::RangeFullyBuffered, crs_matrix::FullyBuffered, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::Dynamic, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::Compressed, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::Buffered, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::FullyBuffered, 0 ) )
{
  typedef CRSMatrix< T, bool > crsmat_range_type;
  typedef CRSMatrix< U, bool > crsmat_basic_type;

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2400;
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix_ranged( simple, nnz );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    REQUIRE( test_util::get_nnz( simple ) == nnz );

    WHEN( "The Range CRSMatrix is assigned by a Basic CRSMatrix")
    {
      crsmat_basic_type matrix( test_util::to_external_crs( simple, nnz ) );
      crsmat_range_type r_matrix;
      r_matrix.assign( matrix );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( r_matrix, simple );
      }
    }

    WHEN( "The Basic CRSMatrix is assigned by a Range CRSMatrix")
    {
      crsmat_range_type r_matrix( test_util::to_external_crs( simple, nnz ) );
      crsmat_basic_type matrix;
      matrix.assign( r_matrix );

      THEN( "It should be identical to the original matrix" )
      {
        test_util::is_identical( matrix, simple );
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG( "Merging two distance indices", "[crsmatrix][bool]",
                   ( ( typename T, typename U, int V /* dummy type unless compile error */ ), T, U, V ),
                   ( crs_matrix::Compressed, crs_matrix::Dynamic, 0 ),
                   ( crs_matrix::Compressed, crs_matrix::Buffered, 0 ),
                   ( crs_matrix::Compressed, crs_matrix::FullyBuffered, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::RangeDynamic, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::RangeBuffered, 0 ),
                   ( crs_matrix::RangeCompressed, crs_matrix::RangeFullyBuffered, 0 ) )
{
  typedef CRSMatrix< T, bool > crsmat_type;
  typedef CRSMatrix< U, bool > crsmat_mutable_type;

  GIVEN( "Two zero CRS matrices" )
  {
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    std::array< std::array< bool, ncols >, nrows > matrix1;
    std::array< std::array< bool, ncols >, nrows > matrix2;
    std::array< std::array< bool, ncols >, nrows > zero;
    test_util::zero_matrix( matrix1 );
    test_util::zero_matrix( matrix2 );
    test_util::zero_matrix( zero );

    WHEN( "They got merged" )
    {
      crsmat_type mat1( test_util::to_external_crs( matrix1 ) );
      crsmat_type mat2( test_util::to_external_crs( matrix2 ) );
      crsmat_type mmat;
      mmat.assign( merge_distance_index< crsmat_mutable_type >( mat1, mat2 ) );

      THEN( "The resulting merged matrix should be a zero matrix" )
      {
        REQUIRE( mat1.nnz() == 0 );
        REQUIRE( mat2.nnz() == 0 );
        REQUIRE( mmat.nnz() == 0 );
        test_util::is_identical( mmat, zero );
      }
    }
  }

  GIVEN( "Two CRS matrices, one is a random matrix and the other is a zero one" )
  {
    constexpr const std::size_t nnz = 4200;
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    std::array< std::array< bool, ncols >, nrows > matrix;
    std::array< std::array< bool, ncols >, nrows > zero;

    test_util::zero_matrix( zero );
    test_util::random_matrix_ranged( matrix, nnz );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    crsmat_type matn( test_util::to_external_crs( matrix ) );
    crsmat_type matz( test_util::to_external_crs( zero ) );
    crsmat_type mmat;

    WHEN( "The random matrix is merged with the zero one" )
    {
      mmat.assign( merge_distance_index< crsmat_mutable_type >( matn, matz ) );

      THEN( "The merged matrix should be identical to the non-zero one" )
      {
        REQUIRE( matz.nnz() == 0 );
        REQUIRE( mmat.nnz() == 4200 );
        test_util::is_identical( mmat, matn );
      }
    }

    WHEN( "They the zero matrix is merged with the random one" )
    {
      mmat.assign( merge_distance_index< crsmat_mutable_type >( matz, matn ) );

      THEN( "The merged matrix should be identical to the non-zero one" )
      {
        REQUIRE( matz.nnz() == 0 );
        REQUIRE( mmat.nnz() == 4200 );
        test_util::is_identical( mmat, matn );
      }
    }
  }

  GIVEN( "Two CRS matrices constructed from two random external matrices" )
  {
    constexpr const std::size_t nnz1 = 5200;
    constexpr const std::size_t nnz2 = 4000;
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    std::array< std::array< bool, ncols >, nrows > matrix1;
    std::array< std::array< bool, ncols >, nrows > matrix2;
    std::array< std::array< bool, ncols >, nrows > merged;

    test_util::random_matrix_ranged( matrix1, nnz1 );
    test_util::random_matrix_ranged( matrix2, nnz2 );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    test_util::zero_matrix( merged );
    for ( std::size_t i = 0; i < nrows; ++i ) {
      for ( std::size_t j = 0; j < ncols; ++j ) {
        if ( matrix1[i][j] == true || matrix2[i][j] == true ) merged[i][j] = true;
      }
    }

    WHEN( "They got merged" )
    {
      crsmat_type mat1( test_util::to_external_crs( matrix1 ) );
      crsmat_type mat2( test_util::to_external_crs( matrix2 ) );
      crsmat_type mmat;
      mmat.assign( merge_distance_index< crsmat_mutable_type >( mat1, mat2 ) );

      THEN( "It should be true whereever either of matrices are true" )
      {
        test_util::is_identical( mmat, merged );
      }
    }
  }

  GIVEN( "Two triangular CRS matrices whose merging yeilds all-ones matrix" )
  {
    constexpr const std::size_t nrows = 200;
    constexpr const std::size_t ncols = 200;
    constexpr const std::size_t offset = 20;
    std::array< std::array< bool, ncols >, nrows > matrix1;
    std::array< std::array< bool, ncols >, nrows > matrix2;
    std::array< std::array< bool, ncols >, nrows > allones;

    test_util::zero_matrix( matrix1 );
    test_util::zero_matrix( matrix2 );
    test_util::zero_matrix( allones );
    for ( std::size_t i = 0; i < nrows ; ++i ) {
      for ( std::size_t j = 0; j < ncols ; ++j ) {
        if ( i + offset < j ) matrix1[i][j] = true;
        else matrix2[i][j] = true;
        allones[i][j] = true;
      }
    }

    WHEN( "They got merged" )
    {
      crsmat_type mat1( test_util::to_external_crs( matrix1 ) );
      crsmat_type mat2( test_util::to_external_crs( matrix2 ) );
      crsmat_type mmat;
      mmat.assign( merge_distance_index< crsmat_mutable_type >( mat1, mat2 ) );

      THEN( "It should yeild an all-ones matrix" )
      {
        test_util::is_identical( mmat, allones );
      }
    }
  }
}

TEMPLATE_SCENARIO_SIG( "Merging two distance index with large dimensions", "[crsmatrix][bool]",
                       ( ( typename T, typename U, int V /* dummy type unless compile error */ ), T, U, V ),
                       ( crs_matrix::Compressed, crs_matrix::Dynamic, 0 ),
                       ( crs_matrix::Compressed, crs_matrix::Buffered, 0 ),
                       ( crs_matrix::Compressed, crs_matrix::FullyBuffered, 0 ),
                       ( crs_matrix::RangeCompressed, crs_matrix::RangeDynamic, 0 ),
                       ( crs_matrix::RangeCompressed, crs_matrix::RangeBuffered, 0 ),
                       ( crs_matrix::RangeCompressed, crs_matrix::RangeFullyBuffered, 0 ) )
{
  typedef CRSMatrix< T, bool, uint16_t, uint32_t > crsmat_type;
  typedef CRSMatrix< U, bool, uint16_t, uint32_t > crsmat_mutable_type;

  GIVEN( "Two CRS matrices with very large dimensions" )
  {
    typedef typename crsmat_type::ordinal_type ordinal_type;

    constexpr const std::size_t nnz1 = 2400;
    constexpr const std::size_t nrows1 = 200;
    constexpr const std::size_t ncols1 = 200;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t nnz2 = 4000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    test_util::random_matrix_ranged( block1, nnz1 );
    test_util::random_matrix_ranged( block2, nnz2 );
    INFO( "Seed for the random number generator: " << rnd::get_iseed() );

    constexpr const std::size_t nrows = std::numeric_limits< ordinal_type >::max() / 2;
    constexpr const std::size_t ncols = nrows;

    auto blk1 = test_util::to_external_crs( block1 );
    auto blk2 = test_util::to_external_crs( block2 );

    auto list = []( auto ...args ) {
      return [=]( auto fn ) { fn( args... ); };
    };

/*
    auto feed_blocks = []( auto xs ) {
      return [=]( auto callback ) {
        auto impl = [=]( auto fn ) {
          return [=]( auto* blk, std::size_t sr, std::size_t sc, auto ...rest ) {  // main
            callback( *blk, sr, sc );
            if constexpr ( sizeof...( rest ) > 0 ) {
              fn( fn )( rest... );
            }
          };
        };
        xs( impl( impl ) );  // apply `main` lambda on all parameters
      };
    };
*/

    auto feed_blocks = []( auto xs ) {
      return [=]( auto callback ) {
        auto impl = [=]( auto* blk, std::size_t sr, std::size_t sc ) -> void {
          callback( *blk, sr, sc );
        };
        xs( impl );
      };
    };

    auto feed_blocks_2 = []( auto xs ) {
      return [=]( auto callback ) {
        auto impl = [=]( auto* blk1, std::size_t sr1, std::size_t sc1,
                         auto* blk2, std::size_t sr2, std::size_t sc2  ) -> void {
          callback( *blk1, sr1, sc1 );
          callback( *blk2, sr2, sc2 );
        };
        xs( impl );
      };
    };

    std::size_t i = nrows - nrows2;
    std::size_t j = ncols - ncols2;
    crsmat_type mat1( nrows, ncols, feed_blocks( list( &blk1, 0, 0 ) ), nnz1 );
    crsmat_type mat2( nrows, ncols, feed_blocks( list( &blk2, i, j ) ), nnz2 );
    crsmat_type merged( nrows, ncols, feed_blocks_2( list( &blk1, 0, 0, &blk2, i, j ) ), nnz1 + nnz2 );

    WHEN( "They got merged" )
    {
      crsmat_type mmat;
      mmat.assign( merge_distance_index< crsmat_mutable_type >( mat1, mat2 ) );

      THEN( "It should be true whereever either of matrices are true" )
      {
        REQUIRE( mat1.nnz() == nnz1 );
        REQUIRE( mat2.nnz() == nnz2 );
        REQUIRE( mmat.nnz() == nnz1 + nnz2 );
        test_util::is_identical_crs( mmat, merged );
      }
    }
  }
}
