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
  get_nnz( std::array< std::array< TValue, NCols >, NRows >& matrix ) {
    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );
    std::size_t nnz = 0;
    for ( std::size_t i = 0; i < nrows; ++i ) {
      for ( std::size_t j = 0; j < ncols; ++j ) {
        if ( access( matrix, i, j ) == 1 ) ++nnz;
      }
    }
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
  zero_matrix( TMatrix& matrix ) {
    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );
    for ( std::size_t i = 0; i < nrows; ++i ) {
      for ( std::size_t j = 0; j < ncols; ++j ) matrix[i][j] = 0;
    }
  }

  template< typename TMatrix >
  inline void
  random_matrix( TMatrix& matrix, std::size_t nnz ) {
    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );

    assert( nnz <= nrows * ncols );

    // initialise the matrix by zero
    zero_matrix( matrix );

    std::size_t i = 0;
    while ( i < nnz ) {
      auto r = random::random_index( nrows );
      auto c = random::random_index( ncols );
      if ( matrix[r][c] == 0 ) {
        matrix[r][c] = 1;
        ++i;
      }
    }
  }

  /* NOTE: The `matrix` cannot be const reference because of Buffered specialisations. */
  template< typename TMatrix, typename TExternalMatrixTraits = pairg::matrixOps >
  inline typename TExternalMatrixTraits::crsMat_t
  to_external_crs( TMatrix& matrix, std::size_t nnz=0 )
  {
    typedef typename TExternalMatrixTraits::size_type size_type;
    typedef typename TExternalMatrixTraits::lno_nnz_view_t lno_nnz_view_t;
    typedef typename TExternalMatrixTraits::scalar_view_t scalar_view_t;
    typedef typename TExternalMatrixTraits::lno_view_t lno_view_t;
    typedef typename TExternalMatrixTraits::crsMat_t crsMat_t;

    auto nrows = numRows( matrix );
    auto ncols = numCols( matrix );
    if ( nnz == 0 ) nnz = get_nnz( matrix );

    lno_nnz_view_t entries( "entries", nnz );
    scalar_view_t values( "values", nnz );
    lno_view_t rowmap( "rowmap", nrows + 1 );

    for ( size_type i = 0; i < nnz; ++i ) values( i ) = 1;  // boolean

    size_type eidx = 0;
    size_type ridx = 0;
    for ( size_type i = 0; i < nrows; ++i ) {
      rowmap[ ridx++ ] = eidx;
      for ( size_type j = 0; j < ncols; ++j ) {
        if ( matrix[i][j] == 1 ) entries[ eidx++ ] = j;
      }
    }
    rowmap[ ridx++ ] = eidx;
    assert( eidx == nnz );
    assert( ridx == nrows + 1 );

    return crsMat_t( "matrix", nrows, ncols, nnz, values, rowmap, entries );
  }

  template< typename TMatrix1, typename TMatrix2 >
  inline void
  fill_block( TMatrix1& matrix, TMatrix2 const& block,
              unsigned int &si, unsigned int &sj )
  {
    for ( std::size_t i = 0; i < block.size(); ++i ) {
      for ( std::size_t j = 0; j < block[0].size(); ++j ) {
        matrix[si + i][sj + j] = block[i][j];
      }
    }
    si += numRows( block );
    sj += numCols( block );
  }

  template< typename TMatrix1, typename TMatrix2 >
  inline void
  is_identical( TMatrix1& matrix1, TMatrix2& matrix2 ) {
    REQUIRE( numRows( matrix1 ) != 0 );
    REQUIRE( numCols( matrix1 ) != 0 );
    REQUIRE( numRows( matrix1 ) == numRows( matrix2 ) );
    REQUIRE( numCols( matrix1 ) == numCols( matrix2 ) );
    REQUIRE( get_nnz( matrix1 ) == get_nnz( matrix2 ) );
    for ( std::size_t i = 0; i < numRows( matrix1 ); ++i ) {
      for ( std::size_t j = 0; j < numCols( matrix1 ); ++j ) {
        REQUIRE( access( matrix1, i, j ) == access( matrix2, i, j ) );
      }
    }
  }
}

TEMPLATE_SCENARIO( "Generic functionality of Boolean CRSMatrices", "[crsmatrix][bool]",
                   crs_matrix::Dynamic,
                   crs_matrix::Compressed,
                   crs_matrix::Buffered,
                   crs_matrix::FullyBuffered )
{
  typedef TestType spec_type;
  typedef CRSMatrix< spec_type, bool > crsmat_type;

  Kokkos::initialize();

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2000;
    constexpr const std::size_t nrows = 1000;
    constexpr const std::size_t ncols = 500;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix( simple, nnz );

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
    constexpr const std::size_t nnz1 = 2000;
    constexpr const std::size_t nrows1 = 1000;
    constexpr const std::size_t ncols1 = 500;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t nnz2 = 1000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1 + ncols2 >, nrows1 + nrows2 > appended;

    test_util::random_matrix( block1, nnz1 );
    test_util::random_matrix( block2, nnz2 );
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
    constexpr const std::size_t nnz1 = 2000;
    constexpr const std::size_t nrows1 = 100;
    constexpr const std::size_t ncols1 = 500;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t znrows = 3;
    constexpr const std::size_t zncols = 6;

    constexpr const std::size_t nnz2 = 1500;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1+zncols+ncols2 >, nrows1+znrows+nrows2 > appended;

    test_util::random_matrix( block1, nnz1 );
    test_util::random_matrix( block2, nnz2 );
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
                   crs_matrix::Compressed )
{
  typedef TestType spec_type;
  typedef CRSMatrix< spec_type, bool > crsmat_type;

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2000;
    constexpr const std::size_t nrows = 1000;
    constexpr const std::size_t ncols = 500;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix( simple, nnz );

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
    constexpr const std::size_t nnz1 = 2000;
    constexpr const std::size_t nrows1 = 1000;
    constexpr const std::size_t ncols1 = 500;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t nnz2 = 1000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1 + ncols2 >, nrows1 + nrows2 > appended;

    test_util::random_matrix( block1, nnz1 );
    test_util::random_matrix( block2, nnz2 );
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
    constexpr const std::size_t nnz1 = 2000;
    constexpr const std::size_t nrows1 = 100;
    constexpr const std::size_t ncols1 = 500;
    std::array< std::array< bool, ncols1 >, nrows1 > block1;

    constexpr const std::size_t znrows = 3;
    constexpr const std::size_t zncols = 6;

    constexpr const std::size_t nnz2 = 1500;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > block2;

    std::array< std::array< bool, ncols1+zncols+ncols2 >, nrows1+znrows+nrows2 > appended;

    test_util::random_matrix( block1, nnz1 );
    test_util::random_matrix( block2, nnz2 );
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

SCENARIO( "Specialised functionalities of Compressed Boolean CRS Matrix", "[crsmatrix][bool]" )
{
  typedef CRSMatrix< crs_matrix::Compressed, bool > crsmat_type;

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2000;
    constexpr const std::size_t nrows = 1000;
    constexpr const std::size_t ncols = 500;
    std::array< std::array< bool, ncols >, nrows > simple;

    test_util::random_matrix( simple, nnz );

    REQUIRE( test_util::get_nnz( simple ) == nnz );

    WHEN( "The CRSMatrix is constructed by a Buffered CRS matrix")
    {
      crsmat_type matrix( test_util::to_external_crs( simple, nnz ),
                          crs_matrix::Buffered() );

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
                          crs_matrix::FullyBuffered() );

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

SCENARIO( "Finalise Kokkos" )
{
  Kokkos::finalize();
}
