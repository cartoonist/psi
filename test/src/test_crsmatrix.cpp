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

TEMPLATE_SCENARIO( "Generic functionality of Boolean CRSMatrix", "[crsmatrix]",
                   crs_matrix::Dynamic, crs_matrix::Succinct )
{
  typedef TestType spec_type;
  typedef CRSMatrix< spec_type > crsmat_type;
  typedef pairg::matrixOps crs_traits_type;  // external CRS traits type

  Kokkos::initialize();

  auto populate_random =
      []( auto& matrix, std::size_t nnz ) {
        auto nrows = matrix.size();
        auto ncols = matrix[0].size();

        for ( std::size_t i = 0; i < nrows; ++i ) {
          for ( std::size_t j = 0; j < ncols; ++j ) matrix[i][j] = 0;
        }

        for ( std::size_t i = 0; i < nnz; ) {
          auto r = random::random_index( nrows );
          auto c = random::random_index( ncols );
          if ( matrix[r][c] == 0 ) {
            matrix[r][c] = 1;
            ++i;
          }
        }
      };

  auto convert_to_crs =
      []( auto const& matrix, std::size_t nnz ) {
        typedef typename crs_traits_type::size_type size_type;
        typedef typename crs_traits_type::lno_nnz_view_t lno_nnz_view_t;
        typedef typename crs_traits_type::scalar_view_t scalar_view_t;
        typedef typename crs_traits_type::lno_view_t lno_view_t;
        typedef typename crs_traits_type::crsMat_t crsMat_t;

        auto nrows = matrix.size();
        auto ncols = matrix[0].size();

        lno_nnz_view_t entries( "entries", nnz );
        scalar_view_t values( "values", nnz );
        lno_view_t rowmap( "rowmap", nrows + 1 );

        for ( size_type i = 0; i < nnz; i++ ) values( i ) = 1;  // boolean

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
      };

  auto append_components =
      []( auto const& comp1, auto const& comp2 ) {
        auto nrows = comp1.size() + comp2.size();
        auto ncols = comp1[0].size() + comp2[0].size();
        std::vector< std::vector< bool > > result;
        result.resize( nrows );
        for ( std::size_t i = 0; i < nrows; ++i ) {
          result[i].resize( ncols );
          for ( std::size_t j = 0; j < ncols; ++j ) result[i][j] = 0;
        }
        for ( std::size_t i = 0; i < comp1.size(); ++i ) {
          for ( std::size_t j = 0; j < comp1[0].size(); ++j ) {
            result[i][j] = comp1[i][j];
          }
        }
        for ( std::size_t i = 0; i < comp2.size(); ++i ) {
          for ( std::size_t j = 0; j < comp2[0].size(); ++j ) {
            result[i + comp1.size()][j + comp1[0].size()] = comp2[i][j];
          }
        }

        return result;
      };

  auto is_identical =
      []( auto const& mat1, auto const& mat2 ) {
        REQUIRE( mat1.numRows() != 0 );
        REQUIRE( mat1.numCols() != 0 );
        REQUIRE( mat1.numRows() == mat2.size() );
        REQUIRE( mat1.numCols() == mat2[0].size() );
        for ( std::size_t i = 0; i < mat1.numRows(); ++i ) {
          for ( std::size_t j = 0; j < mat1.numCols(); ++j ) {
            REQUIRE( mat1( i, j ) == mat2[i][j] );
          }
        }
      };

  GIVEN( "A simple external matrix" )
  {
    constexpr const std::size_t nnz = 2000;
    constexpr const std::size_t nrows = 1000;
    constexpr const std::size_t ncols = 500;
    std::array< std::array< bool, ncols >, nrows > simple;

    populate_random( simple, nnz );

    WHEN( "Constructed by an external CRS matrix" )
    {
      crsmat_type matrix( convert_to_crs( simple, nnz ) );

      THEN( "It should be identical to the original matrix" )
      {
        is_identical( matrix, simple );
      }

      AND_WHEN( "It is written to a file" )
      {
        std::string tmpfpath = get_tmpfile();
        save( matrix, tmpfpath );

        THEN( "It should be identical to the original matrix when loaded" )
        {
          crsmat_type matrix2;
          open( matrix2, tmpfpath );
          is_identical( matrix2, simple );
        }
      }
    }
  }

  GIVEN( "A two external matrix as components" )
  {
    constexpr const std::size_t nnz1 = 2000;
    constexpr const std::size_t nrows1 = 1000;
    constexpr const std::size_t ncols1 = 500;
    std::array< std::array< bool, ncols1 >, nrows1 > comp1;

    constexpr const std::size_t nnz2 = 1000;
    constexpr const std::size_t nrows2 = 400;
    constexpr const std::size_t ncols2 = 400;
    std::array< std::array< bool, ncols2 >, nrows2 > comp2;

    populate_random( comp1, nnz1 );
    populate_random( comp2, nnz2 );

    WHEN( "Constructed by two external CRS components" )
    {
      auto provider =
          [&comp1, &comp2, &nnz1, &nnz2, convert_to_crs]( auto callback ) {
            callback( convert_to_crs( comp1, nnz1 ), 0, 0, 0 );
            callback( convert_to_crs( comp2, nnz2 ), 1000, 500, nnz1 );
          };
      crsmat_type matrix( 1400, 900, 3000, provider );

      THEN( "It should be identical to the original matrix" )
      {
        is_identical( matrix, append_components( comp1, comp2 ) );
      }

      AND_WHEN( "It is written to a file" )
      {
        std::string tmpfpath = get_tmpfile();
        save( matrix, tmpfpath );

        THEN( "It should be identical to the original matrix when loaded" )
        {
          crsmat_type matrix2;
          open( matrix2, tmpfpath );
          is_identical( matrix2, append_components( comp1, comp2 ) );
        }
      }
    }
  }

  Kokkos::finalize();
}
