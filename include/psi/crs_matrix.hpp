/**
 *    @file  crs_matrix.hpp
 *   @brief  CRS matrix class definition.
 *
 *  This header file defines CRS matrix class.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Nov 10, 2020  23:32
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2020, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_CRS_MATRIX_HPP__
#define PSI_CRS_MATRIX_HPP__

#include <cinttypes>
#include <vector>

#include <sdsl/bit_vectors.hpp>

#include "utils.hpp"

namespace psi {
  /* specialisation tags */
  namespace crs_matrix {
    struct Succinct { };
    struct Dynamic { };
  }  /* --- end of namespace crs_matrix --- */

  template< typename TSpec = crs_matrix::Succinct,
            typename TValue = bool,
            typename TOrdinal = uint32_t,
            typename TSize = uint64_t >
  class CRSMatrix;

  template< typename ...TArgs >
  using DynamicCRSMatrix = CRSMatrix< crs_matrix::Dynamic, TArgs... >;

  template< typename ...TArgs >
  using SuccinctCRSMatrix = CRSMatrix< crs_matrix::Succinct, TArgs... >;

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::Dynamic, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef crs_matrix::Dynamic spec_type;
    typedef bool value_type;
    typedef TOrdinal ordinal_type;
    typedef std::vector< ordinal_type > entries_type;
    typedef std::vector< size_type > rowmap_type;
    typedef entries_type mutable_entries_type;
    typedef rowmap_type mutable_rowmap_type;
    typedef TSize size_type;
    /* === LIFECYCLE === */
    CRSMatrix( ) : num_cols( 0 ) { }

    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix&& ext )
      : entries( ext.nnz() ), rowmap( ext.numRows() + 1 ), num_cols( ext.numCols() )
    {
      this->fill_partial( ext );
    }

    /**
     *  NOTE: The `nnz` parameter is considered as an estimation of the number of
     *        non-zero values. The actual size will be determined by the provided
     *        components.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, size_type nnz, TCallback callback )
      : rowmap( nrows + 1 ), num_cols( ncols )
    {
      this->entries.reserve( nnz );
      auto partial_ctor =
          [this]( auto&& mat, ordinal_type srow, ordinal_type scol, size_type snnz ) {
            size_type new_nnz = snnz + mat.nnz();
            this->entries.resize( new_nnz );
            this->fill_partial( mat, srow, scol, snnz );
          };
      callback( partial_ctor );
    }
    /* === OPERATORS === */
    inline value_type
    operator()( ordinal_type i, ordinal_type j ) const
    {
      if ( i >= this->numRows() || j >= this->numCols() ) {
        throw std::runtime_error( "index out of range" );
      }

      return std::binary_search( this->entries.begin() + this->rowmap[ i ],
                                 this->entries.begin() + this->rowmap[ i + 1 ], j );
    }
    /* === METHODS === */
    inline ordinal_type
    numCols( ) const
    {
      return this->num_cols;
    }

    inline ordinal_type
    numRows( ) const
    {
      return this->rowmap.size() - 1;
    }

    inline size_type
    nnz( ) const
    {
      return this->entries.size();
    }

    inline void
    clear( )
    {
      psi::clear( this->entries );
      psi::clear( this->rowmap );
      this->num_cols = 0;
    }

    inline void
    serialize( std::ostream& out ) const
    {
      psi::serialize( out, this->entries );
      psi::serialize( out, this->rowmap );
      psi::serialize( out, static_cast< uint64_t >( this->num_cols ) );
    }

    inline void
    load( std::istream& in )
    {
      this->clear();
      uint64_t ncols;
      psi::deserialize( in, this->entries );
      psi::deserialize( in, this->rowmap );
      psi::deserialize( in, ncols );
      this->num_cols = ncols;
    }
  private:
    /* === DATA MEMBERS === */
    entries_type entries;
    rowmap_type rowmap;
    ordinal_type num_cols;
    /* === METHODS === */
    /**
     *  @brief  Populate a component in the matrix by an external matrix.
     *
     *  @param[in]  ex The external matrix representing the component values.
     *  @param[in]  srow The starting row of the component to be filled.
     *  @param[in]  scol The starting column of the component to be filled.
     *
     *  We call a matrix A is a component matrix when it is composed by smaller
     *  matrices, called components, such that all components are mutually disjoint.
     *  This means all non-zero values in each row or column in the matrix only belong
     *  to one of the components; i.e. no two rows or columns are shared between two or
     *  more components.
     *
     *  When a matrix is a component matrix, it can be gradually populated by components
     *  in parallel. For example, the matrix
     *
     *        | E1   0   0   0   0 |
     *        |  0  E2   0   0   0 |
     *    A = |  0   0  E3   0   0 |
     *        |  0   0   0  E4   0 |
     *        |  0   0   0   0  E5 |
     *
     *  can be populated by filling $E_i$ matrices separately where $i \in [1..5]$. This
     *  method fill a component represented by `ex` matrix starting at [srow][scol].
     *
     *  NOTE: The intermediate matrices are not valid until it completely populated.
     *
     *  NOTE: There should not be any gap between $E_i$s; i.e. the first element of
     *        $E_{i+1}$ should be diagonally adjacent of the last element of $E_i$.
     *
     *  NOTE: This method assumes that the required memory is already allocated.
     */
    template< typename TCrsMatrix >
    inline void
    fill_partial( TCrsMatrix const& ex, ordinal_type srow=0, ordinal_type scol=0, size_type snnz=0 )
    {
      this->rowmap[0] = 0;
      _fill_entries_partial( ex, scol, snnz );
      _fill_rowmap_partial( ex, srow, snnz );
    }

    template< typename TCrsMatrix >
    inline void
    _fill_entries_partial( TCrsMatrix const& ex, ordinal_type scol=0, size_type snnz=0 )
    {
      auto size = ex.graph.entries.extent( 0 );
      auto start = this->entries.begin() + snnz;
      std::copy( ex.graph.entries.data(), ex.graph.entries.data() + size, start );
      std::transform( start, start + size, start,
                      [scol]( ordinal_type e ) { return e + scol; } );
    }

    template< typename TCrsMatrix >
    inline void
    _fill_rowmap_partial( TCrsMatrix const& ex, ordinal_type srow=0, size_type snnz=0 )
    {
      auto size = ex.graph.row_map.extent( 0 );
      auto start = this->rowmap.begin() + srow + 1;
      std::copy( ex.graph.row_map.data() + 1, ex.graph.row_map.data() + size, start );
      std::transform( start, start + size - 1, start,
                      [snnz]( size_type e ) { return e + snnz; } );
    }
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::Succinct, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef crs_matrix::Succinct spec_type;
    typedef bool value_type;
    typedef TOrdinal ordinal_type;
    typedef TSize size_type;
    typedef sdsl::enc_vector< sdsl::coder::elias_delta > entries_type;
    typedef sdsl::bit_vector rowmap_type;
    typedef entries_type::int_vector_type mutable_entries_type;
    typedef rowmap_type mutable_rowmap_type;
    typedef rowmap_type::select_1_type select_type;
    /* === LIFECYCLE === */
    CRSMatrix( ) : num_cols( 0 ) { }

    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix&& ext )
      : num_cols( ext.numCols() )
    {
      mutable_entries_type mut_entries( ext.nnz() );
      mutable_rowmap_type mut_rowmap( ext.numRows() + ext.nnz() );
      this->fill_partial( mut_entries, mut_rowmap, ext );
      this->entries = std::move( mut_entries );
      this->rowmap = std::move( mut_rowmap );
      sdsl::util::init_support( this->select, &this->rowmap );
    }

    /**
     *  NOTE: The `nnz` parameter is considered as an estimation of the number of
     *        non-zero values. The actual size will be determined by the provided
     *        components.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, size_type nnz, TCallback callback )
      : num_cols( ncols )
    {
      mutable_entries_type mut_entries( nnz );
      mutable_rowmap_type mut_rowmap( nnz + nrows );
      auto partial_ctor =
          [this, &mut_entries, &mut_rowmap, nrows]
          ( auto&& mat, ordinal_type srow, ordinal_type scol, size_type snnz ) {
            size_type new_nnz = snnz + mat.nnz();
            mut_entries.resize( new_nnz );
            mut_rowmap.resize( new_nnz + nrows );
            this->fill_partial( mut_entries, mut_rowmap, mat, srow, scol, snnz );
          };
      callback( partial_ctor );
      this->entries = std::move( mut_entries );
      this->rowmap = std::move( mut_rowmap );
      sdsl::util::init_support( this->select, &this->rowmap );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      sdsl::util::init_support( this->select, &this->rowmap );
    }

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      sdsl::util::init_support( this->select, &this->rowmap );
    }
    /* === OPERATORS === */
    inline CRSMatrix&
    operator=( CRSMatrix const& other )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      sdsl::util::init_support( this->select, &this->rowmap );
      return *this;
    }

    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      sdsl::util::init_support( this->select, &this->rowmap );
      return *this;
    }

    inline bool
    operator()( ordinal_type i, ordinal_type j ) const
    {
      if ( i >= this->numRows() || j >= this->numCols() ) {
        throw std::runtime_error( "index out of range" );
      }

      return std::binary_search( this->entries.begin() + this->row_map( i ),
                                 this->entries.begin() + this->row_map( i + 1 ), j );
    }
    /* === METHODS === */
    inline ordinal_type
    numCols( ) const
    {
      return this->num_cols;
    }

    inline ordinal_type
    numRows( ) const
    {
      return this->rowmap.size() - this->entries.size();
    }

    inline size_type
    nnz( ) const
    {
      return this->entries.size();
    }

    inline void
    clear( )
    {
      psi::clear( this->entries );
      sdsl::util::clear( this->rowmap );
      sdsl::util::clear( this->select );
      this->num_cols = 0;
    }

    inline void
    serialize( std::ostream& out ) const
    {
      psi::serialize( out, this->entries );
      this->rowmap.serialize( out );
      psi::serialize( out, static_cast< uint64_t >( this->num_cols ) );
    }

    inline void
    load( std::istream& in )
    {
      this->clear();
      uint64_t ncols;
      psi::deserialize( in, this->entries );
      this->rowmap.load( in );
      psi::deserialize( in, ncols );
      this->num_cols = ncols;
      sdsl::util::init_support( this->select, &this->rowmap );
    }
  private:
    /* === DATA MEMBERS === */
    entries_type entries;
    rowmap_type rowmap;
    select_type select;
    ordinal_type num_cols;
    /* === METHODS === */
    inline size_type
    row_map( ordinal_type i ) const
    {
      if ( i == 0 ) return 0;
      return this->select( i ) + 1 - i /* subtract i number of 1s */;
    }
    /**
     *  @brief  Populate a component in the matrix by an external matrix.
     *
     *  @param[in]  mut_entries The mutable entries array.
     *  @param[in]  mut_rowmap The mutable rowmap bit vector.
     *  @param[in]  ex The external matrix representing the component values.
     *  @param[in]  srow The starting row of the component to be filled.
     *  @param[in]  scol The starting column of the component to be filled.
     *
     *  See the function documentation for `crs_matrix::Dynamic` counterpart.
     *
     *  NOTE: The intermediate matrices are not valid until it completely populated.
     *
     *  NOTE: This method assumes that the required memory is already allocated.
     */
    template< typename TCrsMatrix >
    inline void
    fill_partial( mutable_entries_type& mut_entries, mutable_rowmap_type& mut_rowmap,
                  TCrsMatrix const& ex, ordinal_type srow=0, ordinal_type scol=0, size_type snnz=0 )
    {
      _fill_entries_partial( mut_entries, ex, scol, snnz );
      _fill_rowmap_partial( mut_rowmap, ex, srow, snnz );
    }

    template< typename TCrsMatrix >
    inline void
    _fill_entries_partial( mutable_entries_type& mut_entries, TCrsMatrix const& ex,
                           ordinal_type scol=0, size_type snnz=0 )
    {
      auto size = ex.graph.entries.extent( 0 );
      auto start = mut_entries.begin() + snnz;
      std::copy( ex.graph.entries.data(), ex.graph.entries.data() + size, start );
      std::transform( start, start + size, start,
                      [scol]( ordinal_type e ) { return e + scol; } );
    }

    template< typename TCrsMatrix >
    inline void
    _fill_rowmap_partial( mutable_rowmap_type& mut_rowmap, TCrsMatrix const& ex,
                          ordinal_type srow=0, size_type snnz=0 )
    {
      auto size = ex.graph.row_map.extent( 0 );
      size_type start = snnz + srow;
      for ( ordinal_type i = 0; i < size-1; ++i ) {
        for ( ordinal_type j = 0; j < ex.graph.row_map( i+1 ) - ex.graph.row_map( i ); ++j ) {
          mut_rowmap[ start++ ] = 0;
        }
        mut_rowmap[ start++ ] = 1;
      }
    }
  };
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_CRS_MATRIX_HPP__ --- */
