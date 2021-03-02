/**
 *    @file  crs_matrix.hpp
 *   @brief  CRS matrix template class definition.
 *
 *  This header file defines CRS matrix template class.
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

#include <gum/basic_types.hpp>
#include <gum/iterators.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/enc_vector.hpp>

#include "utils.hpp"


namespace psi {
  /* specialisation tags */
  namespace crs_matrix {
    struct Dynamic { };
    struct Buffered { };
    struct FullyBuffered { };
    struct Compressed { };
  }  /* --- end of namespace crs_matrix --- */

  template< typename TSpec = crs_matrix::Compressed,
            typename TValue = bool,
            typename TOrdinal = uint32_t,
            typename TSize = uint64_t >
  class CRSMatrix;

  template< typename ...TArgs >
  using DynamicCRSMatrix = CRSMatrix< crs_matrix::Dynamic, TArgs... >;

  template< typename ...TArgs >
  using BufferedCRSMatrix = CRSMatrix< crs_matrix::Buffered, TArgs... >;

  template< typename ...TArgs >
  using CompressedCRSMatrix = CRSMatrix< crs_matrix::Compressed, TArgs... >;

  template< typename TNativeCRSMatrix >
  struct make_buffered;

  template< typename TSpec, typename ...TArgs >
  struct make_buffered< CRSMatrix< TSpec, TArgs... > > {
    typedef CRSMatrix< crs_matrix::Buffered, TArgs... > type;
  };

  template< typename TNativeCRSMatrix >
  using make_buffered_t = typename make_buffered< TNativeCRSMatrix >::type;

  template< typename TNativeCRSMatrix >
  struct make_fully_buffered;

  template< typename TSpec, typename ...TArgs >
  struct make_fully_buffered< CRSMatrix< TSpec, TArgs... > > {
    typedef CRSMatrix< crs_matrix::FullyBuffered, TArgs... > type;
  };

  template< typename TNativeCRSMatrix >
  using make_fully_buffered_t = typename make_fully_buffered< TNativeCRSMatrix >::type;

  template< typename TSpec, typename TNativeCRSMatrix >
  struct make_spec;

  template< typename TSpec1, typename TSpec2, typename ...TArgs >
  struct make_spec< TSpec1, CRSMatrix< TSpec2, TArgs... > > {
    typedef CRSMatrix< TSpec1, TArgs... > type;
  };

  template< typename TSpec, typename TNativeCRSMatrix >
  using make_spec_t = typename make_spec< TSpec, TNativeCRSMatrix >::type;

  template< typename TSpec, typename TValue, typename TOrdinal, typename TSize >
  struct CRSMatrixTraits;

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::Dynamic, bool, TOrdinal, TSize > {
    typedef std::vector< TOrdinal > entries_type;
    typedef std::vector< TSize > rowmap_type;
    typedef crs_matrix::Dynamic mutable_spec_type;

    template< typename TValue >
    static inline void
    init( std::vector< TValue >& )
    { /* NOOP */ }

    template< typename T >
    static inline void
    swap( T& a, T& b )
    {
      std::swap( a, b );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::Buffered, bool, TOrdinal, TSize > {
    typedef sdsl::int_vector_buffer< gum::widthof< TOrdinal >::value > entries_type;
    typedef std::vector< TSize > rowmap_type;
    typedef crs_matrix::Buffered mutable_spec_type;

    static inline void
    init( entries_type& entries )
    {
      entries = entries_type( get_tmpfile(), std::ios::out );
    }

    static inline void
    init( rowmap_type& )
    { /* NOOP */ }

    static inline void
    swap( entries_type& entries, entries_type& other )
    {
      entries.swap( other );
    }

    static inline void
    swap( rowmap_type& rowmap, rowmap_type& other )
    {
      std::swap( rowmap, other );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::FullyBuffered, bool, TOrdinal, TSize > {
    typedef sdsl::int_vector_buffer< gum::widthof< TOrdinal >::value > entries_type;
    typedef sdsl::int_vector_buffer< gum::widthof< TSize >::value > rowmap_type;
    typedef crs_matrix::FullyBuffered mutable_spec_type;

    template< typename T >
    static inline void
    init( T& a )
    {
      a = T( get_tmpfile(), std::ios::out );
    }

    template< typename T >
    static inline void
    swap( T& a, T& b )
    {
      a.swap( b );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::Compressed, bool, TOrdinal, TSize > {
    typedef sdsl::coder::elias_delta coder_type;
    typedef sdsl::enc_vector< coder_type > entries_type;
    typedef sdsl::enc_vector< coder_type > rowmap_type;
    typedef crs_matrix::Dynamic mutable_spec_type;

    template< typename T >
    static inline void
    init( T& )
    { /* NOOP */ }

    template< typename T >
    static inline void
    swap( T& a, T& b )
    {
      a.swap( b );
    }
  };

  template< typename TSpec, typename TValue, typename TOrdinal, typename TSize >
  class CRSMatrixBase;

  template< typename TSpec, typename TOrdinal, typename TSize >
  class CRSMatrixBase< TSpec, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef TSpec spec_type;
    typedef bool value_type;
    typedef TOrdinal ordinal_type;
    typedef TSize size_type;
    typedef CRSMatrixTraits< spec_type, value_type, ordinal_type, size_type > traits_type;
    typedef typename traits_type::entries_type entries_type;
    typedef typename traits_type::rowmap_type rowmap_type;
    typedef typename traits_type::mutable_spec_type mutable_spec_type;
    /* === FRIENDSHIP === */
    template< typename TSpec2, typename TValue2, typename TOrdinal2, typename TSize2 >
    friend class CRSMatrixBase;
    /* === LIFECYCLE === */
    CRSMatrixBase( ordinal_type ncols=0 ) : num_cols( ncols )
    {
      traits_type::init( this->entries );
      traits_type::init( this->rowmap );
    }
    /* === OPERATORS === */
    /**
     *  NOTE: The entries in each row must be sorted in ascending order.
     */
    inline value_type
    operator()( ordinal_type i, ordinal_type j ) const
    {
      assert( i < this->numRows() && j < this->numCols() );

      return std::binary_search( this->entries.begin() + this->rowmap[ i ],
                                 this->entries.begin() + this->rowmap[ i + 1 ], j );
    }

    /**
     *  NOTE: The entries in each row must be sorted in ascending order.
     *  NOTE: Buffered containers used in Buffered specialisations do not have constant
     *        iterator.
     */
    inline value_type
    operator()( ordinal_type i, ordinal_type j )
    {
      assert( i < this->numRows() && j < this->numCols() );

      return std::binary_search( this->entries.begin() + this->rowmap[ i ],
                                 this->entries.begin() + this->rowmap[ i + 1 ], j );
    }
    /* === METHODS === */
    inline value_type
    at( ordinal_type i, ordinal_type j ) const
    {
      if ( i >= this->numRows() || j >= this->numCols() ) {
        throw std::runtime_error( "index out of range" );
      }

      this->operator()( i, j );
    }

    /**
     *  NOTE: The `operator()` overload in Buffered specialisations is not constant.
     */
    inline value_type
    at( ordinal_type i, ordinal_type j )
    {
      if ( i >= this->numRows() || j >= this->numCols() ) {
        throw std::runtime_error( "index out of range" );
      }

      this->operator()( i, j );
    }

    template< typename TCrsMatrixBase >
    inline void
    assign( TCrsMatrixBase&& matrix )
    {
      typedef std::decay_t< decltype( matrix.entries ) > ext_entries_type;
      typedef typename ext_entries_type::value_type ext_value_type;

      gum::RandomAccessNonConstProxyContainer< ext_entries_type, uint64_t > proxy_entries(
          &matrix.entries,
          []( ext_value_type e ) { return static_cast< uint64_t >( e ); } );
      this->entries = proxy_entries;
      this->rowmap = matrix.rowmap;
      this->num_cols = matrix.num_cols;
    }

    template< typename TCrsMatrixBase >
    inline void
    swap( TCrsMatrixBase&& matrix )
    {
      traits_type::swap( this->entries, matrix.entries );
      traits_type::swap( this->rowmap, matrix.rowmap );
      this->num_cols = matrix.num_cols;
    }

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

    inline ordinal_type
    entry( size_type i ) const
    {
      return this->entries[ i ];
    }

    /**
     *  NOTE: Buffered containers used in Buffered specialisations do not have constant
     *        iterator.
     */
    inline ordinal_type
    entry( size_type i )
    {
      return this->entries[ i ];
    }

    inline size_type
    rowMap( ordinal_type i ) const
    {
      return this->rowmap[ i ];
    }

    /**
     *  NOTE: Buffered containers used in Buffered specialisations do not have constant
     *        iterator.
     */
    inline size_type
    rowMap( ordinal_type i )
    {
      return this->rowmap[ i ];
    }

    inline void
    reserve( size_type nnz )
    {
      psi::reserve( this->entries, nnz );
    }

    inline void
    shrink_to_fit( )
    {
      psi::shrink_to_fit( this->entries );
      psi::shrink_to_fit( this->rowmap );
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
  protected:
    /* === DATA MEMBERS === */
    entries_type entries;
    rowmap_type rowmap;
    ordinal_type num_cols;
    /* === METHODS === */
    /**
     *  @brief Construct the matrix by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    inline void
    build( TCrsMatrix&& ext )
    {
      psi::resize( this->entries, ext.nnz() );
      psi::resize( this->rowmap, ext.numRows() + 1 );
      this->fill_partial( ext );
    }

    /**
     *  @brief Construct the matrix block by block
     *
     *  @param[in] nrows Number of rows.
     *  @param[in] ncols Number of columns.
     *  @param[in] callback The function providing the matrix blocks.
     *  @param[in] nnz An estimation of number of non-zero values for memory reservation.
     *
     *  See `fill_partial` function documentation block.
     */
    template< typename TCallback >
    inline void
    build( ordinal_type nrows, ordinal_type ncols, TCallback callback, size_type nnz=0 )
    {
      if ( nnz ) this->reserve( nnz );
      psi::resize( this->rowmap, nrows + 1 );
      ordinal_type lrow = 0;  /* last filled row index */
      auto partial_ctor =
          [this, &lrow]( auto&& mat, ordinal_type srow=0, ordinal_type scol=0 ) {
            psi::resize( this->entries, this->nnz() + mat.nnz() );
            if ( lrow < srow  ) lrow = this->fill_zero_rows( lrow, srow );
            assert( lrow == srow );
            this->fill_partial( mat, srow, scol );
            lrow += mat.numRows();
          };
      callback( partial_ctor );
    }

    /**
     *  @brief  Populate a block in the matrix by an external matrix.
     *
     *  @param[in]  ex The external `KokkosSparse::CrsMatrix`-like matrix.
     *  @param[in]  srow The starting row of the block to be filled.
     *  @param[in]  scol The starting column of the block to be filled.
     *
     *  We call a matrix A is a block matrix when it is composed by smaller matrices,
     *  called blocks, such that all blocks are mutually disjoint. This means all
     *  non-zero values in each row or column in the matrix only belong to one of the
     *  blocks; i.e. no two rows or columns are shared between two or more blocks.
     *
     *  When a matrix is a block matrix, it can be gradually populated by blocks in
     *  parallel. For example, the matrix
     *
     *        | E1   0   0   0   0 |
     *        |  0  E2   0   0   0 |
     *    A = |  0   0  E3   0   0 |
     *        |  0   0   0  E4   0 |
     *        |  0   0   0   0  E5 |
     *
     *  can be populated by filling $E_i$ matrices separately where $i \in [1..5]$. This
     *  method fill a block represented by `ex` matrix starting at [srow][scol].
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
    fill_partial( TCrsMatrix const& ex, ordinal_type srow=0, ordinal_type scol=0 )
    {
      this->rowmap[0] = 0;
      _fill_entries_partial( ex, scol );
      _fill_rowmap_partial( ex, srow );
    }

    template< typename TCrsMatrix >
    inline void
    _fill_entries_partial( TCrsMatrix const& ex, ordinal_type scol=0 )
    {
      auto size = ex.graph.entries.extent( 0 );
      auto start = this->entries.end() - size;
      std::transform( ex.graph.entries.data(), ex.graph.entries.data() + size, start,
                      [scol]( ordinal_type e ) { return e + scol; } );
    }

    template< typename TCrsMatrix >
    inline void
    _fill_rowmap_partial( TCrsMatrix const& ex, ordinal_type srow=0 )
    {
      auto size = ex.graph.row_map.extent( 0 );
      auto start = this->rowmap.begin() + srow + 1;
      size_type snnz = this->nnz() - ex.nnz();  /* Assumed entries are already resized */
      std::transform( ex.graph.row_map.data() + 1, ex.graph.row_map.data() + size, start,
                      [snnz]( size_type e ) { return e + snnz; } );
    }

    inline ordinal_type
    fill_zero_rows( ordinal_type srow, ordinal_type erow )
    {
      if ( srow == 0 ) this->rowmap[0] = 0;
      auto last = this->rowmap.begin() + srow;
      auto end = this->rowmap.begin() + erow + 1;
      auto value = *last;
      auto start = last + 1;
      while ( start != end ) *start++ = value;
      return start - this->rowmap.begin() - 1;  /* return erow; */
    }
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::Dynamic, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::Dynamic, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::Dynamic, bool, TOrdinal, TSize > base_type;
    typedef typename base_type::spec_type spec_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::ordinal_type ordinal_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::entries_type entries_type;
    typedef typename base_type::rowmap_type rowmap_type;
    typedef typename base_type::mutable_spec_type mutable_spec_type;
    /* === FRIENDSHIP === */
    template< typename TSpec2, typename TValue2, typename TOrdinal2, typename TSize2 >
    friend class CRSMatrix;
    /* === LIFECYCLE === */
    CRSMatrix( ) : base_type( ) { }

    CRSMatrix( ordinal_type ncols, entries_type e_entries, rowmap_type e_rowmap )
      : base_type( ncols )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix const& ext )
      : base_type( ext.numCols() )
    {
      base_type::build( ext );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz=0 )
      : base_type( ncols )
    {
      base_type::build( nrows, ncols, callback, nnz );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other )
      : CRSMatrix( )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
    }

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
      : CRSMatrix( )
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
    }
    /* === OPERATORS === */
    /* copy assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix const& other )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      return *this;
    }

    /* move assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      return *this;
    }
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::Buffered, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::Buffered, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::Buffered, bool, TOrdinal, TSize > base_type;
    typedef typename base_type::spec_type spec_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::ordinal_type ordinal_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::entries_type entries_type;
    typedef typename base_type::rowmap_type rowmap_type;
    typedef typename base_type::mutable_spec_type mutable_spec_type;
    /* === FRIENDSHIP === */
    template< typename TSpec2, typename TValue2, typename TOrdinal2, typename TSize2 >
    friend class CRSMatrix;
    /* === LIFECYCLE === */
    CRSMatrix( ) : base_type( ) { }

    CRSMatrix( ordinal_type ncols, entries_type&& e_entries, rowmap_type e_rowmap )
      : base_type( ncols )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix const& ext )
      : base_type( ext.numCols() )
    {
      base_type::build( ext );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz=0 )
      : base_type( ncols )
    {
      base_type::build( nrows, ncols, callback, nnz );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other ) = delete;

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
      : CRSMatrix( )
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
    }
    /* === OPERATORS === */
    /* copy assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix const& other ) = delete;

    /* move assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      return *this;
    }
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::FullyBuffered, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::FullyBuffered, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::FullyBuffered, bool, TOrdinal, TSize > base_type;
    typedef typename base_type::spec_type spec_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::ordinal_type ordinal_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::entries_type entries_type;
    typedef typename base_type::rowmap_type rowmap_type;
    typedef typename base_type::mutable_spec_type mutable_spec_type;
    /* === FRIENDSHIP === */
    template< typename TSpec2, typename TValue2, typename TOrdinal2, typename TSize2 >
    friend class CRSMatrix;
    /* === LIFECYCLE === */
    CRSMatrix( ) : base_type( ) { }

    CRSMatrix( ordinal_type ncols, entries_type&& e_entries, rowmap_type&& e_rowmap )
      : base_type( ncols )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix const& ext )
      : base_type( ext.numCols() )
    {
      base_type::build( ext );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz=0 )
      : base_type( ncols )
    {
      base_type::build( nrows, ncols, callback, nnz );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other ) = delete;

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
      : CRSMatrix( )
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
    }
    /* === OPERATORS === */
    /* copy assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix const& other ) = delete;

    /* move assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      return *this;
    }
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::Compressed, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::Compressed, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::Compressed, bool, TOrdinal, TSize > base_type;
    typedef typename base_type::spec_type spec_type;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::ordinal_type ordinal_type;
    typedef typename base_type::size_type size_type;
    typedef typename base_type::traits_type traits_type;
    typedef typename base_type::entries_type entries_type;
    typedef typename base_type::rowmap_type rowmap_type;
    typedef typename base_type::mutable_spec_type mutable_spec_type;
    /* === FRIENDSHIP === */
    template< typename TSpec2, typename TValue2, typename TOrdinal2, typename TSize2 >
    friend class CRSMatrix;
    /* === LIFECYCLE === */
    CRSMatrix( ) : base_type( ) { }

    CRSMatrix( ordinal_type ncols, entries_type e_entries, rowmap_type e_rowmap )
      : base_type( ncols )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix, typename TMutableSpec=mutable_spec_type >
    CRSMatrix( TCrsMatrix const& ext, TMutableSpec={} )
      : base_type( )
    {
      typedef CRSMatrix< TMutableSpec, value_type, ordinal_type, size_type > mutable_type;
      mutable_type mut_mat( ext );
      this->assign( mut_mat );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback, typename TMutableSpec=mutable_spec_type >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz=0, TMutableSpec={} )
      : base_type( )
    {
      typedef CRSMatrix< TMutableSpec, value_type, ordinal_type, size_type > mutable_type;
      mutable_type mut_mat( nrows, ncols, callback, nnz );
      this->assign( mut_mat );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other )
      : CRSMatrix( )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
    }

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
      : CRSMatrix( )
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
    }
    /* === OPERATORS === */
    /* copy assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix const& other )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      return *this;
    }

    /* move assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      return *this;
    }
  };
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_CRS_MATRIX_HPP__ --- */
