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
#include <type_traits>

#include <gum/basic_types.hpp>
#include <gum/iterators.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector_buffer.hpp>
#include <sdsl/enc_vector.hpp>

#include "utils.hpp"


namespace psi {
  /* specialisation tags */
  namespace crs_matrix {
    // Group specialisation tags
    struct BasicGroup { };
    struct RangeGroup { };

    // Basic specialisation tags
    struct Dynamic { };
    struct Buffered { };
    struct FullyBuffered { };
    struct Compressed { };

    // Range specialisation tags
    struct RangeDynamic { };
    struct RangeBuffered { };
    struct RangeFullyBuffered { };
    struct RangeCompressed { };

    // Group meta-function
    template< typename TSpec >
    struct Group {
      typedef BasicGroup type;
    };

    template< >
    struct Group< RangeDynamic > {
      typedef RangeGroup type;
    };

    template< >
    struct Group< RangeBuffered > {
      typedef RangeGroup type;
    };

    template< >
    struct Group< RangeFullyBuffered > {
      typedef RangeGroup type;
    };

    template< >
    struct Group< RangeCompressed > {
      typedef RangeGroup type;
    };

    /* === Specialised helper functions === */

    template< typename TIter, typename TOrdinal >
    inline bool
    binary_search( TIter begin, TIter end, TOrdinal key, RangeGroup /*tag*/ )
    {
      auto first = begin;
      auto len = end - first;

      while ( len > 0 ) {
        auto half = len / 2;
        auto middle = first;
        middle += half;
        if ( *middle < key ) {
          first = middle;
          ++first;
          len = len - half - 1;
	      } else len = half;
      }
      auto idx = first - begin;
      if ( idx % 2 ) --first;
      return first != end && *first <= key;
    }

    /**
     *  NOTE: The input distance indices are passed by non-const references, since
     *        containers in const Buffered specialisations cannot be iterated.
     */
    template< typename TEntries, typename TRowmap >
    inline typename TRowmap::value_type
    nnz( TEntries& entries, TRowmap& rowmap, RangeGroup /*tag*/ )
    {
      typedef typename TEntries::value_type ordinal_type;
      typedef typename TRowmap::value_type size_type;

      size_type nnz_counter = 0;
      size_type idx = 0;
      for ( ordinal_type row_i = 1; row_i < rowmap.size(); ++row_i ) {
        for ( ; idx < rowmap[ row_i ]; idx += 2 ) {
          assert( idx+1 < rowmap[ row_i ] );
          nnz_counter += entries[ idx + 1 ] - entries[ idx ] + 1;
        }
      }
      return nnz_counter;
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TEntries,
              typename TRowmap,
              typename TFromEntriesIter,
              typename TFromRowmapIter >
    inline void
    append( TEntries& entries,
            TRowmap& rowmap,
            TFromEntriesIter e_entries_begin,
            TFromRowmapIter e_rowmap_begin,
            TFromRowmapIter e_rowmap_end,
            RangeGroup /*to*/,
            BasicGroup /*from*/,
            typename TEntries::value_type scol=0,
            typename TEntries::value_type srow=0 )
    {
      typedef typename TRowmap::value_type size_type;

      size_type idx = 0;
      size_type last_value_pp = 0;  // last value + 1: since there is zero values.
      rowmap[ 0 ] = 0;
      auto rowmap_itr = rowmap.begin() + srow + 1;
      ++e_rowmap_begin;  // skip the first element which always is zero
      while ( e_rowmap_begin != e_rowmap_end ) {
        while ( idx < *e_rowmap_begin ) {
          if ( last_value_pp &&
               static_cast< size_type >( *e_entries_begin ) == last_value_pp ) last_value_pp++;
          else {
            if ( last_value_pp ) entries.push_back( last_value_pp - 1 + scol );
            entries.push_back( *e_entries_begin + scol );
            last_value_pp = *e_entries_begin + 1;
          }
          ++idx;
          ++e_entries_begin;
        }
        if ( last_value_pp ) {
          entries.push_back( last_value_pp - 1 + scol );
          last_value_pp = 0;
        }
        *rowmap_itr++ = entries.size();
        ++e_rowmap_begin;
      }
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TEntries,
              typename TRowmap,
              typename TFromEntriesIter,
              typename TFromRowmapIter >
    inline void
    append( TEntries& entries,
            TRowmap& rowmap,
            TFromEntriesIter e_entries_begin,
            TFromRowmapIter e_rowmap_begin,
            TFromRowmapIter e_rowmap_end,
            BasicGroup /*to*/,
            RangeGroup /*from*/,
            typename TEntries::value_type scol=0,
            typename TEntries::value_type srow=0 )
    {
      typedef typename TRowmap::value_type size_type;

      size_type idx = 0;
      rowmap[ 0 ] = 0;
      auto rowmap_itr = rowmap.begin() + srow + 1;
      ++e_rowmap_begin;  // skip the first element which always is zero
      while ( e_rowmap_begin != e_rowmap_end ) {
        while ( idx < *e_rowmap_begin ) {
          assert( idx+1 < *e_rowmap_begin );
          for ( auto elem = *e_entries_begin; elem <= *( e_entries_begin + 1 ); ++elem ) {
            entries.push_back( elem + scol );
          }
          idx += 2;
          e_entries_begin += 2;
        }
        *rowmap_itr++ = entries.size();
        ++e_rowmap_begin;
      }
    }

    template< typename TEntries, typename TRowmap, typename TCrsMatrix, typename TOrdinal >
    inline void
    fill_all_partial( TEntries& entries, TRowmap& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0, BasicGroup={} /*tag*/ )
    {
      typedef typename TRowmap::value_type size_type;

      const size_type snnz = entries.size();
      auto ent_size = ex.graph.entries.extent( 0 );
      assert( ent_size == ex.nnz() );
      psi::resize( entries, snnz + ent_size );
      auto ent_start = entries.begin() + snnz;
      std::transform( ex.graph.entries.data(), ex.graph.entries.data() + ent_size,
                      ent_start, [scol]( TOrdinal e ) { return e + scol; } );

      auto rowmap_size = ex.graph.row_map.extent( 0 );
      auto row_start = rowmap.begin() + srow + 1;
      std::transform( ex.graph.row_map.data() + 1, ex.graph.row_map.data() + rowmap_size,
                      row_start, [snnz]( size_type e ) { return e + snnz; } );
    }

    template< typename TEntries, typename TRowmap, typename TCrsMatrix, typename TOrdinal >
    inline void
    fill_all_partial( TEntries& entries, TRowmap& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0, RangeGroup={} /*tag*/ )
    {
      auto rowmap_size = ex.graph.row_map.extent( 0 );
      append( entries, rowmap, ex.graph.entries.data(), ex.graph.row_map.data(),
              ex.graph.row_map.data() + rowmap_size, RangeGroup{}, BasicGroup{}, scol, srow );
    }
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

  template< typename ...TArgs >
  using RangeDynamicCRSMatrix = CRSMatrix< crs_matrix::RangeDynamic, TArgs... >;

  template< typename ...TArgs >
  using RangeBufferedCRSMatrix = CRSMatrix< crs_matrix::RangeBuffered, TArgs... >;

  template< typename ...TArgs >
  using RangeFullyBufferedCRSMatrix = CRSMatrix< crs_matrix::RangeFullyBuffered, TArgs... >;

  template< typename ...TArgs >
  using RangeCompressedCRSMatrix = CRSMatrix< crs_matrix::RangeCompressed, TArgs... >;

  template< typename TNativeCRSMatrix, typename TGroup >
  struct _make_buffered_impl;

  template< typename TSpec, typename ...TArgs >
  struct _make_buffered_impl< CRSMatrix< TSpec, TArgs... >, crs_matrix::RangeGroup > {
  private:
    typedef typename crs_matrix::Group< CRSMatrix< TSpec, TArgs... > >::type group_type;
    static_assert( std::is_same< group_type, crs_matrix::RangeGroup >::value, "should be in Range group" );
  public:
    typedef CRSMatrix< crs_matrix::RangeBuffered, TArgs... > type;
  };

  template< typename TSpec, typename ...TArgs >
  struct _make_buffered_impl< CRSMatrix< TSpec, TArgs... >, crs_matrix::BasicGroup > {
  private:
    typedef typename crs_matrix::Group< CRSMatrix< TSpec, TArgs... > >::type group_type;
    static_assert( std::is_same< group_type, crs_matrix::BasicGroup >::value, "should be in Basic group" );
  public:
    typedef CRSMatrix< crs_matrix::Buffered, TArgs... > type;
  };

  template< typename TNativeCRSMatrix >
  struct make_buffered {
  private:
    typedef typename crs_matrix::Group< TNativeCRSMatrix >::type group_type;
  public:
    typedef typename _make_buffered_impl< TNativeCRSMatrix, group_type >::type type;
  };

  template< typename TNativeCRSMatrix >
  using make_buffered_t = typename make_buffered< TNativeCRSMatrix >::type;

  template< typename TNativeCRSMatrix, typename TGroup >
  struct _make_fully_buffered_impl;

  template< typename TSpec, typename ...TArgs >
  struct _make_fully_buffered_impl< CRSMatrix< TSpec, TArgs... >, crs_matrix::RangeGroup > {
  private:
    typedef typename crs_matrix::Group< CRSMatrix< TSpec, TArgs... > >::type group_type;
    static_assert( std::is_same< group_type, crs_matrix::RangeGroup >::value, "should be in Range group" );
  public:
    typedef CRSMatrix< crs_matrix::RangeFullyBuffered, TArgs... > type;
  };

  template< typename TSpec, typename ...TArgs >
  struct _make_fully_buffered_impl< CRSMatrix< TSpec, TArgs... >, crs_matrix::BasicGroup > {
  private:
    typedef typename crs_matrix::Group< CRSMatrix< TSpec, TArgs... > >::type group_type;
    static_assert( std::is_same< group_type, crs_matrix::BasicGroup >::value, "should be in Basic group" );
  public:
    typedef CRSMatrix< crs_matrix::FullyBuffered, TArgs... > type;
  };

  template< typename TNativeCRSMatrix >
  struct make_fully_buffered {
  private:
    typedef typename crs_matrix::Group< TNativeCRSMatrix >::type group_type;
  public:
    typedef typename _make_fully_buffered_impl< TNativeCRSMatrix, group_type >::type type;
  };

  template< typename TNativeCRSMatrix >
  using make_fully_buffered_t = typename make_fully_buffered< TNativeCRSMatrix >::type;

  template< typename TNativeCRSMatrix >
  struct make_range;

  template< typename ...TArgs >
  struct make_range< CRSMatrix< crs_matrix::Dynamic, TArgs... > > {
    typedef CRSMatrix< crs_matrix::RangeDynamic, TArgs... > type;
  };

  template< typename ...TArgs >
  struct make_range< CRSMatrix< crs_matrix::Buffered, TArgs... > > {
    typedef CRSMatrix< crs_matrix::RangeBuffered, TArgs... > type;
  };

  template< typename ...TArgs >
  struct make_range< CRSMatrix< crs_matrix::FullyBuffered, TArgs... > > {
    typedef CRSMatrix< crs_matrix::RangeFullyBuffered, TArgs... > type;
  };

  template< typename ...TArgs >
  struct make_range< CRSMatrix< crs_matrix::Compressed, TArgs... > > {
    typedef CRSMatrix< crs_matrix::RangeCompressed, TArgs... > type;
  };

  template< typename TNativeCRSMatrix >
  using make_range_t = typename make_range< TNativeCRSMatrix >::type;

  template< typename TNativeCRSMatrix >
  struct make_basic;

  template< typename ...TArgs >
  struct make_basic< CRSMatrix< crs_matrix::RangeDynamic, TArgs... > > {
    typedef CRSMatrix< crs_matrix::Dynamic, TArgs... > type;
  };

  template< typename ...TArgs >
  struct make_basic< CRSMatrix< crs_matrix::RangeBuffered, TArgs... > > {
    typedef CRSMatrix< crs_matrix::Buffered, TArgs... > type;
  };

  template< typename ...TArgs >
  struct make_basic< CRSMatrix< crs_matrix::RangeFullyBuffered, TArgs... > > {
    typedef CRSMatrix< crs_matrix::FullyBuffered, TArgs... > type;
  };

  template< typename ...TArgs >
  struct make_basic< CRSMatrix< crs_matrix::RangeCompressed, TArgs... > > {
    typedef CRSMatrix< crs_matrix::Compressed, TArgs... > type;
  };

  template< typename TNativeCRSMatrix >
  using make_basic_t = typename make_basic< TNativeCRSMatrix >::type;

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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return std::binary_search( begin, end, key );
    }

    static inline TSize
    nnz( entries_type const& entries, rowmap_type const& )
    {
      return entries.size();
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      crs_matrix::append( entries, rowmap, e_entries.begin(), e_rowmap.begin(), e_rowmap.end(),
                          crs_matrix::BasicGroup{}, crs_matrix::RangeGroup{} );
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      crs_matrix::fill_all_partial( entries, rowmap, ex, scol, srow, crs_matrix::BasicGroup{} );
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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return std::binary_search( begin, end, key );
    }

    static inline TSize
    nnz( entries_type const& entries, rowmap_type const& )
    {
      return entries.size();
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      crs_matrix::append( entries, rowmap, e_entries.begin(), e_rowmap.begin(), e_rowmap.end(),
                          crs_matrix::BasicGroup{}, crs_matrix::RangeGroup{} );
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      crs_matrix::fill_all_partial( entries, rowmap, ex, scol, srow, crs_matrix::BasicGroup{} );
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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return std::binary_search( begin, end, key );
    }

    static inline TSize
    nnz( entries_type const& entries, rowmap_type const& )
    {
      return entries.size();
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      crs_matrix::append( entries, rowmap, e_entries.begin(), e_rowmap.begin(), e_rowmap.end(),
                          crs_matrix::BasicGroup{}, crs_matrix::RangeGroup{} );
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      crs_matrix::fill_all_partial( entries, rowmap, ex, scol, srow, crs_matrix::BasicGroup{} );
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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return std::binary_search( begin, end, key );
    }

    static inline TSize
    nnz( entries_type const& entries, rowmap_type const& )
    {
      return entries.size();
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      typedef CRSMatrixTraits< mutable_spec_type, bool, TOrdinal, TSize > mutable_trait_type;
      typedef typename mutable_trait_type::entries_type m_entries_type;
      typedef typename mutable_trait_type::rowmap_type m_rowmap_type;

      m_entries_type m_entries;
      m_rowmap_type m_rowmap;
      mutable_trait_type::init( m_entries );
      mutable_trait_type::init( m_rowmap );
      m_rowmap.resize( e_rowmap.size() );
      m_rowmap[ 0 ] = 0;
      mutable_trait_type::assign( m_entries, m_rowmap, e_entries, e_rowmap, crs_matrix::RangeGroup{} );

      gum::RandomAccessNonConstProxyContainer< m_entries_type, uint64_t > proxy_entries(
          &m_entries,
          []( typename m_entries_type::value_type e ) { return static_cast< uint64_t >( e ); } );
      assign( entries, rowmap, proxy_entries, m_rowmap, crs_matrix::BasicGroup{} );
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      throw std::runtime_error( "A Compressed Basic CRS cannot be modified" );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::RangeDynamic, bool, TOrdinal, TSize > {
    typedef std::vector< TOrdinal > entries_type;
    typedef std::vector< TSize > rowmap_type;
    typedef crs_matrix::RangeDynamic mutable_spec_type;

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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return crs_matrix::binary_search( begin, end, key, crs_matrix::RangeGroup{} );
    }

    static inline TSize
    nnz( entries_type const& entries, rowmap_type const& rowmap )
    {
      return crs_matrix::nnz( entries, rowmap, crs_matrix::RangeGroup{} );
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      crs_matrix::append( entries, rowmap, e_entries.begin(), e_rowmap.begin(), e_rowmap.end(),
                          crs_matrix::RangeGroup{}, crs_matrix::BasicGroup{} );
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      crs_matrix::fill_all_partial( entries, rowmap, ex, scol, srow, crs_matrix::RangeGroup{} );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::RangeBuffered, bool, TOrdinal, TSize > {
    typedef sdsl::int_vector_buffer< gum::widthof< TOrdinal >::value > entries_type;
    typedef std::vector< TSize > rowmap_type;
    typedef crs_matrix::RangeBuffered mutable_spec_type;

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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return crs_matrix::binary_search( begin, end, key, crs_matrix::RangeGroup{} );
    }

    static inline TSize
    nnz( entries_type& entries, rowmap_type& rowmap )
    {
      return crs_matrix::nnz( entries, rowmap, crs_matrix::RangeGroup{} );
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      crs_matrix::append( entries, rowmap, e_entries.begin(), e_rowmap.begin(), e_rowmap.end(),
                          crs_matrix::RangeGroup{}, crs_matrix::BasicGroup{} );
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      crs_matrix::fill_all_partial( entries, rowmap, ex, scol, srow, crs_matrix::RangeGroup{} );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::RangeFullyBuffered, bool, TOrdinal, TSize > {
    typedef sdsl::int_vector_buffer< gum::widthof< TOrdinal >::value > entries_type;
    typedef sdsl::int_vector_buffer< gum::widthof< TSize >::value > rowmap_type;
    typedef crs_matrix::RangeFullyBuffered mutable_spec_type;

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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return crs_matrix::binary_search( begin, end, key, crs_matrix::RangeGroup{} );
    }

    static inline TSize
    nnz( entries_type& entries, rowmap_type& rowmap )
    {
      return crs_matrix::nnz( entries, rowmap, crs_matrix::RangeGroup{} );
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      crs_matrix::append( entries, rowmap, e_entries.begin(), e_rowmap.begin(), e_rowmap.end(),
                          crs_matrix::RangeGroup{}, crs_matrix::BasicGroup{} );
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      crs_matrix::fill_all_partial( entries, rowmap, ex, scol, srow, crs_matrix::RangeGroup{} );
    }
  };

  template< typename TOrdinal, typename TSize >
  struct CRSMatrixTraits< crs_matrix::RangeCompressed, bool, TOrdinal, TSize > {
    typedef sdsl::coder::elias_delta coder_type;
    typedef sdsl::enc_vector< coder_type > entries_type;
    typedef sdsl::enc_vector< coder_type > rowmap_type;
    typedef crs_matrix::RangeDynamic mutable_spec_type;

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

    template< typename TIter >
    static inline bool  /* value_type */
    binary_search( TIter begin, TIter end, TOrdinal key )
    {
      return crs_matrix::binary_search( begin, end, key, crs_matrix::RangeGroup{} );
    }

    static inline TSize
    nnz( entries_type const& entries, rowmap_type const& rowmap )
    {
      return crs_matrix::nnz( entries, rowmap, crs_matrix::RangeGroup{} );
    }

    /**
     *  NOTE: This function assumes that `rowmap` is an already allocated array of
     *        length `numRows`+1 with `rowmap[0] == 0`.
     */
    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::BasicGroup )
    {
      typedef CRSMatrixTraits< mutable_spec_type, bool, TOrdinal, TSize > mutable_trait_type;
      typedef typename mutable_trait_type::entries_type m_entries_type;
      typedef typename mutable_trait_type::rowmap_type m_rowmap_type;

      m_entries_type m_entries;
      m_rowmap_type m_rowmap;
      mutable_trait_type::init( m_entries );
      mutable_trait_type::init( m_rowmap );
      m_rowmap.resize( e_rowmap.size() );
      m_rowmap[ 0 ] = 0;
      mutable_trait_type::assign( m_entries, m_rowmap, e_entries, e_rowmap, crs_matrix::BasicGroup{} );

      gum::RandomAccessNonConstProxyContainer< m_entries_type, uint64_t > proxy_entries(
          &m_entries,
          []( typename m_entries_type::value_type e ) { return static_cast< uint64_t >( e ); } );
      assign( entries, rowmap, proxy_entries, m_rowmap, crs_matrix::RangeGroup{} );
    }

    template< typename TExtEntries, typename TExtRowmap >
    static inline void
    assign( entries_type& entries, rowmap_type& rowmap, TExtEntries& e_entries,
            TExtRowmap& e_rowmap, crs_matrix::RangeGroup )
    {
      entries = e_entries;
      rowmap = e_rowmap;
    }

    template< typename TCrsMatrix >
    static inline void
    fill_all_partial( entries_type& entries, rowmap_type& rowmap, TCrsMatrix const& ex,
                      TOrdinal scol=0, TOrdinal srow=0 )
    {
      throw std::runtime_error( "A Compressed Range CRS cannot be modified" );
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

      return traits_type::binary_search( this->entries.begin() + this->rowmap[ i ],
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

      return traits_type::binary_search( this->entries.begin() + this->rowmap[ i ],
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
      typedef typename std::decay_t< TCrsMatrixBase >::spec_type ext_spec_type;
      typedef typename crs_matrix::Group< ext_spec_type >::type ext_group_type;

      gum::RandomAccessNonConstProxyContainer< ext_entries_type, uint64_t > proxy_entries(
          &matrix.entries,
          []( ext_value_type e ) { return static_cast< uint64_t >( e ); } );
      psi::resize( this->rowmap, matrix.numRows() + 1 );
      psi::clear( this->entries );
      traits_type::assign( this->entries, this->rowmap, proxy_entries, matrix.rowmap,
                           ext_group_type{} );
      this->num_cols = matrix.num_cols;
    }

    template< typename TCrsMatrixBase >
    inline void
    swap( TCrsMatrixBase&& matrix )
    {
      traits_type::swap( this->entries, matrix.entries );
      traits_type::swap( this->rowmap, matrix.rowmap );
      std::swap( this->num_cols, matrix.num_cols );
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
      return traits_type::nnz( this->entries, this->rowmap );
    }

    /**
     *  NOTE: Buffered containers used in Buffered specialisations need to be non-const.
     */
    inline size_type
    nnz( )
    {
      return traits_type::nnz( this->entries, this->rowmap );
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
    reserve( size_type nnz_est )
    {
      psi::reserve( this->entries, nnz_est );
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
     *
     *  @return the number of non-zero values in the final matrix.
     */
    template< typename TCrsMatrix >
    inline size_type
    build( TCrsMatrix&& ext )
    {
      psi::resize( this->rowmap, ext.numRows() + 1 );
      this->fill_partial( ext );
      return ext.nnz();
    }

    /**
     *  @brief Construct the matrix block by block
     *
     *  @param[in] nrows Number of rows.
     *  @param[in] ncols Number of columns.
     *  @param[in] callback The function providing the matrix blocks.
     *  @param[in] nnz_est An estimation of number of non-zero values for memory reservation.
     *
     *  @return the number of non-zero values in the final matrix.
     *
     *  See `fill_partial` function documentation block.
     */
    template< typename TCallback >
    inline size_type
    build( ordinal_type nrows, ordinal_type ncols, TCallback callback, size_type nnz_est=0 )
    {
      if ( nnz_est ) this->reserve( nnz_est );
      psi::resize( this->rowmap, nrows + 1 );
      ordinal_type lrow = 0;  /* last filled row index */
      size_type tnnz = 0;
      auto partial_ctor =
          [this, &lrow, &tnnz]( auto&& mat, ordinal_type srow=0, ordinal_type scol=0 ) {
            if ( lrow < srow ) lrow = this->fill_zero_rows( lrow, srow );
            assert( lrow == srow );
            this->fill_partial( mat, srow, scol );
            lrow += mat.numRows();
            tnnz += mat.nnz();
          };
      callback( partial_ctor );
      return tnnz;
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
      traits_type::fill_all_partial( this->entries, this->rowmap, ex, scol, srow );
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
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0 )
      : base_type( ncols )
    {
      base_type::build( nrows, ncols, callback, nnz_est );
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
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0 )
      : base_type( ncols )
    {
      base_type::build( nrows, ncols, callback, nnz_est );
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
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0 )
      : base_type( ncols )
    {
      base_type::build( nrows, ncols, callback, nnz_est );
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
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback, typename TMutableSpec=mutable_spec_type >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0, TMutableSpec={} )
      : base_type( )
    {
      typedef CRSMatrix< TMutableSpec, value_type, ordinal_type, size_type > mutable_type;
      mutable_type mut_mat( nrows, ncols, callback, nnz_est );
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

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::RangeDynamic, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::RangeDynamic, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::RangeDynamic, bool, TOrdinal, TSize > base_type;
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
    CRSMatrix( ) : base_type( ), m_nnz( 0 ) { }

    CRSMatrix( ordinal_type ncols, entries_type e_entries, rowmap_type e_rowmap,
               size_type e_nnz=0 )
      : base_type( ncols ), m_nnz( e_nnz )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
      if ( this->m_nnz == 0 ) this->m_nnz = base_type::nnz();
      assert( this->entries.size() == 0 || this->m_nnz != 0 );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix const& ext )
      : base_type( ext.numCols() ), m_nnz( ext.nnz() )
    {
      base_type::build( ext );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0 )
      : base_type( ncols )
    {
      this->m_nnz = base_type::build( nrows, ncols, callback, nnz_est );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other )
      : CRSMatrix( )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
    }

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
      : CRSMatrix( )
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
    }
    /* === OPERATORS === */
    /* copy assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix const& other )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
      return *this;
    }

    /* move assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
      return *this;
    }
    /* === METHODS === */
    template< typename TCrsMatrix >
    inline void
    assign( TCrsMatrix&& matrix )
    {
      base_type::assign( matrix );
      this->m_nnz = matrix.nnz();
    }

    inline size_type
    nnz() const
    {
      return this->m_nnz;
    }

    inline void
    serialize( std::ostream& out ) const
    {
      base_type::serialize( out );
      psi::serialize( out, static_cast< uint64_t >( this->m_nnz ) );
    }

    inline void
    load( std::istream& in )
    {
      base_type::load( in );
      uint64_t dnnz;
      psi::deserialize( in, dnnz );
      this->m_nnz = dnnz;
    }
  private:
    /* === DATA MEMBER === */
    size_type m_nnz;
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::RangeBuffered, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::RangeBuffered, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::RangeBuffered, bool, TOrdinal, TSize > base_type;
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
    CRSMatrix( ) : base_type( ), m_nnz( 0 ) { }

    CRSMatrix( ordinal_type ncols, entries_type&& e_entries, rowmap_type e_rowmap,
               size_type e_nnz=0 )
      : base_type( ncols ), m_nnz( e_nnz )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
      if ( this->m_nnz == 0 ) this->m_nnz = base_type::nnz();
      assert( this->entries.size() == 0 || this->m_nnz != 0 );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix const& ext )
      : base_type( ext.numCols() ), m_nnz( ext.nnz() )
    {
      base_type::build( ext );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0 )
      : base_type( ncols )
    {
      this->m_nnz = base_type::build( nrows, ncols, callback, nnz_est );
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
      this->m_nnz = other.m_nnz;
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
      this->m_nnz = other.m_nnz;
      return *this;
    }
    /* === METHODS === */
    template< typename TCrsMatrix >
    inline void
    assign( TCrsMatrix&& matrix )
    {
      base_type::assign( matrix );
      this->m_nnz = matrix.nnz();
    }

    inline size_type
    nnz( )
    {
      return this->m_nnz;
    }

    inline void
    serialize( std::ostream& out ) const
    {
      base_type::serialize( out );
      psi::serialize( out, static_cast< uint64_t >( this->m_nnz ) );
    }

    inline void
    load( std::istream& in )
    {
      base_type::load( in );
      uint64_t dnnz;
      psi::deserialize( in, dnnz );
      this->m_nnz = dnnz;
    }
  private:
    /* === DATA MEMBER === */
    size_type m_nnz;
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::RangeFullyBuffered, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::RangeFullyBuffered, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::RangeFullyBuffered, bool, TOrdinal, TSize > base_type;
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
    CRSMatrix( ) : base_type( ), m_nnz( 0 ) { }

    CRSMatrix( ordinal_type ncols, entries_type&& e_entries, rowmap_type&& e_rowmap,
               size_type e_nnz=0 )
      : base_type( ncols ), m_nnz( e_nnz )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
      if ( this->m_nnz == 0 ) this->m_nnz = base_type::nnz();
      assert( this->entries.size() == 0 || this->m_nnz != 0 );
    }

    /**
     *  @brief Construct by an external `KokkosSparse::CrsMatrix`-like matrix.
     */
    template< typename TCrsMatrix >
    CRSMatrix( TCrsMatrix const& ext )
      : base_type( ext.numCols() ), m_nnz( ext.nnz() )
    {
      base_type::build( ext );
    }

    /**
     *  @brief Construct by a `KokkosSparse::CrsMatrix`-like matrix part by part.
     *
     *  @param[in] nrows
     *  @param[in] ncols
     *  @param[in] callback
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0 )
      : base_type( ncols )
    {
      this->m_nnz = base_type::build( nrows, ncols, callback, nnz_est );
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
      this->m_nnz = other.m_nnz;
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
      this->m_nnz = other.m_nnz;
      return *this;
    }
    /* === METHODS === */
    template< typename TCrsMatrix >
    inline void
    assign( TCrsMatrix&& matrix )
    {
      base_type::assign( matrix );
      this->m_nnz = matrix.nnz();
    }

    inline size_type
    nnz( )
    {
      return this->m_nnz;
    }

    inline void
    serialize( std::ostream& out ) const
    {
      base_type::serialize( out );
      psi::serialize( out, static_cast< uint64_t >( this->m_nnz ) );
    }

    inline void
    load( std::istream& in )
    {
      base_type::load( in );
      uint64_t dnnz;
      psi::deserialize( in, dnnz );
      this->m_nnz = dnnz;
    }
  private:
    /* === DATA MEMBER === */
    size_type m_nnz;
  };

  template< typename TOrdinal, typename TSize >
  class CRSMatrix< crs_matrix::RangeCompressed, bool, TOrdinal, TSize >
    : public CRSMatrixBase< crs_matrix::RangeCompressed, bool, TOrdinal, TSize > {
  public:
    /* === TYPE MEMBERS === */
    typedef CRSMatrixBase< crs_matrix::RangeCompressed, bool, TOrdinal, TSize > base_type;
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
    CRSMatrix( ) : base_type( ), m_nnz( 0 ) { }

    CRSMatrix( ordinal_type ncols, entries_type e_entries, rowmap_type e_rowmap,
               size_type e_nnz=0 )
      : base_type( ncols ), m_nnz( e_nnz )
    {
      this->entries = std::move( e_entries );
      this->rowmap = std::move( e_rowmap );
      if ( this->m_nnz == 0 ) this->m_nnz = base_type::nnz();
      assert( this->entries.size() == 0 || this->m_nnz != 0 );
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
     *  @param[in] nnz_est An estimation of the number of non-zero values for reserving memory.
     */
    template< typename TCallback, typename TMutableSpec=mutable_spec_type >
    CRSMatrix( ordinal_type nrows, ordinal_type ncols, TCallback callback,
               size_type nnz_est=0, TMutableSpec={} )
      : base_type( )
    {
      typedef CRSMatrix< TMutableSpec, value_type, ordinal_type, size_type > mutable_type;
      mutable_type mut_mat( nrows, ncols, callback, nnz_est );
      this->assign( mut_mat );
    }

    /* copy constructor */
    CRSMatrix( CRSMatrix const& other )
      : CRSMatrix( )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
    }

    /* move constructor */
    CRSMatrix( CRSMatrix&& other ) noexcept
      : CRSMatrix( )
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
    }
    /* === OPERATORS === */
    /* copy assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix const& other )
    {
      this->entries = other.entries;
      this->rowmap = other.rowmap;
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
      return *this;
    }

    /* move assignment operator */
    inline CRSMatrix&
    operator=( CRSMatrix&& other ) noexcept
    {
      this->entries = std::move( other.entries );
      this->rowmap = std::move( other.rowmap );
      this->num_cols = other.num_cols;
      this->m_nnz = other.m_nnz;
      return *this;
    }
    /* === METHODS === */
    template< typename TCrsMatrix >
    inline void
    assign( TCrsMatrix&& matrix )
    {
      base_type::assign( matrix );
      this->m_nnz = matrix.nnz();
    }

    inline size_type
    nnz() const
    {
      return this->m_nnz;
    }

    inline void
    serialize( std::ostream& out ) const
    {
      base_type::serialize( out );
      psi::serialize( out, static_cast< uint64_t >( this->m_nnz ) );
    }

    inline void
    load( std::istream& in )
    {
      base_type::load( in );
      uint64_t dnnz;
      psi::deserialize( in, dnnz );
      this->m_nnz = dnnz;
    }
  private:
    /* === DATA MEMBER === */
    size_type m_nnz;
  };

  /* === interface functions === */

  /**
   *  @brief  Merge two distance indices (Range Group)
   *
   *  @param  dindex1 first distance index
   *  @param  dindex2 second distance index
   *  @return a mutable merged distance index of type `TMutableCRSMatrix`
   *
   *  NOTE: The resulting mutable matrix can be assigned to a immutable compressed
   *        matrix afterwards.
   *
   *  NOTE: The input distance indices are passed by non-const references, since
   *        containers in const Buffered specialisations cannot be iterated.
   */
  template< typename TMutableCRSMatrix, typename TCRSMatrix >
  inline TMutableCRSMatrix
  merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
                        crs_matrix::RangeGroup /* tag */ )
  {
    typedef TMutableCRSMatrix crsmat_mutable_type;
    typedef TCRSMatrix crsmat_type;
    typedef typename crsmat_type::ordinal_type ordinal_type;
    typedef typename crsmat_type::size_type size_type;

    typename crsmat_mutable_type::entries_type entries;
    typename crsmat_mutable_type::rowmap_type rowmap;
    crsmat_mutable_type::base_type::traits_type::init( entries );
    crsmat_mutable_type::base_type::traits_type::init( rowmap );
    size_type cursor1 = 0;    // current entry index in the first distance index
    size_type cursor2 = 0;    // current entry index in the second distance index
    size_type end1;  // last entry index of the row in the first distance index
    size_type end2;  // last entry index of the row in the second distance index

    assert( dindex1.numRows() == dindex2.numRows() );
    assert( dindex1.numCols() == dindex2.numCols() );

    auto nof_rows = dindex1.numRows();
    auto nof_cols = dindex1.numCols();

    /* NOTE: this function assumes that `cursors` are in the entry ranges. */
    auto fetch_min_and_adv =
      []( TCRSMatrix& dindex1, TCRSMatrix& dindex2, size_type& cursor1, size_type& cursor2 )
      {
        ordinal_type l1 = dindex1.entry( cursor1 );
        ordinal_type u1 = dindex1.entry( cursor1 + 1 );
        ordinal_type l2 = dindex2.entry( cursor2 );
        ordinal_type u2 = dindex2.entry( cursor2 + 1 );

        if ( l1 < l2 || ( l1 == l2 && u1 < u2 ) ) {
          ++cursor1; ++cursor1;
          return std::make_pair( l1, u1 );
        }
        else {
          ++cursor2; ++cursor2;
          return std::make_pair( l2, u2 );
        }
      };

    auto merge_and_adv =
      []( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
          size_type& cursor1, size_type end1,
          size_type& cursor2, size_type end2,
          ordinal_type l, ordinal_type u )
      {
        ordinal_type l1 = std::numeric_limits< ordinal_type >::max();
        ordinal_type u1 = std::numeric_limits< ordinal_type >::max();
        ordinal_type l2 = std::numeric_limits< ordinal_type >::max();
        ordinal_type u2 = std::numeric_limits< ordinal_type >::max();

        if ( cursor1 < end1 ) {
          l1 = dindex1.entry( cursor1 );
          u1 = dindex1.entry( cursor1 + 1 );
        }

        if ( cursor2 < end2 ) {
          l2 = dindex2.entry( cursor2 );
          u2 = dindex2.entry( cursor2 + 1 );
        }

        while ( true ) {
          if ( u + 1 >= l1 ) {
            l = std::min( l, l1 );
            u = std::max( u, u1 );
            ++cursor1; ++cursor1;
            if ( cursor1 < end1 ) {
              l1 = dindex1.entry( cursor1 );
              u1 = dindex1.entry( cursor1 + 1 );
            }
            else l1 = std::numeric_limits< ordinal_type >::max();
          }
          else if ( u + 1 >= l2 ) {
            l = std::min( l, l2 );
            u = std::max( u, u2 );
            ++cursor2; ++cursor2;
            if ( cursor2 < end2 ) {
              l2 = dindex2.entry( cursor2 );
              u2 = dindex2.entry( cursor2 + 1 );
            }
            else l2 = std::numeric_limits< ordinal_type >::max();
          }
          else break;
        }

        return std::make_pair( l, u );
      };

    size_type c_nnz = 0;
    ordinal_type l = 0;
    ordinal_type u = 0;
    for ( ordinal_type nrow = 0; nrow < nof_rows; ++nrow ) {
      rowmap.push_back( entries.size() );

      end1 = dindex1.rowMap( nrow + 1 );
      end2 = dindex2.rowMap( nrow + 1 );
      l = 0;
      u = 0;
      while ( cursor1 < end1 ) {
        if ( cursor2 >= end2 ) {
          while ( cursor1 < end1 ) {
            l = dindex1.entry( cursor1++ );
            u = dindex1.entry( cursor1++ );
            entries.push_back( l );
            entries.push_back( u );
            c_nnz += u - l + 1;
          }
          break;
        }
        std::tie( l, u ) = fetch_min_and_adv( dindex1, dindex2, cursor1, cursor2 );
        std::tie( l, u ) =
          merge_and_adv( dindex1, dindex2, cursor1, end1, cursor2, end2, l, u );
        entries.push_back( l );
        entries.push_back( u );
        c_nnz += u - l + 1;
      }
      while ( cursor2 < end2 ) {
        l = dindex2.entry( cursor2++ );
        u = dindex2.entry( cursor2++ );
        entries.push_back( l );
        entries.push_back( u );
        c_nnz += u - l + 1;
      }
    }
    rowmap.push_back( entries.size() );

    return crsmat_mutable_type( nof_cols, std::move( entries ), std::move( rowmap ), c_nnz );
  }

  /**
   *  @brief  Merge two distance indices (Basic Group)
   *
   *  @param  dindex1 first distance index
   *  @param  dindex2 second distance index
   *  @return a mutable merged distance index of type `TMutableCRSMatrix`
   *
   *  NOTE: The resulting mutable matrix can be assigned to a immutable compressed
   *        matrix afterwards.
   *
   *  NOTE: The input distance indices are passed by non-const references, since
   *        containers in const Buffered specialisations cannot be iterated.
   */
  template< typename TMutableCRSMatrix, typename TCRSMatrix >
  inline TMutableCRSMatrix
  merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
                        crs_matrix::BasicGroup /* tag */ )
  {
    typedef TMutableCRSMatrix crsmat_mutable_type;
    typedef TCRSMatrix crsmat_type;
    typedef typename crsmat_type::ordinal_type ordinal_type;
    typedef typename crsmat_type::size_type size_type;

    typename crsmat_mutable_type::entries_type entries;
    typename crsmat_mutable_type::rowmap_type rowmap;
    crsmat_mutable_type::base_type::traits_type::init( entries );
    crsmat_mutable_type::base_type::traits_type::init( rowmap );
    size_type cursor1 = 0;    // current entry index in the first distance index
    size_type cursor2 = 0;    // current entry index in the second distance index
    size_type end1;  // last entry index of the row in the first distance index
    size_type end2;  // last entry index of the row in the second distance index

    assert( dindex1.numRows() == dindex2.numRows() );
    assert( dindex1.numCols() == dindex2.numCols() );

    auto nof_rows = dindex1.numRows();
    auto nof_cols = dindex1.numCols();

    for ( ordinal_type nrow = 0; nrow < nof_rows; ++nrow ) {
      rowmap.push_back( entries.size() );

      end1 = dindex1.rowMap( nrow + 1 );
      end2 = dindex2.rowMap( nrow + 1 );
      while ( cursor1 < end1 ) {
        if ( cursor2 >= end2 ) {
          while ( cursor1 < end1 ) entries.push_back( dindex1.entry( cursor1++ ) );
          break;
        }
        if ( dindex2.entry( cursor2 ) < dindex1.entry( cursor1 ) ) {
          entries.push_back( dindex2.entry( cursor2++ ) );
        } else {
          entries.push_back( dindex1.entry( cursor1 ) );
          if ( dindex2.entry( cursor2 ) == dindex1.entry( cursor1 ) ) ++cursor2;
          ++cursor1;
        }
      }
      while ( cursor2 < end2 ) entries.push_back( dindex2.entry( cursor2++ ) );
    }
    rowmap.push_back( entries.size() );
    assert( cursor1 == dindex1.nnz() );
    assert( cursor2 == dindex2.nnz() );

    return crsmat_mutable_type( nof_cols, std::move( entries ), std::move( rowmap ) );
  }

  template< typename TMutableCRSMatrix, typename TCRSMatrix >
  inline TMutableCRSMatrix
  merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2 )
  {
    return merge_distance_index< TMutableCRSMatrix >(
      dindex1, dindex2, typename crs_matrix::Group< typename TCRSMatrix::spec_type >::type{} );
  }
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_CRS_MATRIX_HPP__ --- */
