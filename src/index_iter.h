/**
 *    @file  index_iter.h
 *   @brief  Custom iterators for SeqAn indexes.
 *
 *  The header file includes custom iterator template classes for SeqAn indexes.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Mar 01, 2017  14:27
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  PSI_INDEX_ITER_H__
#define  PSI_INDEX_ITER_H__

#include <vector>
#include <functional>
#include <algorithm>
#include <type_traits>

#include "index.h"
#include "seed.h"

namespace grem {
  /**
   *  @brief  Top-down index iterator tag with finer traversal in virtual suffix tree.
   *
   *  @tparam  TSpec The type for further specialization of the `TopDownFine` tag type.
   *
   *  The `TSpec` can `seqan::Preorder` or `seqan::ParentLinks<>`. It specifies internal
   *  SeqAn iterator of IndexIter class.
   */
  template< typename TSpec = seqan::Preorder >
    struct TopDownFine;

  /**
   *  @brief  Generic index iterator template class.
   *
   *  @tparam  TIndex The index type.
   *  @tparam  TSpec The type for further specialization of the `IndexIter` type.
   */
  template< typename TIndex, typename TSpec >
    class IndexIter;

  /* Forwards  ------------------------------------------------------------------- */
  template< typename TIndex, typename TSpec >
      bool
    go_down( IndexIter< TIndex, TopDownFine< TSpec > > &iterator,
        typename seqan::Value< TIndex >::Type c);

  template< typename TIndex, typename TSpec >
      bool
    go_down_on_edge( IndexIter< TIndex, TopDownFine< TSpec > > &iterator,
        typename seqan::Value< TIndex >::Type c );

  template< typename TIndex, typename TSpec >
      bool
    is_root( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

  template< typename TIndex, typename TSpec >
      bool
    go_down( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

  template< typename TIndex >
      bool
    go_up( IndexIter< TIndex, TopDownFine< seqan::ParentLinks<> > > &iterator );

  template< typename TIndex, typename TSpec >
      bool
    go_right( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

  template< typename TIndex, typename TSpec >
      typename seqan::Value< TIndex >::Type
    parent_edge_label( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

  template< typename TIndex, typename TSpec >
      typename seqan::Size< TIndex >::Type
    rep_length( const IndexIter< TIndex, TopDownFine< TSpec > >& iterator );
  /* END OF Forwards  ------------------------------------------------------------ */

  /**
   *  @brief  IndexIter specialization for `TopDownFine` iterator tag.
   *
   *  This specialization is a wrapper class for `Iterator< TopDown< TSpec > >::Type`
   *  allowing to traverse the virtual suffix tree finer; i.e. one character at a time
   *  which cannot be done by regular `TopDown` iterator.
   */
  template< typename TIndex, typename TSpec >
    class IndexIter< TIndex, TopDownFine< TSpec > >
    {
      /* ====================  FRIENDSHIP    ======================================= */
      friend bool
        go_down< TIndex, TSpec >( IndexIter< TIndex, TopDownFine< TSpec > > &iterator,
            typename seqan::Value< TIndex >::Type c );

      friend bool
        go_down_on_edge< TIndex, TSpec >(
            IndexIter< TIndex, TopDownFine< TSpec > > &iterator,
            typename seqan::Value< TIndex >::Type c );

      friend bool
        is_root< TIndex, TSpec >( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

      friend bool
        go_down< TIndex, TSpec >( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

      friend bool
        go_up< TIndex >(
            IndexIter< TIndex, TopDownFine< seqan::ParentLinks<> > > &iterator );

      friend bool
        go_right< TIndex, TSpec >( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

      friend typename seqan::Value< TIndex >::Type
        parent_edge_label< TIndex, TSpec >(
            IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

      friend typename seqan::Size< TIndex >::Type
        rep_length< TIndex, TSpec >(
            const IndexIter< TIndex, TopDownFine< TSpec > >& iterator );

      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename seqan::Iterator< TIndex, seqan::TopDown< TSpec > >::Type base_type;
        typedef typename seqan::Container< base_type >::Type container_type;
        /* ====================  LIFECYCLE     ======================================= */
        IndexIter( TIndex &index ) :                             /* constructor */
          iter_(index), boffset(0) { }

        /* ====================  ACCESSORS     ======================================= */

        /**
         *  @brief  getter function for iter_.
         */
        inline const typename seqan::Iterator< TIndex, seqan::TopDown< TSpec > >::Type &
          get_iter_(  ) const
          {
            return iter_;
          }  /* -----  end of method IndexIter::get_iter_  ----- */

        /**
         *  @brief  getter function for boffset.
         */
        inline unsigned int
          get_boffset(  ) const
          {
            return boffset;
          }  /* -----  end of method IndexIter::get_boffset  ----- */

      private:
        /* ====================  DATA MEMBERS  ======================================= */

        /** @brief Internal regular `TopDown< TSpec >` iterator. */
        base_type iter_;
        /**
         *  @brief  Backward offset.
         *
         * If `iterator` is an instance of this class,
         * @code
         * repr(iterator) = repr(iter_)[0:-boffset]
         * @endcode
         * where `iter_` is the internal `TopDown` iterator, and `repr` is the
         * representative string of the given index iterator.
         */
        unsigned int boffset;
    };  /* ----------  end of template class IndexIter  ---------- */

  /**
   *  @brief  Go down in virtual suffix tree by the given character character.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @param[in]      c `iterator` goes down the edge beginning with `c`.
   *  @return `true` if the edge or path to go down exists, otherwise `false`.
   *
   *  A wrapper function for `seqan::goDown` allowing to go down the virtual suffix
   *  tree finer; i.e. one character at a time. A regular `TopDown` iterator cannot
   *  perform the same because the virtual suffix tree is a compressed trie so that
   *  the new representative string of iterator after going down might be extended
   *  by more than one character.
   */
  template< typename TIndex, typename TSpec >
      bool
    go_down( IndexIter< TIndex, TopDownFine< TSpec > > &iterator,
        typename seqan::Value< TIndex >::Type c )
    {
      // :TODO:Thu Mar 02 05:21:\@cartoonist: Handle "N" outside of the function.
      if ( iterator.boffset == 0 ) {          // iterator points a node.
        // Go down using `seqan::goDown` function updating the internal iterator.
        if ( seqan::goDown( iterator.iter_, c ) ) {
          // Set `boffset` such that it points to the first char of the parent edge.
          iterator.boffset = parentEdgeLength( iterator.iter_ ) - 1;
          return true;
        } else {                            // iterator cannot go down further.
          return false;
        }
      } else {                              // iterator points a char on an edge.
        return go_down_on_edge( iterator, c );
      }
    }  /* -----  end of template function go_down  ----- */

  /**
   *  @brief  Go down when iterator points to the given character on an edge.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @param[in]      c `iterator` goes down if next character on the edge is `c`.
   *  @return `true` if the path to go down exists, otherwise `false`.
   *
   *  Looking up the next character on the edge, go down by updating `boffset` if
   *  the next char on the edge is `c`. Internal iterator does not change in this
   *  function.
   */
  template< typename TIndex, typename TSpec >
      bool
    go_down_on_edge( IndexIter< TIndex, TopDownFine< TSpec > > &iterator,
        typename seqan::Value< TIndex >::Type c )
    {
      auto const &parent_edge_str = parentEdgeLabel( iterator.iter_ );
      auto const &parent_edge_len = parentEdgeLength( iterator.iter_ );
      auto const &next_char = parent_edge_str[ parent_edge_len - iterator.boffset ];
      if ( c == next_char ) {
        --iterator.boffset;
        return true;
      } else {
        return false;
      }
    }  /* -----  end of function go_down_on_edge  ----- */

  /**
   *  @brief  Check whether the index iterator points to the root node.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @return `true` if it points to the root node, otherwise `false`.
   *
   *  This function is just a wrapper for `seqan::isRoot` function which is called on
   *  internal index iterator.
   */
  template< typename TIndex, typename TSpec >
      bool
    is_root( IndexIter< TIndex, TopDownFine< TSpec > > &iterator )
    {
      if ( iterator.boffset != 0 ) {      // iterator points an edge.
        return false;
      }
      else {
        return seqan::isRoot ( iterator.iter_ );
      }
    }  /* -----  end of template function is_root  ----- */

  /**
   *  @brief  Go down in virtual suffix tree by one character.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @return `true` if the edge or path to go down exists, otherwise `false`.
   *
   *  Such as `go_down(iterator, c)` but go down the virtual suffix tree in preorder.
   */
  template< typename TIndex, typename TSpec >
      bool
    go_down( IndexIter< TIndex, TopDownFine< TSpec > > &iterator )
    {
      // :TODO:Thu Mar 02 05:21:\@cartoonist: Handle "N" outside of the function.
      if ( iterator.boffset == 0 ) {          // iterator points a node.
        // Go down using `seqan::goDown` function updating the internal iterator.
        if ( seqan::goDown(iterator.iter_) ) {
          // Set `boffset` such that it points to the first char of the parent edge.
          iterator.boffset = parentEdgeLength( iterator.iter_ ) - 1;
          return true;
        } else {                            // iterator cannot go down further.
          return false;
        }
      } else {                              // iterator points a char on an edge.
        --iterator.boffset;
        return true;
      }
    }  /* -----  end of template function go_down  ----- */

  /**
   *  @brief  Go up in virtual suffix tree by one character.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @return `true` when go up is possible, otherwise (root node) `false`.
   *
   *  A wrapper function for `seqan::goUp` allowing to go up the virtual suffix tree
   *  finer; i.e. one character at a time. This function is only specialized for
   *  `TopDownFine< ParentLinks<> >` iterator with history.
   */
  template< typename TIndex >
      bool
    go_up ( IndexIter< TIndex, TopDownFine< seqan::ParentLinks<> > > &iterator )
    {
      if ( is_root (iterator) ) {
        return false;
      }

      if ( iterator.boffset == parentEdgeLength(iterator.iter_) - 1 ) {  // first char on edge?
        // Go up using `seqan::goUp` function updating the internal iterator.
        if ( seqan::goUp ( iterator.iter_ ) ) {
          // Set `boffset` such that it points to the last char of the parent edge.
          iterator.boffset = 0;
          return true;
        }
        else {
          assert( false );
          // SHOULD NOT BE REACHED
          // Since any non-root node should be able to go up.
          return false;
        }
      }
      else {
        ++iterator.boffset;
        return true;
      }
    }  /* -----  end of template function go_up  ----- */

  /**
   *  @brief  Go right (a sibling node) in the virtual suffix tree.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @return `true` when go right is possible, otherwise `false`.
   *
   *  This function is a wrapper for `seqan::goRight`.
   */
  template< typename TIndex, typename TSpec >
      bool
    go_right( IndexIter< TIndex, TopDownFine< TSpec > > &iterator )
    {
      if ( iterator.boffset == parentEdgeLength( iterator.iter_ ) - 1 ) {  // iterator points a node?
        // Go right using `seqan::goRight` function updating the internal iterator.
        if ( seqan::goRight ( iterator.iter_ ) ) {
          // Set `boffset` such that it points to the first char of the parent edge.
          iterator.boffset = parentEdgeLength ( iterator.iter_ ) - 1;
          return true;
        }
        else {
          return false;
        }
      }
      else {
        return false;
      }
    }  /* -----  end of template function go_right  ----- */

  /**
   *  @brief  Get the parent edge label.
   *
   *  @param  iterator The iterator of virtual suffix tree.
   *  @return The label of the incoming (parent) edge which is a one character in
   *  `TopDownFine` iterators.
   *
   *  This function is a wrapper for `seqan::parentEdgeLabel` method.
   */
  template< typename TIndex, typename TSpec >
      typename seqan::Value< TIndex >::Type
    parent_edge_label( IndexIter< TIndex, TopDownFine< TSpec > > &iterator )
    {
      auto const &parent_edge_str = parentEdgeLabel ( iterator.iter_ );
      auto const &parent_edge_len = parentEdgeLength ( iterator.iter_ );
      return parent_edge_str [ parent_edge_len - iterator.boffset - 1 ];
    }  /* -----  end of template function parent_edge_label  ----- */

  /**
   *  @brief  Get length of representative string.
   *
   *  @param  iterator The iterator of virtual suffix tree.
   *  @return The length of the representative string of the given iterator.
   *
   *  This is a wrapper function for `seqan::repLength` method.
   */
  template< typename TIndex, typename TSpec >
      typename seqan::Size< TIndex >::Type
    rep_length( const IndexIter< TIndex, TopDownFine< TSpec > >& iterator )
    {
      return repLength( iterator.iter_ ) - iterator.boffset;
    }  /* -----  end of template function rep_length  ----- */
  /* Typedefs  ------------------------------------------------------------------- */

  /**
   *  @brief  Shorter alias for `Iterator<>::Type` syntax.
   */
  template< typename TIndex, typename TSpec >
    using TIndexIter = typename seqan::Iterator< TIndex, TSpec >::Type;

  template< typename TIter >
    using TIterIndex = typename seqan::Container< TIter >::Type;

  template< typename TIter >
    using TIterText = typename seqan::Fibre< TIterIndex< TIter >, seqan::FibreText >::Type;

  template< typename TIter >
    using TIterRawText = typename seqan::Value< TIterText< TIter > >::Type;

  /* END OF Typedefs  ------------------------------------------------------------ */

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Iterator type specialization for `grem::TopDownFine< TSpec >` tag.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template< typename TIndex, typename TSpec >
    class Iterator< TIndex, grem::TopDownFine< TSpec > >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef grem::IndexIter< TIndex, grem::TopDownFine< TSpec > > Type;
    };  /* ----------  end of template class Iterator  ---------- */

  template< typename TIndex, typename TSpec >
    class Container< grem::IndexIter< TIndex, TSpec > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename grem::IndexIter< TIndex, TSpec >::container_type Type;
    };  /* ----------  end of template class Container  ---------- */

}  /* -----  end of namespace seqan  ----- */

namespace grem {

  /* Typedefs  ------------------------------------------------------------------- */
  template< typename TIndex, typename TSpec = seqan::Preorder >
    using TFineIndexIter =
      typename seqan::Iterator< TIndex, TopDownFine< TSpec > >::Type;
  /* END OF Typedefs  ------------------------------------------------------------ */

  /* Index interator interface functions  ---------------------------------------- */

  template< typename TIndex, typename TIterSpec >
      inline bool
    go_right_stree( seqan::Iter< TIndex, TIterSpec >& iter )
    {
      return goRight( iter );
    }

  template< typename TText, typename TIterSpec >
      inline bool
    go_right_stree( seqan::Iter< seqan::Index< TText, grem::CBiFMIndex >, TIterSpec >& iter )
    {
      return goRight( iter, seqan::Rev() );
    }

  template< typename TIndex, typename TIterSpec >
      inline bool
    go_down_stree( seqan::Iter< TIndex, TIterSpec >& iter )
    {
      return goDown( iter );
    }

  template< typename TText, typename TIterSpec >
      inline bool
    go_down_stree( seqan::Iter< seqan::Index< TText, grem::CBiFMIndex >, TIterSpec >& iter )
    {
      return goDown( iter, seqan::Rev() );
    }

  template< typename TIndex, typename TIterSpec, typename TPattern >
      inline bool
    go_down_stree( seqan::Iter< TIndex, TIterSpec >& iter, TPattern&& p )
    {
      return goDown( iter, p );
    }

  template< typename TText, typename TIterSpec, typename TPattern >
      inline bool
    go_down_stree( seqan::Iter< seqan::Index< TText, grem::CBiFMIndex >, TIterSpec >& iter,
        TPattern&& p )
    {
      return goDown( iter, p, seqan::Rev() );
    }

  template< typename TIndex, typename TIterSpec >
      inline auto
    parent_edge_char_stree( const seqan::Iter< TIndex, TIterSpec >& iter )
    {
      return parentEdgeLabel( iter )[0];
    }

  template< typename TText, typename TIterSpec >
      inline auto
    parent_edge_char_stree( const seqan::Iter< seqan::Index< TText, grem::CBiFMIndex >, TIterSpec >& iter )
    {
      return parentEdgeLabel( iter, seqan::Rev() );
    }

  template< typename TIndex, typename TIterSpec >
      inline unsigned int
    parent_edge_len_stree( const seqan::Iter< TIndex, TIterSpec >& iter )
    {
      return length( parentEdgeLabel( iter ) );
    }

  template< typename TText, typename TIterSpec >
      inline unsigned int
    parent_edge_len_stree( const seqan::Iter< seqan::Index< TText, grem::CBiFMIndex >, TIterSpec >& iter )
    {
      if ( isRoot( iter ) ) return 0;
      return length( parentEdgeLabel( iter, seqan::Rev() ) );
    }

  template< typename TIndex, typename TIterSpec >
      inline auto
    get_occurrences_stree( const seqan::Iter< TIndex, TIterSpec >& iter )
    {
      return getOccurrences( iter );
    }

  template< typename TText, typename TIterSpec >
      inline auto
    get_occurrences_stree( const seqan::Iter< seqan::Index< TText, grem::CBiFMIndex >, TIterSpec >& iter )
    {
      return getOccurrences( iter, seqan::Rev() );
    }

  /* END OF Index interator interface functions  --------------------------------- */

  template< typename TIter >
      inline bool
    next_kmer( TIter& itr, unsigned int& cp_len, unsigned int k )
    {
      cp_len = repLength( itr );
      do {
        if ( repLength( itr ) >= k ||
            !go_down_stree( itr ) ||
            parent_edge_char_stree( itr ) == 'N' ) {
          if ( !go_right_stree( itr ) ) {
            while ( goUp( itr ) && !go_right_stree( itr ) );
          }
          cp_len = std::min( cp_len,
              static_cast< unsigned int >( repLength( itr ) - parent_edge_len_stree( itr ) ) );
        }

        if ( repLength( itr ) >= k ) {
          return true;
        }
      } while ( !isRoot( itr ) );        /* -----  end do-while  ----- */

      return false;
    }

  template< typename TIter >
      inline typename seqan::Size< TIter >::Type
    upto_prefix( TIter& itr, const unsigned int& cp_len )
    {
      while( repLength( itr ) > cp_len ) goUp( itr );
      return repLength( itr );
    }

  template< typename TOccurrence1, typename TOccurrence2, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    _add_seed( TOccurrence1 oc1, TOccurrence2 oc2,
        const TRecords1* rec1, const TRecords2* rec2, TCallback callback )
    {
      Seed<> hit;
      hit.node_id = position_to_id( *rec1, oc1 );
      hit.node_offset = position_to_offset( *rec1, oc1 );
      hit.read_id = position_to_id( *rec2, oc2.i1 );
      hit.read_offset = position_to_offset( *rec2, oc2 );

      callback( hit );
    }

  template< typename TOccurrence >
      inline TOccurrence
    _map_occurrences( TOccurrence const& oc, unsigned int k, Forward )
    {
      return TOccurrence( oc );
    }

  template< typename TOccurrence >
      inline TOccurrence
    _map_occurrences( TOccurrence const& oc, unsigned int k, Reversed )
    {
      return TOccurrence( oc.i1, oc.i2 + k - 1 );  /**< @brief End position of the occurrence. */
    }

  template< typename TIter1, typename TIter2, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    _add_occurrences( TIter1& itr1, TIter2& itr2, const TRecords1* rec1,
        const TRecords2* rec2, unsigned int k, TCallback callback )
    {
      typedef typename Direction< TRecords1 >::Type TPathDir;

      using seqan::length;

      const auto& occurrences1 = get_occurrences_stree( itr1 );
      const auto& occurrences2 = get_occurrences_stree( itr2 );

      for ( unsigned int i = 0; i < length( occurrences1 ); ++i ) {
        for ( unsigned int j = 0; j < length( occurrences2 ); ++j ) {
          auto oc = _map_occurrences( occurrences1[i], k, TPathDir() );
          _add_seed( oc, occurrences2[j], rec1, rec2, callback );
        }
      }
    }

  template< typename TIter1, typename TIter2, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    _kmer_exact_match_impl( TIter1& fst_itr, TIter2& snd_itr, const TRecords1* rec1,
        const TRecords2* rec2, unsigned int k, bool swapped, TCallback callback )
    {
      unsigned int cp_len;
      while ( next_kmer( fst_itr, cp_len, k ) ) {
        auto&& s = upto_prefix( snd_itr, cp_len );

        if ( go_down_stree( snd_itr, infix( representative( fst_itr ), s, k ) ) ) {
          if ( !swapped ) _add_occurrences( fst_itr, snd_itr, rec1, rec2, k, callback );
          else _add_occurrences( snd_itr, fst_itr, rec1, rec2, k, callback );
        }
      }
    }

  template< typename TText1, typename TText2, typename TIndexSpec1, typename TIndexSpec2, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    kmer_exact_matches( seqan::Index< TText1, TIndexSpec1 >& fst,
        seqan::Index< TText2, TIndexSpec2 >& snd,
        const TRecords1* rec1,
        const TRecords2* rec2,
        unsigned int k,
        TCallback callback )
    {
      static_assert( ( is_fmindex< TIndexSpec1 >::value &&
            std::is_same< typename Direction< TRecords1 >::Type, Reversed >::value ) ||
          ( !is_fmindex< TIndexSpec1 >::value &&
            std::is_same< typename Direction< TRecords1 >::Type, Forward >::value ),
          "The paths direction and the path index used are not compatible." );

      if ( k == 0 ) return;

      auto fst_len = length( indexRawText( fst ) );
      auto snd_len = length( indexRawText( snd ) );

      typedef seqan::Index< TText1, TIndexSpec1 > TIndex1;
      typedef seqan::Index< TText2, TIndexSpec2 > TIndex2;
      typedef seqan::TopDown< seqan::ParentLinks<> > TIterSpec;
      TIndexIter< TIndex1, TIterSpec > fst_itr( fst );
      TIndexIter< TIndex2, TIterSpec > snd_itr( snd );

      if ( fst_len <= snd_len ) {
        _kmer_exact_match_impl( fst_itr, snd_itr, rec1, rec2, k, false, callback );
      }
      else {
        _kmer_exact_match_impl( snd_itr, fst_itr, rec1, rec2, k, true, callback );
      }
    }

  template< typename TIndex, typename TIterSpec >
      inline typename seqan::Size< IndexIter< TIndex, TopDownFine< TIterSpec > > >::Type
    upto_prefix( IndexIter< TIndex, TopDownFine< TIterSpec > >& itr,
        const unsigned int& cp_len )
    {
      while( rep_length( itr ) > cp_len ) go_up( itr );
      return rep_length( itr );
    }

  template< typename TIndex1, typename TIndex2, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    kmer_exact_matches( IndexIter< TIndex1, TopDownFine< seqan::ParentLinks<> > >& fst_itr,
        IndexIter< TIndex2, TopDownFine< seqan::ParentLinks<> > >& snd_itr,
        const TRecords1* rec1,
        const TRecords2* rec2,
        unsigned int k,
        TCallback callback )
    {
      static_assert( ( is_fmindex< typename seqan::Spec< TIndex1 >::Type >::value &&
            std::is_same< typename Direction< TRecords1 >::Type, Reversed >::value ) ||
          ( !is_fmindex< typename seqan::Spec< TIndex1 >::Type >::value &&
            std::is_same< typename Direction< TRecords1 >::Type, Forward >::value ),
          "The paths direction and the path index used are not compatible." );

      if ( k == 0 ) return;

      seqan::DnaString seed;                   // seed = A..(k)..A
      for ( unsigned int i = 0; i < k; ++i ) appendValue( seed, 'A' );

      int s = 0;
      do {
        unsigned int plen;
        bool fst_agreed = true;
        bool snd_agreed = true;

        upto_prefix( fst_itr, s );
        upto_prefix( snd_itr, s );
        for ( plen = s; fst_agreed && snd_agreed && plen < k; ++plen ) {
          fst_agreed = go_down( fst_itr, seed[plen] );
          snd_agreed = go_down( snd_itr, seed[plen] );
        }
        if ( fst_agreed && snd_agreed ) {
          _add_occurrences( fst_itr.get_iter_(), snd_itr.get_iter_(), rec1, rec2, k, callback );
          plen = 0;
        }
        s = increment_kmer( seed, plen );
      } while ( s >= 0 );
    }

  template< typename TIndex, typename TRecords1, typename TRecordsIter, typename TCallback >
      inline void
    kmer_exact_matches( TIndex& paths_index, const TRecords1* pathset,
        TRecordsIter& seeds_itr, TCallback callback )
    {
      static_assert( std::is_same< typename Direction< TRecords1 >::Type, Forward >::value,
          "The paths should be forward sequences." );

      seqan::Finder< TIndex > paths_finder( paths_index );
      while ( !at_end( seeds_itr ) ) {
        while ( find( paths_finder, *seeds_itr ) ) {
          _add_seed( beginPosition( paths_finder ), get_position( seeds_itr ), pathset,
              seeds_itr.get_records_ptr(), callback );
        }
        ++seeds_itr;
        clear( paths_finder );
      }
    }

  template< typename TIndex, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    all_exact_matches( TIndex& paths_index, const TRecords1* pathset,
        const TRecords2* reads, TCallback callback )
    {
      typedef typename TRecords2::TSize size_type;
      typedef typename TRecords2::TStringSetPosition pos_type;

      static_assert( std::is_same< typename Direction< TRecords1 >::Type, Forward >::value,
          "The paths should be forward sequences." );

      seqan::Finder< TIndex > paths_finder( paths_index );
      for ( size_type i = 0; i < length( *reads ); ++i ) {
        while ( find( paths_finder, reads->str[i] ) ) {
          _add_seed( beginPosition( paths_finder ), pos_type( { i, 0 } ), pathset, reads, callback );
        }
        clear( paths_finder );
      }
    }
}  /* -----  end of namespace grem  ----- */

#endif  /* --- #ifndef PSI_INDEX_ITER_H__ --- */
