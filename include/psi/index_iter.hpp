/**
 *    @file  index_iter.hpp
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

#ifndef  PSI_INDEX_ITER_HPP__
#define  PSI_INDEX_ITER_HPP__

#include <vector>
#include <functional>
#include <algorithm>
#include <type_traits>

#include "index.hpp"
#include "seed.hpp"

namespace psi {
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
      void
    go_root( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

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

      friend void
        go_root< TIndex, TSpec >( IndexIter< TIndex, TopDownFine< TSpec > > &iterator );

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
   *  @brief  IndexIter TopDownFine class specialization for psi::FMIndex.
   *
   *  The index iterator supports trie traversal by default.
   */
  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
    class IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >
    : public seqan::Iter< seqan::Index< TText, psi::FMIndex< TWT, TDens, TInvDens > >, seqan::TopDown< TSpec > >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TSpec spec_type;
        typedef seqan::Index< TText, psi::FMIndex< TWT, TDens, TInvDens > > index_type;
        typedef index_type container_type;
        typedef seqan::Iter< index_type, seqan::TopDown< TSpec > > base_type;
        /* ====================  LIFECYCLE     ======================================= */
        IndexIter( index_type const* index_ptr )
          : base_type( index_ptr ) { }
        IndexIter( index_type const& index )
          : base_type( index ) { }

          inline base_type const&
        get_iter_(  ) const
        {
          return *this;
        }
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

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    go_down( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >& iterator,
        typename seqan::Value< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > > >::Type c )
    {
      return goDown( iterator, c );
    }

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

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    is_root( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >& iterator )
    {
      return isRoot( iterator );
    }

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

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline bool
    go_up( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< seqan::ParentLinks<> > >& iterator )
    {
      return goUp( iterator );
    }

  template< typename TIndex, typename TSpec >
      inline void
    go_root( IndexIter< TIndex, TopDownFine< TSpec > > &iterator )
    {
      seqan::goRoot( iterator.iter_ );
      iterator.boffset = 0;
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline void
    go_root( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >& iterator )
    {
      goRoot( iterator );
    }

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

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline bool
    go_right( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >& iterator )
    {
      return goRight( iterator );
    }

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

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename seqan::Value< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > > >::Type
    parent_edge_label( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >& iterator )
    {
      return parentEdgeLabel( iterator );
    }

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

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens, typename TSpec >
      inline typename seqan::Size< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > > >::Type
    rep_length( IndexIter< seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >, TopDownFine< TSpec > >& iterator )
    {
      return repLength( iterator );
    }

  template< typename TIndex, typename TSpec >
      typename seqan::Size< TIndex >::Type
    count_occurrences( const IndexIter< TIndex, TopDownFine< TSpec > >& iterator )
    {
      return countOccurrences( iterator.get_iter_() );
    }

  template< typename TIndex, typename TSpec >
      inline auto
    get_occurrences( const IndexIter< TIndex, TopDownFine< TSpec > >& iterator )
    {
      return getOccurrences( iterator.get_iter_() );
    }
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

}  /* --- end of namespace psi --- */

namespace seqan {
  /**
   *  @brief  Iterator type specialization for `psi::TopDownFine< TSpec >` tag.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template< typename TIndex, typename TSpec >
    class Iterator< TIndex, psi::TopDownFine< TSpec > >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef psi::IndexIter< TIndex, psi::TopDownFine< TSpec > > Type;
    };  /* ----------  end of template class Iterator  ---------- */

  template< typename TIndex, typename TSpec >
    class Container< psi::IndexIter< TIndex, TSpec > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename psi::IndexIter< TIndex, TSpec >::container_type Type;
    };  /* ----------  end of template class Container  ---------- */

}  /* -----  end of namespace seqan  ----- */

namespace psi {

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
    go_right_stree( seqan::Iter< seqan::Index< TText, psi::CBiFMIndex >, TIterSpec >& iter )
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
    go_down_stree( seqan::Iter< seqan::Index< TText, psi::CBiFMIndex >, TIterSpec >& iter )
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
    go_down_stree( seqan::Iter< seqan::Index< TText, psi::CBiFMIndex >, TIterSpec >& iter,
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
    parent_edge_char_stree( const seqan::Iter< seqan::Index< TText, psi::CBiFMIndex >, TIterSpec >& iter )
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
    parent_edge_len_stree( const seqan::Iter< seqan::Index< TText, psi::CBiFMIndex >, TIterSpec >& iter )
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
    get_occurrences_stree( const seqan::Iter< seqan::Index< TText, psi::CBiFMIndex >, TIterSpec >& iter )
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
    upto_prefix( TIter& itr, unsigned int cp_len )
    {
      while( repLength( itr ) > cp_len ) goUp( itr );
      return repLength( itr );
    }

  template< typename TOccurrence1, typename TOccurrence2, typename TRecords1, typename TRecords2, typename TCallback >
      inline void
    _add_seed( TOccurrence1 oc1, TOccurrence2 oc2,
        const TRecords1* rec1, const TRecords2* rec2, unsigned int len, unsigned int gocc,
        TCallback callback )
    {
      Seed<> hit;
      hit.node_id = position_to_id( *rec1, oc1 );
      hit.node_offset = position_to_offset( *rec1, oc1 );
      hit.read_id = position_to_id( *rec2, oc2.i1 );
      hit.read_offset = position_to_offset( *rec2, oc2 );
      hit.match_len = len;
      hit.gocc = gocc;

      callback( hit );
    }

  template< typename TOccurrence >
      inline TOccurrence
    _map_occurrences( TOccurrence const& oc, unsigned int k, Forward )
    {
      return TOccurrence( oc );
    }

  /**
   *  @brief  Adjust occurrence positions based on the text direction of the underlying index.
   *
   *  @param oc  occurrence position of a pattern.
   *  @param len the length of the occurrence.
   *  @param tag Reverse direction tag.
   *
   *  The occurrence position of a pattern in the path set would be a pair `(i, o)`
   *  where the pattern occurs in path with id `i` at offset `o`. If the path sequences
   *  are reversed, the position points at the end of pattern occurrence in reversed
   *  sequence; e.g.:
   *
   *  The pattern 'ttc' is found in the reversed sequence:
   *         0123 456 7890123
   *         acga ctt taggtcc
   *
   *  The reported occurrence offset would be 4 while the actual offset would be 7; i.e.
   *
   *     offset_{fwd,start} = | sequence | - offset_{rev,end} + 1
   *     offset_{fwd,start} = 14 - 6 - 1 = 7
   *
   *  where
   *
   *     offset_{rev,end} = offset_{rev,start} + | pattern | - 1
   *     offset_{rev,end} = 4 + 3 - 1 = 6
   *
   *  The first part are calculated in `position_to_offset` interface functions since
   *  the length of the sequence is only known there. The second part is done in this
   *  function.
   *
   *  NOTE: If the text direction is forward, nothing needs to be done.
   */
  template< typename TOccurrence >
      inline TOccurrence
    _map_occurrences( TOccurrence const& oc, unsigned int len, Reversed /* tag */ )
    {
      return TOccurrence( oc.i1, oc.i2 + len - 1 );  /**< @brief End position of the occurrence. */
    }

  /**
   * NOTE: itr1 should be index iterator of the genome.
   */
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
          _add_seed( oc, occurrences2[j], rec1, rec2, k, length( occurrences1 ), callback );
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
      inline void
    upto_prefix( IndexIter< TIndex, TopDownFine< TIterSpec > >& itr,
        unsigned int cp_len )
    {
      auto rlen = rep_length( itr );
      if ( rlen <= cp_len ) return;
      while ( rlen-- != cp_len ) go_up( itr );
    }

  template< typename TIndex1, typename TIndex2, typename TRecords1, typename TRecords2,
            typename TCallback, typename TStats = std::function< void( std::size_t, bool ) > >
      inline void
    kmer_exact_matches( IndexIter< TIndex1, TopDownFine< seqan::ParentLinks<> > >& fst_itr,
        IndexIter< TIndex2, TopDownFine< seqan::ParentLinks<> > >& snd_itr,
        const TRecords1* rec1,
        const TRecords2* rec2,
        unsigned int k,
        TCallback callback,
        unsigned int gocc_threshold = 0,
        TStats collect_stats=[]( std::size_t, bool )->void{} )
    {
      static_assert( ( is_fmindex< typename seqan::Spec< TIndex1 >::Type >::value &&
            std::is_same< typename Direction< TRecords1 >::Type, Reversed >::value ) ||
          ( !is_fmindex< typename seqan::Spec< TIndex1 >::Type >::value &&
            std::is_same< typename Direction< TRecords1 >::Type, Forward >::value ),
          "The paths direction and the path index used are not compatible." );

      if ( k == 0 ) return;
      if ( gocc_threshold == 0 ) {
        gocc_threshold = std::numeric_limits< decltype( gocc_threshold ) >::max();
      }

      seqan::DnaString seed;                   // seed = A..(k)..A
      for ( unsigned int i = 0; i < k; ++i ) appendValue( seed, 'A' );

      unsigned int plen = 0;
      do {
        upto_prefix( fst_itr, plen );
        upto_prefix( snd_itr, plen );
        for ( ; plen < k; ++plen ) {
          if ( !go_down( fst_itr, seed[plen] ) ) break;
          if ( !go_down( snd_itr, seed[plen] ) ) break;
        }
        if ( plen == k ) {
          auto count = count_occurrences( fst_itr );
          if ( count <= gocc_threshold ) {
            collect_stats( count, false );
            _add_occurrences( fst_itr.get_iter_(), snd_itr.get_iter_(), rec1, rec2, k, callback );
          } else collect_stats( count, true );
          --plen;
        }
        plen = increment_kmer( seed, plen, true );
      } while ( plen + 1 > 0 );
    }

  template< typename TString, typename TIndex, typename TSpec, typename TRecords, typename TCallback >
      inline void
    find_mems( TString const& pattern,
               IndexIter< TIndex, TopDownFine< TSpec > >& idx_itr,
               const TRecords* pathset,
               unsigned int minlen,
               unsigned int context,
               TCallback callback,
               unsigned int gocc_threshold = 0,
               bool find_all=true )
    {
      typedef typename Direction< TRecords >::Type TPathDir;

      if ( gocc_threshold == 0 ) {
        gocc_threshold = std::numeric_limits< decltype( gocc_threshold ) >::max();
      }

      unsigned int start = 0;
      unsigned int plen = 0;
      bool has_hit = false;
      while( start + plen < pattern.size() ) {
        if ( plen >= minlen && count_occurrences( idx_itr ) <= gocc_threshold ) {
          using seqan::length;

          has_hit = true;
          auto const& occs = get_occurrences( idx_itr );
          Seed<> hit;
          for ( unsigned int i = 0; i < length( occs ); ++i ) {
            auto oc = _map_occurrences( occs[i], plen, TPathDir() );
            hit.node_id = position_to_id( *pathset, oc );
            hit.node_offset = position_to_offset( *pathset, oc );
            hit.read_offset = start;
            hit.match_len = plen;
            hit.gocc = length( occs );
            callback( hit );
          }
          if ( !find_all ) break;
        }
        if ( has_hit /*|| plen > context*/ ||
             pattern[ start + plen ] == 'N' ||
             !go_down( idx_itr, pattern[ start + plen ] ) ) {
          go_root( idx_itr );
          start = start + plen + 1;
          plen = 0;
          has_hit = false;
          continue;
        }
        ++plen;
      }
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
          // TODO: set `gocc` value for the hit
          _add_seed( beginPosition( paths_finder ), get_position( seeds_itr ), pathset,
              seeds_itr.get_records_ptr(), length( paths_finder ), 0, callback );
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
          // TODO: set `gocc` value for the hit
          _add_seed( beginPosition( paths_finder ), pos_type( { i, 0 } ), pathset, reads,
              length( paths_finder ), 0, callback );
        }
        clear( paths_finder );
      }
    }
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_INDEX_ITER_HPP__ --- */
