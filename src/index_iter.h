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

#ifndef  INDEX_ITERATOR_H__
#define  INDEX_ITERATOR_H__

#include <seqan/index.h>
#include <seqan/seeds.h>

namespace grem {
  /**
   *  @brief  Top-down index iterator tag with finer traversal in virtual suffix tree.
   *
   *  @tparam  TSpec The type for further specialization of the `TopDownFine` tag type.
   *
   *  The `TSpec` can `seqan::Preorder` or `seqan::ParentLinks<>`. It specifies internal
   *  SeqAn iterator of IndexIter class.
   */
  template < typename TSpec = seqan::Preorder >
    struct TopDownFine;

  /**
   *  @brief  Generic index iterator template class.
   *
   *  @tparam  TIndex The index type.
   *  @tparam  TSpec The type for further specialization of the `IndexIter` type.
   */
  template < typename TIndex, typename TSpec >
    class IndexIter;

  /* Forwards  ------------------------------------------------------------------- */
  template < typename TIndex, typename TSpec >
    bool go_down (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator,
        typename seqan::Value < TIndex >::Type c);

  template < typename TIndex, typename TSpec >
    bool go_down_on_edge (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator,
        typename seqan::Value < TIndex >::Type c);

  template < typename TIndex, typename TSpec >
    bool is_root (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

  template < typename TIndex, typename TSpec >
    bool go_down (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

  template < typename TIndex >
    bool go_up (
        IndexIter < TIndex, TopDownFine < seqan::ParentLinks<> > > &iterator);

  template < typename TIndex, typename TSpec >
    bool go_right (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

  template < typename TIndex, typename TSpec >
    typename seqan::Value < TIndex >::Type parent_edge_label (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator);
  /* END OF Forwards  ------------------------------------------------------------ */

  /**
   *  @brief  IndexIter specialization for `TopDownFine` iterator tag.
   *
   *  This specialization is a wrapper class for `Iterator < TopDown < TSpec > >::Type`
   *  allowing to traverse the virtual suffix tree finer; i.e. one character at a time
   *  which cannot be done by regular `TopDown` iterator.
   */
  template < typename TIndex, typename TSpec >
    class IndexIter < TIndex, TopDownFine < TSpec > >
    {
      /* ====================  FRIENDSHIP    ======================================= */
      friend bool
        go_down < TIndex, TSpec > (
            IndexIter < TIndex, TopDownFine < TSpec > > &iterator,
            typename seqan::Value < TIndex >::Type c);

      friend bool
        go_down_on_edge < TIndex, TSpec > (
            IndexIter < TIndex, TopDownFine < TSpec > > &iterator,
            typename seqan::Value < TIndex >::Type c);

      friend bool
        is_root < TIndex, TSpec > (
            IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

      friend bool
        go_down < TIndex, TSpec > (
            IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

      friend bool
        go_up < TIndex > (
            IndexIter < TIndex, TopDownFine < seqan::ParentLinks<> > > &iterator);

      friend bool
        go_right < TIndex, TSpec > (
            IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

      friend typename seqan::Value < TIndex >::Type
        parent_edge_label < TIndex, TSpec > (
            IndexIter < TIndex, TopDownFine < TSpec > > &iterator);

      public:
        /* ====================  LIFECYCLE     ======================================= */
        IndexIter ( TIndex &index ) :                             /* constructor */
          iter_(index), boffset(0) { }

        /* ====================  ACCESSORS     ======================================= */

        /**
         *  @brief  getter function for iter_.
         */
        inline const typename seqan::Iterator < TIndex, seqan::TopDown< TSpec > >::Type &
          get_iter_ (  ) const
          {
            return iter_;
          }  /* -----  end of method IndexIter::get_iter_  ----- */

        /**
         *  @brief  getter function for boffset.
         */
        inline unsigned int
          get_boffset (  ) const
          {
            return boffset;
          }  /* -----  end of method IndexIter::get_boffset  ----- */

      private:
        /* ====================  DATA MEMBERS  ======================================= */

        /** @brief Internal regular `TopDown< TSpec >` iterator. */
        typename seqan::Iterator < TIndex, seqan::TopDown< TSpec > >::Type iter_;
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
  template < typename TIndex, typename TSpec >
    bool go_down (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator,
        typename seqan::Value < TIndex >::Type c)
  {
    // :TODO:Thu Mar 02 05:21:\@cartoonist: Handle "N" outside of the function.
    if (iterator.boffset == 0) {          // iterator points a node.
      // Go down using `seqan::goDown` function updating the internal iterator.
      if (seqan::goDown(iterator.iter_, c)) {
        // Set `boffset` such that it points to the first char of the parent edge.
        iterator.boffset = parentEdgeLength(iterator.iter_) - 1;
        return true;
      } else {                            // iterator cannot go down further.
        return false;
      }
    } else {                              // iterator points a char on an edge.
      return go_down_on_edge(iterator, c);
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
  template < typename TIndex, typename TSpec >
    bool go_down_on_edge (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator,
        typename seqan::Value < TIndex >::Type c)
  {
    auto const &parent_edge_str = parentEdgeLabel(iterator.iter_);
    auto const &parent_edge_len = parentEdgeLength(iterator.iter_);
    auto const &next_char = parent_edge_str [ parent_edge_len - iterator.boffset ];
      if (c == next_char) {
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
  template < typename TIndex, typename TSpec >
    bool is_root (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator)
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
  template < typename TIndex, typename TSpec >
    bool go_down (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator)
  {
    // :TODO:Thu Mar 02 05:21:\@cartoonist: Handle "N" outside of the function.
    if (iterator.boffset == 0) {          // iterator points a node.
      // Go down using `seqan::goDown` function updating the internal iterator.
      if (seqan::goDown(iterator.iter_)) {
        // Set `boffset` such that it points to the first char of the parent edge.
        iterator.boffset = parentEdgeLength(iterator.iter_) - 1;
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
   *  `TopDownFine < ParentLinks<> >` iterator with history.
   */
  template < typename TIndex >
    bool go_up (
        IndexIter < TIndex, TopDownFine < seqan::ParentLinks<> > > &iterator)
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
  template < typename TIndex, typename TSpec >
    bool go_right (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator)
    {
      if ( iterator.boffset == parentEdgeLength(iterator.iter_) - 1 ) {  // iterator points a node?
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
  template < typename TIndex, typename TSpec >
    typename seqan::Value < TIndex >::Type parent_edge_label (
        IndexIter < TIndex, TopDownFine < TSpec > > &iterator)
    {
      auto const &parent_edge_str = parentEdgeLabel ( iterator.iter_ );
      auto const &parent_edge_len = parentEdgeLength ( iterator.iter_ );
      return parent_edge_str [ parent_edge_len - iterator.boffset - 1 ];
    }  /* -----  end of template function parent_edge_label  ----- */

  /* Typedefs  ------------------------------------------------------------------- */

  /**
   *  @brief  Shorter alias for `Iterator<>::Type` syntax.
   */
  template< typename TIndex, typename TSpec >
    using TIndexIterator = typename seqan::Iterator< TIndex, TSpec >::Type;

  /* END OF Typedefs  ------------------------------------------------------------ */

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Iterator type specialization for `grem::TopDownFine< TSpec >` tag.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template < typename TIndex, typename TSpec >
    class Iterator < TIndex, grem::TopDownFine < TSpec > >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef grem::IndexIter < TIndex, grem::TopDownFine < TSpec > > Type;
    };  /* ----------  end of template class Iterator  ---------- */

}  /* -----  end of namespace seqan  ----- */

namespace grem {

  /* Typedefs  ------------------------------------------------------------------- */
  template < typename TIndex, typename TSpec = seqan::Preorder >
    using TFineIterator =
      typename seqan::Iterator < TIndex, TopDownFine < TSpec > >::Type;
  /* END OF Typedefs  ------------------------------------------------------------ */

  /* Top-down fine interator meta-function declarations  ------------------------- */

  /**
   *  @brief  Find k-mer exact matches between two texts.
   *
   *  @param[out]  seeds The list of k-mer exact matches found as a set of seeds.
   *  @param[in,out]  first First text's top-down fine index iterator.
   *  @param[in,out]  second Second text's top-down fine index iterator.
   *  @param[in]  k The length of k-mers; i.e. `k`.
   *
   *  This function finds the exact k-mer matches in the given two texts by using top-
   *  down fine index iterator. It traverses the index iterator in DFS approach, and
   *  reports representative sequence of nodes in level k which are supported by both
   *  texts' indexes.
   */
  template < typename TIndex1, typename TIndex2 >
    void kmer_exact_matches ( seqan::SeedSet < seqan::Seed < seqan::Simple > > &seeds,
        typename seqan::Iterator < TIndex1, TopDownFine< seqan::ParentLinks<> > >::Type &first,
        typename seqan::Iterator < TIndex2, TopDownFine< seqan::ParentLinks<> > >::Type &second,
        unsigned int k)
    {
      if ( k == 0 ) return;

      unsigned int level = 0;
      bool second_agrees = true;

      do {
        if ( level == k ) {
          typedef typename seqan::SAValue< TIndex1 >::Type TSAValue1;
          typedef typename seqan::SAValue< TIndex2 >::Type TSAValue2;
          seqan::String<TSAValue1> saPositions1 = getOccurrences( first.get_iter_() );
          seqan::String<TSAValue2> saPositions2 = getOccurrences( second.get_iter_() );
          for (unsigned i = 0; i < length(saPositions1); ++i) {
            for (unsigned j = 0; j < length(saPositions2); ++j) {
              // :TODO:Wed Mar 08 10:01:\@cartoonist: typdef SimpleSeed in seed.h?
              seqan::Seed < seqan::Simple > hit;
              seqan::setBeginPositionH ( hit, saPositions1[i].i1 );
              seqan::setEndPositionH ( hit, saPositions1[i].i2 );
              seqan::setBeginPositionV ( hit, saPositions2[j].i1 );
              seqan::setEndPositionV ( hit, saPositions2[j].i2 );

              seqan::addSeed (seeds, std::move(hit), seqan::Single());
            }
          }
        }

        bool right = false;
        bool down = false;
        if ( ( level == k || !second_agrees || !(down = go_down ( first )) )
            && !(right = go_right ( first )) ) {
          do {
            go_up ( second );
            --level;
          } while ( go_up ( first ) && !(right = go_right ( first )) );
        }

        if ( right && second_agrees ) {
          go_up ( second );
          --level;
        }

        second_agrees = true;
        if ( right || down ) {
          if ( parent_edge_label ( first ) == 'N' ) {
            second_agrees = false;
          }
          else {
            second_agrees = go_down ( second, parent_edge_label ( first ) );
          }
          if ( second_agrees ) ++level;
        }
      } while ( !is_root( first ) );
    }  /* -----  end of template function kmer_exact_matches  ----- */

  /* END OF Top-down fine interator meta-function declarations  ------------------ */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef INDEX_ITERATOR_H__  ----- */
