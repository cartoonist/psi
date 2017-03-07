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

namespace grem {
  /**
   *  @brief  Top-down index iterator trait with finer traversal in virtual suffix tree.
   */
  template < typename TSpec = void >
    struct FinePreorder;

  /**
   *  @brief  Generic index iterator template class.
   *
   *  @tparam  TIndex The index type.
   *  @tparam  TSpec The type for further specialization of the IndexIter type.
   */
  template < typename TIndex, typename TSpec >
    class IndexIter;

  /* Forwards  ------------------------------------------------------------------- */
  template < typename TIndex >
    bool go_down (
        grem::IndexIter < TIndex, seqan::TopDown < FinePreorder<> > > &iterator,
        typename seqan::Value < TIndex >::Type c);
  template < typename TIndex >
    bool go_down_on_edge (
        grem::IndexIter < TIndex, seqan::TopDown < FinePreorder<> > > &iterator,
        typename seqan::Value < TIndex >::Type c);
  /* END OF Forwards  ------------------------------------------------------------ */

  /**
   *  @brief  IndexIter template specialization for `TopDown<grem::FinePreorder<>>`
   *          trait.
   *
   *  This specialization is a wrapper class for `Iterator < TopDown <> >` allowing to go
   *  down the virtual suffix tree finer - i.e. one character at a time - which cannot
   *  be done by regular `TopDown` iterator.
   */
  template < typename TIndex >
    class IndexIter < TIndex, seqan::TopDown < FinePreorder <> > >
    {
      friend bool
        go_down < TIndex > (
            grem::IndexIter < TIndex, seqan::TopDown < FinePreorder <> > > &iterator,
            typename seqan::Value < TIndex >::Type c);

      friend bool
        go_down_on_edge < TIndex > (
            grem::IndexIter < TIndex, seqan::TopDown < FinePreorder <> > > &iterator,
            typename seqan::Value < TIndex >::Type c);

      public:
        /* ====================  LIFECYCLE     ======================================= */
        IndexIter ( TIndex index ) :                             /* constructor */
          iter_(index), boffset(0) { }

        /* ====================  ACCESSORS     ======================================= */

        /**
         *  @brief  getter function for iter_.
         */
        inline const typename seqan::Iterator < TIndex, seqan::TopDown<> >::Type &
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

        /** @brief Internal regular `TopDown` iterator. */
        typename seqan::Iterator < TIndex, seqan::TopDown<> >::Type iter_;
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
   *  @brief  Go down in virtual suffix tree by one character.
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
  template < typename TIndex >
    bool go_down (
        grem::IndexIter < TIndex, seqan::TopDown < FinePreorder <> > > &iterator,
        typename seqan::Value < TIndex >::Type c)
  {
    // :TODO:Thu Mar 02 05:21:\@cartoonist: Handle "N" outside of the function.
    if (iterator.boffset == 0) {          // iterator points a node.
      // Go down using `seqan::goDown` function updating internal iterator.
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
   *  @brief  Go down when iterator points to a char on an edge.
   *
   *  @param[in,out]  iterator The iterator of virtual suffix tree.
   *  @param[in]      c `iterator` goes down if next character on the edge is `c`.
   *  @return `true` if the path to go down exists, otherwise `false`.
   *
   *  Looking up the next character on the edge, go down by updating `boffset` if
   *  the next char on the edge is `c`. Internal iterator does not change in this
   *  function.
   */
  template < typename TIndex >
    bool go_down_on_edge (
        grem::IndexIter < TIndex, seqan::TopDown < FinePreorder <> > > &iterator,
        typename seqan::Value < TIndex >::Type c)
  {
    auto const &parent_edge_label = parentEdgeLabel(iterator.iter_);
    auto const &parent_edge_length = parentEdgeLength(iterator.iter_);
    auto const &next_char = parent_edge_label[ parent_edge_length - iterator.boffset ];
      if (c == next_char) {
        --iterator.boffset;
        return true;
      } else {
        return false;
      }
  }  /* -----  end of function go_down_on_edge  ----- */

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
   *  @brief  Iterator type specialization for `TopDown < grem::FinePreorder<> >` trait.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template < typename TIndex >
    class Iterator < TIndex, TopDown < grem::FinePreorder <> > >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename grem::IndexIter < TIndex, TopDown < grem::FinePreorder <> > >::Type Type;
    };  /* ----------  end of template class Iterator  ---------- */

}  /* -----  end of namespace seqan  ----- */

#endif  /* ----- #ifndef INDEX_ITERATOR_H__  ----- */
