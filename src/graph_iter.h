/**
 *    @file  graph_iter.h
 *   @brief  Generic graph iterator template classes.
 *
 *  This is a header file for generic graph iterators. All definitions/declarations are
 *  generic and they should be specialized for a specific graph type; e.g. VarGraph.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Mar 03, 2017  11:28
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  GRAPH_ITER_H__
#define  GRAPH_ITER_H__

namespace grem {

  /* FORWARDS  ------------------------------------------------------------------- */

  template < typename TGraph, typename TSpec >
    class GraphIter;

  /* Graph iterator meta-functions. */
  template < typename TGraph, typename TSpec >                      // at_end
    bool at_end ( GraphIter<TGraph, TSpec> &it );
  template < typename TGraph, typename TSpec >                      // begin
    GraphIter < TGraph, TSpec > begin ( const TGraph &g, typename TSpec::Value start=0 );
  template < typename TGraph, typename TSpec >                      // level
    typename TSpec::Level level ( GraphIter< TGraph, TSpec > &it );

  /* Graph iterator tags. */
  template < typename TGraph, typename TSpec >
    struct BFSIter;
  template < typename TGraph, typename TSpec >
    struct BacktrackerIter;
  template < typename TGraph, typename TSpec >
    struct HaplotyperIter;

  /* END OF FORWARDS  ------------------------------------------------------------ */

  /**
   *  @brief  Graph iterator template class.
   *
   *  @tparam  TGraph The graph type.
   *  @tparam  TSpec  Specify the traversal algorithm.
   *
   *  Graph iterator template class for traversing a variation graph. The iterator will
   *  go forward by operator ++ and backward by operator -- (if supported by iterator
   *  type). They can be initialized and check if it reached to the end by `begin` and
   *  `end` meta-functions, respectively.
   */
  // :TODO:Thu Mar 16 17:04:\@cartoonist: Fix the problem of determining the end of an
  //   iterator.
  template < typename TGraph, typename TSpec >
    class GraphIter
    {
      /* ====================  FRIENDSHIP    ======================================= */

      /* Common meta-functions. */
      friend bool
        at_end < TGraph, TSpec > ( GraphIter < TGraph, TSpec > &it );
      friend GraphIter<TGraph, TSpec>
        begin < TGraph, TSpec > ( const TGraph &g, typename TSpec::Value start );

      /* Meta-functions specific for BFS graph iterator. */
      friend typename TSpec::Level
        level < TGraph, TSpec > ( GraphIter< TGraph, TSpec > &it );

      public:
        /* ====================  LIFECYCLE     ======================================= */

        GraphIter ( const TGraph *graph, typename TSpec::Value start=0 )  // constructor
        {
          // Use [implicit] move assignment operator for initializing the iterator.
          *this = begin<TGraph, TSpec>(*graph, start);
        }

        /**
         *  @overload
         */
        GraphIter ( const TGraph &vargraph, typename TSpec::Value start=0 ) :
          GraphIter ( &vargraph, start ) { }

        /* ====================  OPERATORS     ======================================= */
        typename TSpec::Value operator* ( ) { return this->itr_value; }
        GraphIter &operator++ ( );
        GraphIter &operator-- ( );

      private:
        /* ====================  METHODS       ======================================= */

        GraphIter() : vargraph_ptr(nullptr) {}

        /* ====================  DATA MEMBERS  ======================================= */

        const TGraph *vargraph_ptr;                  /**< @brief Pointer to graph. */
        typename TSpec::Value itr_value;             /**< @brief Iter. current value. */
        typename TSpec::TContainer visiting_buffer;  /**< @brief Visiting buffer. */
        typename TSpec::TSet visited;                /**< @brief Visited set. */
    };  /* ----------  end of template class GraphIter  ---------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef GRAPH_ITER_H__  ----- */
