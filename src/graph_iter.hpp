/**
 *    @file  graph_iter.hpp
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

#ifndef  PSI_GRAPH_ITER_HPP__
#define  PSI_GRAPH_ITER_HPP__

#include <seqan/basic.h>

#include "base.hpp"


namespace grem {

  /* FORWARDS  ------------------------------------------------------------------- */

  template< typename TGraph, typename TSpec >
    class GraphIter;

  /* GraphIter strategies. */
  struct BFSStrategy;
  struct DFSStrategy;
  struct BacktrackerStrategy;

  struct Random;
  struct Global;
  struct Local;
  template< typename TSpec >
    struct HaplotyperStrategy;

  /* GraphIter specialization tags. */
  typedef seqan::Tag< BFSStrategy > BFS;
  typedef seqan::Tag< DFSStrategy > DFS;
  typedef seqan::Tag< BacktrackerStrategy > Backtracker;
  template< typename TSpec = Global >
    using Haplotyper = seqan::Tag< HaplotyperStrategy< TSpec > >;

  /* Graph iterator traits. */
  template< typename TGraph, typename TSpec >
    struct GraphIterTraits;

  /* Graph iterator interface functions. */
  template< typename TGraph, typename TSpec >                      // at_end
      bool
    at_end( GraphIter<TGraph, TSpec>& it );
  template< typename TGraph, typename TSpec >                      // begin
      void
    begin( GraphIter< TGraph, TSpec >& it, const TGraph* g,
        typename GraphIter< TGraph, TSpec >::value_type start=0,
        typename GraphIter< TGraph, TSpec >::parameter_type p=GraphIter< TGraph, TSpec >::TTraits::param_default );
  template< typename TGraph, typename TSpec >
      void
    go_begin( GraphIter<TGraph, TSpec>& it,
        typename GraphIter< TGraph, TSpec >::value_type start=0,
        typename GraphIter< TGraph, TSpec >::parameter_type p=GraphIter< TGraph, TSpec >::TTraits::param_default );
  template< typename TGraph, typename TSpec >    // level
      typename GraphIter< TGraph, TSpec >::level_type
    level( GraphIter< TGraph, TSpec >& it );

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
   *  `end` interface functions, respectively.
   */
  template< typename TGraph, typename TSpec >
    class GraphIter
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef GraphIterTraits< TGraph, TSpec > TTraits;
        typedef typename TTraits::Value value_type;
        typedef typename TTraits::Level level_type;
        typedef typename TTraits::TContainer container_type;
        typedef typename TTraits::TSet set_type;
        typedef typename TTraits::TState state_type;
        typedef typename TTraits::TParameter parameter_type;
        /* ====================  DATA MEMBERS  ======================================= */
        bool raise_on_end;          /**< @brief Throw an exception if it hits the end. */
        /* ====================  ACCESSORS     ======================================= */
          inline const TGraph*
        get_vargraph( ) const
        {
          return this->vargraph_ptr;
        }
        /* ====================  INTERFACE FUNCTIONS  ================================ */
        /* Common interface functions. */
        friend bool
          at_end< TGraph, TSpec >( GraphIter< TGraph, TSpec >& it );
        friend void
          begin< TGraph, TSpec >( GraphIter< TGraph, TSpec >& it, const TGraph* g, value_type start,
              parameter_type p );
        friend void
          go_begin< TGraph, TSpec >( GraphIter<TGraph, TSpec>& it, value_type start,
              parameter_type p );
        /* Interface functions specific for BFS graph iterator and Haplotyper. */
        friend level_type
          level< TGraph, TSpec >( GraphIter< TGraph, TSpec >& it );
        /* ====================  LIFECYCLE     ======================================= */
        GraphIter( const TGraph* graph, value_type start=0, parameter_type p=TTraits::param_default )
          : GraphIter( )
        {
          grem::begin<TGraph, TSpec>( *this, graph, start, std::move( p ) );
        }

        /**
         *  @overload
         */
        GraphIter( const TGraph& vargraph, value_type start=0, parameter_type p=TTraits::param_default ) :
          GraphIter( &vargraph, start, std::move( p ) ) { }
        /* ====================  OPERATORS     ======================================= */
        inline value_type operator*( )
        {
          if ( at_end( *this ) && this->raise_on_end ) {
            throw std::range_error( "The iterator has reached the end." );
          }
          return this->itr_value;
        }
        GraphIter& operator++( );  // prefix: go forward.
        GraphIter& operator--( );  // prefix: go backward (reset) considering last traverse as history.
        GraphIter& operator--( int );  // postfix: go backward (reset) without side effects.
        template< typename T >      // traversal history query.
          bool operator[]( const T& param );
      private:
        /* ====================  METHODS       ======================================= */
        /* Private constructor. */
        GraphIter() : raise_on_end( false ), vargraph_ptr( nullptr ) { }
        /* Internal methods. */
        value_type next_unvisited( );
        void set_setback( );
        /* ====================  DATA MEMBERS  ======================================= */
        const TGraph* vargraph_ptr;      /**< @brief Pointer to graph. */
        value_type itr_value;            /**< @brief Iter. current value. */
        container_type visiting_buffer;  /**< @brief Visiting buffer. */
        set_type visited;                /**< @brief Visited set. */
        state_type state;                /**< @brief Special-purpose vars. */
        parameter_type param;            /**< @brief Input parameter. */
    };  /* ----------  end of template class GraphIter  ---------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* --- #ifndef PSI_GRAPH_ITER_HPP__  --- */
