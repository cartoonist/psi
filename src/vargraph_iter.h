/**
 *    @file  vargraph_iter.h
 *   @brief  Variation graph iterators.
 *
 *  Implementation of different algorithms for traversing a variation graph.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Jan 18, 2017  14:40
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef VARGRAPH_ITER_H__
#define VARGRAPH_ITER_H__

#include <unordered_set>
#include <deque>
#include <utility>

#include "vargraph.h"

namespace grem {
  /* FORWARDS  ------------------------------------------------------------------- */

  /* Generic iterator template class. */
  template<typename TContainer, typename TSpec>
    class Iterator;

  /* Graph iterator meta-functions. */
  template < typename TGraph, typename TSpec >                      // at_end
    bool at_end ( Iterator<TGraph, TSpec> it );
  template < typename TGraph, typename TSpec >                      // begin
    Iterator< TGraph, TSpec > begin ( const TGraph &g, typename TSpec::Value start=0 );
  template < typename TGraph, typename TSpec >                      // level
    typename TSpec::Level level ( Iterator< TGraph, TSpec > &it );

  /* END OF FORWARDS  ------------------------------------------------------------ */

  /**
   *  @brief  Breadth-first search graph iterator trait.
   */
  template < typename TSpec = void >
    struct BFS
    {
      typedef VarGraph::NodeID Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Level > > TContainer;

      struct pair_hash
      {
        inline std::size_t operator()(const std::pair< Value, Level > & v) const
        {
          return std::hash< Value >()(v.first);
        }
      };

      struct pair_pred
      {
        inline bool operator()
          (const std::pair< Value, Level > & v,
           const std::pair< Value, Level > & u)
          const noexcept
        {
          if (v.first == u.first) return true;
          else return false;
        }
      };

      typedef std::unordered_set< TContainer::value_type, pair_hash, pair_pred > TSet;
    };

  /**
   *  @brief  Graph iterator template class.
   *
   *  @tparam  TGraph The graph.
   *  @tparam  TSpec  Specify the traversal algorithm.
   *
   *  This class is a template specialization of general Iterator template class of
   *  containers for variation graphs.
   */
  template < typename TGraph, typename TSpec >
    class Iterator
    {
      /* ====================  FRIENDSHIP    ======================================= */

      /* Common meta-functions. */
      friend bool
        at_end < TGraph, TSpec > ( Iterator < TGraph, TSpec > it );
      friend Iterator<TGraph, TSpec>
        begin < TGraph, TSpec > ( const TGraph &g, typename TSpec::Value start );

      /* Meta-functions specific for BFS graph iterator. */
      friend typename TSpec::Level
        level < TGraph, TSpec > ( Iterator< TGraph, TSpec > &it );

      public:
        /* ====================  LIFECYCLE     ======================================= */

        Iterator ( const TGraph *graph, typename TSpec::Value start=0 )  // constructor
        {
          // Use [implicit] move assignment operator for initializing the iterator.
          *this = begin<TGraph, TSpec>(*graph, start);
        }

        /**
         *  @overload
         */
        Iterator ( const TGraph &vargraph, typename TSpec::Value start=0 ) :
          Iterator ( &vargraph, start ) { }

        /* ====================  OPERATORS     ======================================= */
        typename TSpec::Value operator* ( ) { return this->itr_value; }
        Iterator &operator++ ( );
        Iterator &operator-- ( );

      protected:
        /* ====================  METHODS       ======================================= */

        Iterator() : vargraph_ptr(nullptr) {}

        /* ====================  DATA MEMBERS  ======================================= */

        const TGraph *vargraph_ptr;                  /**< @brief Pointer to graph. */
        typename TSpec::Value itr_value;             /**< @brief Iter. current value. */
        typename TSpec::TContainer visiting_buffer;  /**< @brief Visiting buffer. */
        typename TSpec::TSet visited;                /**< @brief Visited set. */
    };  /* ----------  end of template class Iterator  ---------- */

  /* BFS template specialization  ------------------------------------------------ */

  /* Meta-functions specialization. */

  template < >
    bool at_end ( Iterator < VarGraph, BFS<> > it )
    {
      return it.visiting_buffer.empty();
    }

  template < >
    Iterator< VarGraph, BFS<> > begin ( const VarGraph &g, BFS<>::Value start )
    {
      Iterator < VarGraph, BFS<> > begin_itr;
      BFS<>::Value start_node_id;
      if (start != 0) start_node_id = start;
      else start_node_id = g.node_at(0).id();

      begin_itr.vargraph_ptr = &g;
      begin_itr.visiting_buffer.push_back(std::make_pair(start_node_id, 0));
      begin_itr.visited.insert(std::make_pair(start_node_id, 0));
      begin_itr.itr_value = begin_itr.visiting_buffer.front().first;

      return begin_itr;
    }

  template < >
    BFS <>::Level level( Iterator < VarGraph, BFS <> > & it )
    {
      if ( !it.visiting_buffer.empty() ) {
        return it.visiting_buffer.front().second;
      }

      return -1;
    }

  /* Member functions specialization. */

  template < >
    Iterator<VarGraph, BFS<>> &
    Iterator<VarGraph, BFS<>>::operator++ ( )
    {
      BFS<>::Level plevel = level(*this);
      if (this->vargraph_ptr->has_fwd_edge(this->itr_value))
      {
        auto edges = this->vargraph_ptr->fwd_edges(this->itr_value);
        for (auto it = edges.begin(); it != edges.end(); ++it)
        {
          BFS<>::Value adj_node = (*it)->to();
          if (visited.find(std::make_pair(adj_node, 0)) == // level doesn't matter (see
              visited.end())                               //   method `pair_pred`).
          {
            this->visiting_buffer.push_back(
                std::make_pair(adj_node, plevel + 1));
            if ( this->vargraph_ptr->is_merge (adj_node) ) {  // Just add merges for efficiency.
              this->visited.insert(std::make_pair(adj_node, plevel + 1));
            }
          }
        }
      }
      if ( !this->visiting_buffer.empty() ) {
        this->visiting_buffer.pop_front();
        this->itr_value = this->visiting_buffer.front().first;
      }

      return *this;
    }

  /* END OF BFS template specialization  ----------------------------------------- */

}  /* -----  end of namespace grem  ----- */
#endif  /* ----- #ifndef VARGRAPH_ITER_H__  ----- */
