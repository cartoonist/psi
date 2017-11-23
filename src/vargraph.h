/**
 *    @file  vargraph.h
 *   @brief  VarGraph class definition.
 *
 *  This header file contains VarGraph, Paths, and VarGraph iterators class definitions
 *  and interface functions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Nov 11, 2016  01:08
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef VARGRAPH_H__
#define VARGRAPH_H__

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <unordered_set>
#include <deque>
#include <utility>
#include <functional>
#include <random>

#ifdef HAVE_MEMCPY
#  define MACRO_STACK
#  pragma push_macro("HAVE_MEMCPY")
#  undef HAVE_MEMCPY
#endif

#include <xg/xg.hpp>

#ifdef MACRO_STACK
#  pragma pop_macro("HAVE_MEMCPY")
#  undef MACRO_STACK
#endif

#include "path.h"
#include "graph_iter.h"


namespace grem
{
  class VarGraph : public xg::XG
  {
    public:
      // typedefs
      typedef vg::Node node_type;                               /**< @brief Node type. */
      typedef decltype( vg::Node().id() ) nodeid_type;          /**< @brief Node ID type. */
      typedef std::size_t rank_type;                            /**< @brief Node ID type. */
      typedef decltype( vg::Position().offset() ) offset_type;  /**< @brief Node offset type. */

      // Constructors
      VarGraph( void ) : XG( ) { }

      VarGraph( std::ifstream &ifs, bool fmt_xg = true )
      {
        if ( fmt_xg ) {
          load( ifs );
        }
        else {
          from_stream( ifs );
        }
      }

      VarGraph( vg::Graph &vg_graph ) : XG( vg_graph )
      { }

      // Public methods
        inline bool
      is_branch ( nodeid_type node_id ) const
      {
        if ( this->edges_from( node_id ).size() > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_branch  ----- */

        inline bool
      is_merge ( nodeid_type node_id ) const
      {
        if ( this->edges_to( node_id ).size() > 1 ) {
          return true;
        }
        return false;
      }  /* -----  end of method is_merge  ----- */

        inline bool
      has_edges_from( nodeid_type node_id ) const
      {
        return this->edges_from( node_id ).size() != 0;
      }

        inline bool
      has_edges_to( nodeid_type node_id ) const
      {
        return this->edges_to( node_id ).size() != 0;
      }
  };

  /* Graph interface functions  ------------------------------------------------ */

  /**
   *  @brief  Get the node ID of an ajacent node randomly.
   *
   *  @param  vargraph The variation graph.
   *  @param  node_id The ID of the node whose an adjacent node should be returned.
   *  @return an adjacent node ID if available any; otherwise 0 -- an invalid node ID.
   *
   *  Picking one of the adjacent nodes by generating a pseudo-random number with
   *  uniform distribution in [0, out-degree(v)].
   */
    inline VarGraph::nodeid_type
  get_random_adjacent ( const VarGraph &vargraph, VarGraph::nodeid_type node_id )
  {
    if ( !vargraph.has_edges_from( node_id ) ) {
      return 0;
    }

    auto fwd_edges = vargraph.edges_from(node_id);

    std::random_device rd;  // Will be used to obtain a seed for the random no. engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, fwd_edges.size() - 1);

    return fwd_edges.at(dis(gen)).to();
  }  /* -----  end of function get_random_adjacent  ----- */

  /**
   *  @brief  Get the ID of an adjacent node with least coverage.
   *
   *  @param  vargraph The variation graph.
   *  @param  node_id The ID ot the node whose an adjacent node should be returned.
   *  @param  paths_set A set of paths as a container of `Path`.
   *  @return the ID of the one of adjacent nodes with least coverage. If there are
   *          multiple nodes with minimum coverage, it returns one of them. If all nodes
   *          are covered equally or no forward edge exists, it returns zero.
   *
   *  Calculate coverage for all adjacent nodes and find the smallest one.
   */
  template< typename TContainer >
      inline VarGraph::nodeid_type
    least_covered_adjacent( const VarGraph& vargraph, VarGraph::nodeid_type node_id,
        const TContainer& paths_set )
    {
      if ( vargraph.has_edges_from( node_id ) ) {
        auto fwd_edges = vargraph.edges_from( node_id );

        VarGraph::nodeid_type first_adj_id = ( *fwd_edges.begin() ).to();

        VarGraph::nodeid_type lc_node_id = first_adj_id;
        unsigned int lc = get_path_coverage( first_adj_id, paths_set );
        bool equally_covered = true;

        for ( auto e_itr = ++fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
          const VarGraph::nodeid_type& next_node = ( *e_itr ).to();
          unsigned int next_node_cov = get_path_coverage( next_node, paths_set );

          if ( equally_covered && lc != next_node_cov ) {
            equally_covered = false;
          }

          if ( next_node_cov < lc ) {
            lc = next_node_cov;
            lc_node_id = next_node;
          }
        }

        if ( !equally_covered ) {
          return lc_node_id;
        }
      }

      return 0;
    }  /* -----  end of template function least_covered_adjacent  ----- */

  /* END OF graph interface functions  ----------------------------------------- */

  /* GRAPH ITERATORS  ============================================================ */

  /* Traits template specialization  ------------------------------------------- */

  /**
   *  @brief  Breadth-first search graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BFS for VarGraph.
   */
  template< >
    struct GraphIterTraits < VarGraph, BFS >
    {
      typedef VarGraph::nodeid_type Value;
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
      typedef struct {
        std::size_t lb_visited_rank;  /**< @brief lower-bound for rank of visited nodes. */
      } TState;
    };

  /**
   *  @brief  Backtracker graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BacktrackerIter for VarGraph.
   */
  template< >
    struct GraphIterTraits< VarGraph, Backtracker > {
      typedef VarGraph::nodeid_type Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Value > > TContainer;
      typedef void* TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        Value buffer;                           /**< @brief Buffer node ID. 0=nothing */
        bool end;                               /**< @brief End flag. */
      } TState;
    };  /* ----------  end of struct Backtracker  ---------- */

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag HaplotyperIter for VarGraph.
   */
  template< >
    struct GraphIterTraits< VarGraph, Haplotyper > {
      typedef VarGraph::nodeid_type Value;
      typedef std::deque< Value > TContainer;
      /**< @brief Set of visited paths. */
      typedef std::vector< Path< VarGraph, Compact > > TSet;
      typedef TSet::size_type Level;        /**< @brief No. of selected path so far. */
      typedef struct {
        Value start;                        /**< @brief Start node ID. */
        bool end;                           /**< @brief End flag. */
        Path< VarGraph, Compact > current_path;
        unsigned int setback;
      } TState;
    };  /* ----------  end of struct HaplotyperIter  ---------- */

  /* END OF tags template specialization  -------------------------------------- */

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Iterator template class specialization for VarGraph.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template< typename TSpec >
    class Iterator< grem::VarGraph, TSpec >
    {
      public:
        typedef grem::GraphIter< grem::VarGraph, TSpec > Type;
        /* ====================  TYPEDEFS      ======================================= */
    };  /* ----------  end of template class Iterator  ---------- */

  template< typename TSpec >
    struct Value< grem::GraphIter< grem::VarGraph, TSpec > > {
      using Type = typename grem::GraphIter< grem::VarGraph, TSpec >::TTraits::Value;
    };

  template< typename T >
    struct Level;

  template< typename TSpec >
    struct Level< grem::GraphIter< grem::VarGraph, TSpec > > {
      using Type = typename grem::GraphIter< grem::VarGraph, TSpec >::TTraits::Level;
    };
}  /* -----  end of namespace seqan  ----- */

namespace grem {
  /* VarGraph iterators template specialization  --------------------------------- */

  /* BFS template specialization  ------------------------------------------------ */

  /* Internal functions specialization. */

  /**
   *  @brief  Search the graph for next unvisited node.
   *
   *  @return the node ID of the next visited node or `0` if all are visited.
   *
   *  The lower-bound for visited rank is also updated so that any node with smaller
   *  rank is visited. So in order to find the next unvisited node, we should search
   *  among nodes with higher ranks.
   */
  template< >
      inline typename seqan::Value< GraphIter< VarGraph, BFS > >::Type
    GraphIter< VarGraph, BFS >::next_unvisited( )
    {
      TTraits::Value next_node = 0;
      while ( this->state.lb_visited_rank <= this->vargraph_ptr->max_node_rank() ) {
        next_node = this->vargraph_ptr->rank_to_id( this->state.lb_visited_rank );
        if ( !(*this)[ next_node ] ) break;
        ++this->state.lb_visited_rank;
      }
      if ( this->state.lb_visited_rank <= this->vargraph_ptr->max_node_rank() ) {
        return next_node;
      }
      else {
        return 0;
      }
    }  /* -----  end of template function get_next_unvisited  ----- */

  /* Interface functions specialization. */

  template< >
      inline bool
    at_end( GraphIter< VarGraph, BFS >& it )
    {
      return it.visiting_buffer.empty();
    }

  template< >
      inline GraphIter< VarGraph, BFS >
    begin( const VarGraph& g,
        typename seqan::Value< GraphIter< VarGraph, BFS > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, BFS > >::Type TValue;

      GraphIter< VarGraph, BFS > begin_itr;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.state.lb_visited_rank = 1;
      if ( g.id_to_rank( start_node_id ) == 1 ) {
        ++begin_itr.state.lb_visited_rank;
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.visiting_buffer.push_back( std::make_pair( start_node_id, 0 ) );
      begin_itr.visited.insert( std::make_pair( start_node_id, 0 ) );
      begin_itr.itr_value = begin_itr.visiting_buffer.front().first;

      return begin_itr;
    }

  template< >
      inline void
    go_begin( GraphIter< VarGraph, BFS >& it,
        typename seqan::Value< GraphIter< VarGraph, BFS > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, BFS > >::Type TValue;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.vargraph_ptr->rank_to_id( 1 );
      }

      it.state.lb_visited_rank = 1;
      if ( it.vargraph_ptr->id_to_rank( start_node_id ) == 1 ) {
        ++it.state.lb_visited_rank;
      }

      it.visiting_buffer.clear();
      it.visiting_buffer.push_back( std::make_pair( start_node_id, 0 ) );
      it.visited.clear();
      it.visited.insert( std::make_pair( start_node_id, 0 ) );
      it.itr_value = it.visiting_buffer.front().first;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, BFS > >::Type
    level( GraphIter< VarGraph, BFS >& it )
    {
      if ( !it.visiting_buffer.empty() ) {
        return it.visiting_buffer.front().second;
      }

      throw std::runtime_error( "invalid level query on the end of iterator." );
    }

  /* Member functions specialization. */

  template< >
      inline GraphIter< VarGraph, BFS >&
    GraphIter< VarGraph, BFS >::operator++( )
    {
      typedef typename seqan::Level< GraphIter< VarGraph, BFS > >::Type TLevel;
      typedef typename seqan::Value< GraphIter< VarGraph, BFS > >::Type TValue;

      if ( at_end( *this ) && this->raise_on_end ) {
          throw std::range_error( "The iterator has reached the end." );
      }

      if ( this->visiting_buffer.empty() ) return *this;

      TLevel plevel = level( *this );
      auto edges = this->vargraph_ptr->edges_from( this->itr_value );
      for ( auto it = edges.begin(); it != edges.end(); ++it ) {
        TValue adj_node = (*it).to();
        if ( !(*this)[ adj_node ] ) {
          this->visiting_buffer.push_back( std::make_pair( adj_node, plevel + 1 ) );
          this->visited.insert( std::make_pair( adj_node, plevel + 1 ) );
        }
      }
      this->visiting_buffer.pop_front();

      if ( !this->visiting_buffer.empty() ) {
        this->itr_value = this->visiting_buffer.front().first;
      }
      else {
        this->itr_value = this->next_unvisited();
        if ( this->itr_value != 0 ) {
          this->visiting_buffer.push_back( std::make_pair( this->itr_value, 0 ) );
          this->visited.insert( std::make_pair( this->itr_value, 0 ) );
        }
      }

      if ( this->itr_value != 0 &&
          this->state.lb_visited_rank == this->vargraph_ptr->id_to_rank( this->itr_value ) ) {
        ++this->state.lb_visited_rank;
      }

      return *this;
    }

  /**
   *  @brief  Check whether a node ID is visited by the BFS iterator or not.
   *
   *  @param  id The ID of the query node.
   *  @return `true` if node is visited by BFS; otherwise `false`.
   *
   *  It queries visited set for the node ID. Level doesn't matter so it is set to zero
   *  (see method `pair_pred`).
   */
  template< >
  template< typename TId >
      inline bool
    GraphIter< VarGraph, BFS >::operator[]( const TId& id )
    {
      return this->visited.find( std::make_pair( id, 0 ) ) != this->visited.end();
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator[]  ----- */

  /* END OF BFS template specialization  ----------------------------------------- */

  /* Backtracker template specialization  ---------------------------------------- */

  /* Interface functions specialization. */
  template< >
      inline bool
    at_end( GraphIter< VarGraph, Backtracker >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline GraphIter< VarGraph, Backtracker >
    begin( const VarGraph& g,
        typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type TValue;

      GraphIter< VarGraph, Backtracker > begin_itr;
      TValue start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.itr_value = start_node_id;
      begin_itr.state.buffer = 0;
      begin_itr.state.end = false;
      begin_itr.state.start = start_node_id;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Backtracker >& it,
        typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type TValue;

      TValue start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.state.start;  // Re-use stored start node.
      }

      it.itr_value = start_node_id;
      it.state.buffer = 0;  // Re-set buffer.
      it.state.end = false;  // Re-set at-end flag.
      it.visiting_buffer.clear();
    }  /* -----  end of template function go_begin  ----- */

  /* Member functions specialization. */

  template< >
      inline GraphIter< VarGraph, Backtracker >&
    GraphIter< VarGraph, Backtracker >::operator++( )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Backtracker > >::Type TValue;

      if ( this->state.end && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( this->state.buffer != 0 ) {                             // Any node buffered?
        this->itr_value = this->state.buffer;                      // Use it.
        this->state.buffer = 0;                                    // Clear up buffer.
      }
      else {                                                  // else
        TValue cnode_id = this->itr_value;
        if ( this->vargraph_ptr->has_edges_from( cnode_id ) ) {  // Any forward edge?
          // Go forward.
          this->itr_value = this->vargraph_ptr->edges_from( cnode_id ).at( 0 ).to();
          // On each branch nodes enqueue other branches for traversing later.
          auto edges = this->vargraph_ptr->edges_from( cnode_id );
          for ( int i = edges.size() - 1; i >= 1; --i ) {
            this->visiting_buffer.push_back( std::make_pair( cnode_id, edges[i].to() ) );
          }
        }
        else {
          this->state.end = true;
        }
      }

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Backtracker >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Backtracker >&
    GraphIter< VarGraph, Backtracker >::operator--( )
    {
      if ( this->state.buffer != 0 ) {                             // Any node buffered?
        while (                // Remove all buffered branches of the current node.
            !this->visiting_buffer.empty() &&
            this->visiting_buffer.back().first == this->itr_value ) {
          this->visiting_buffer.pop_back();
        }
        this->state.buffer = 0;
      }

      if ( !this->visiting_buffer.empty() ) {                 // Go back in buffer.
        this->itr_value = this->visiting_buffer.back().first;
        this->state.buffer = this->visiting_buffer.back().second;
        this->visiting_buffer.pop_back();
        this->state.end = false;  // Reset at-end flag.
      }
      else {
        this->state.end = true;  // Set at-end flag.
      }

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Backtracker >::operator--  ----- */

  /* END OF Backtracker template specialization  --------------------------------- */

  /* Haplotyper template specialization  ----------------------------------------- */

  /* Interface functions specialization. */
  template< >
      inline bool
    at_end( GraphIter< VarGraph, Haplotyper >& it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper >
    begin( const VarGraph& g,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type TValue;

      GraphIter< VarGraph, Haplotyper > begin_itr;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.itr_value = start_node_id;
      begin_itr.state.start = start_node_id;
      begin_itr.state.end = false;
      add_node( begin_itr.state.current_path, begin_itr.itr_value );
      begin_itr.state.setback = 0;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  template< >
      inline void
    go_begin( GraphIter< VarGraph, Haplotyper >& it,
        typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type start )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type TValue;

      TValue start_node_id;
      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.state.start;  // Re-use start node.
      }

      it.itr_value = start_node_id;
      it.state.start = start_node_id;
      it.visiting_buffer.clear();
      it.state.end = false;  // Re-set at-end flag.
      it.visited.clear();
      clear( it.state.current_path );
      add_node( it.state.current_path, it.itr_value );
      it.state.setback = 0;
    }  /* -----  end of template function go_begin  ----- */

  template< >
      inline typename seqan::Level< GraphIter< VarGraph, Haplotyper > >::Type
    level( GraphIter< VarGraph, Haplotyper >& it )
    {
      return it.visited.size();
    }

  /* Member functions specialization. */

  template< >
      inline void
    GraphIter< VarGraph, Haplotyper >::set_setback( )
    {
      this->state.setback = (( this->visited.size() == 0 /* first path */||
                               this->visited.size() % 2 /* odd */) ?
                             this->visited.size() : this->visited.size() + 1 );
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::set_setback  ----- */

  /**
   *  A setback path is a sequence of last 's' nodes of currently generating haplotype.
   *  An unvisited setback path is a path that does not occur as a subset of any
   *  previously generated haplotypes. This function search in adjacent nodes set for
   *  a node that together with 's-1' previously selected nodes forms an unvisited
   *  setback in order to cover more k-mers from all paths in the graph; i.e. generating
   *  more diverse haplotypes.
   */
  template< >
      inline GraphIter< VarGraph, Haplotyper >&
    GraphIter< VarGraph, Haplotyper >::operator++( )
    {
      typedef typename seqan::Value< GraphIter< VarGraph, Haplotyper > >::Type TValue;

      if ( this->state.end && this->raise_on_end ) {
        throw std::range_error( "The iterator has reached the end." );
      }

      if ( !this->vargraph_ptr->has_edges_from( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      if ( this->state.setback != 0 &&
          this->visiting_buffer.size() >= this->state.setback ) {
        this->visiting_buffer.pop_front();
      }

      TValue next_candidate = 0;
      auto fwd_edges = this->vargraph_ptr->edges_from( this->itr_value );
      if ( this->state.setback == 0 || fwd_edges.size() == 1 ) {
        next_candidate = fwd_edges[0].to();
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        for ( auto e : fwd_edges ) {
          this->visiting_buffer.push_back( e.to() );
          if ( (*this)[ this->visiting_buffer ] ) {  // Visited?
            this->visiting_buffer.pop_back();
            continue;                             // Next edge.
          }
          this->visiting_buffer.pop_back();       // No change to the iterator state.
          next_candidate = e.to();                // Found!
        }
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( next_candidate == 0 ) {
        next_candidate = least_covered_adjacent( *this->vargraph_ptr,
            this->itr_value, this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( next_candidate == 0 ) {
        next_candidate =
          get_random_adjacent( ( *this->vargraph_ptr ),  this->itr_value );
      }

      this->itr_value = next_candidate;
      if ( this->state.setback != 0 ) {
        this->visiting_buffer.push_back( this->itr_value );
      }
      add_node( this->state.current_path, this->itr_value );

      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator++  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper >&
    GraphIter< VarGraph, Haplotyper >::operator--( int )
    {
      this->itr_value = this->state.start;    // Reset the iterator to the start node.
      this->visiting_buffer.clear();
      if ( this->state.setback != 0 ) {
        this->visiting_buffer.push_back( this->itr_value );
      }
      this->state.end = false;                // Reset at-end flag.
      clear( this->state.current_path );
      add_node( this->state.current_path, this->itr_value );
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator--  ----- */

  template< >
      inline GraphIter< VarGraph, Haplotyper >&
    GraphIter< VarGraph, Haplotyper >::operator--( )
    {
      this->visited.push_back( this->state.current_path );
      this->set_setback();
      (*this)--;
      return *this;
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator--  ----- */

  /**
   *  @brief  Check if the given path is present in paths generated so far.
   *
   *  @param  path A container of node IDs indicating nodes in a path.
   *  @return `true` if the path is present; `false` otherwise.
   *
   *  Check whether the given path is generated before or not.
   */
  template< >
  template< typename TContainer >
      inline bool
    GraphIter< VarGraph, Haplotyper >::operator[]( const TContainer& path )
    {
      return covered_by( path.begin(), path.end(), this->visited );
    }  /* -----  end of method GraphIter< VarGraph, Haplotyper >::operator[]  ----- */

  /* END OF Haplotyper template specialization  ---------------------------------- */

  /* Haplotyper iterator interface functions  ------------------------------------ */

  /**
   *  @brief  Simulate a unique haplotype.
   *
   *  @param[out]  haplotype The simulated haplotype as a Path.
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @param[in]  tries Number of tries if the generated haplotype is not unique.
   *
   *  This function gets a Haplotyper graph iterator and generate a unique haplotype
   *  if available. The input Haplotyper iterator stores required information of the
   *  previous simulated haplotypes for which the iterator is used. So, in order to
   *  simulate multiple unique haplotypes use the same iterator as the input. It tries
   *  `tries` times to generated a unique haplotype.
   */
    inline void
  get_uniq_full_haplotype( Path< VarGraph >& haplotype,
      typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
      int tries=0 )
  {
    do {
      while ( !at_end( iter ) ) {
        add_node( haplotype, *iter );
        ++iter;
      }
      if ( tries-- && iter[ haplotype.get_nodes() ] ) {
        iter--;  // discard the traversed path and reset the Haplotyper iterator.
      }
      else {
        --iter;  // save the traversed path and reset the Haplotyper iterator.
        break;
      }
      /* trying again */
    } while ( true );
  }

  /**
   *  @brief  Add a simulated haplotype to the given pathset.
   *
   *  @param  paths The pathset.
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @param[in]  tries Number of tries if the generated haplotype is not unique.
   *
   *  @overload for `PathSet`.
   */
  template< typename TPathSet >
      inline void
    get_uniq_full_haplotype( TPathSet& paths,
        typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
        int tries=0 )
    {
      Path< VarGraph > haplotype( iter.get_vargraph() );
      get_uniq_full_haplotype( haplotype, iter, tries );
      paths.add_path( std::move( haplotype ) );
    }

  template< typename TSpec >
      inline void
    extend_to_k( Path< VarGraph, TSpec >& path,
        typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
        unsigned int k )
    {
      while ( path.get_sequence_len() < k ) {
        add_node( path, *iter );
        ++iter;
      }
    }

    inline void
  move_forward( Path< VarGraph, Dynamic >& segment,
      typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
      unsigned int k )
  {
    add_node( segment, *iter );
    ++iter;
    while ( segment.get_sequence_len() -
        iter.get_vargraph()->node_length( segment.get_nodes().front() ) >= k ) {
      pop_front( segment );
    }
  }

  template< typename TPathSet >
      inline void
    get_uniq_patches( TPathSet& paths,
        typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
        unsigned int k )
    {
      iter.raise_on_end = true;
      Path< VarGraph > patch( iter.get_vargraph() );
      Path< VarGraph, Dynamic > frontier( iter.get_vargraph() );
      try {
        while ( true ) {
          // Bootstrap.
          extend_to_k( frontier, iter, k );
          // Check the next patch distance to merge with previous patch if is less than k.
          if ( length( patch ) > 0 && iter[ frontier.get_nodes() ] ) {
            paths.add_path( patch );
            clear( patch );
          }
          // Search for a patch of length k that is not covered by visited paths of `iter`.
          while ( iter[ patch.get_nodes() ] ) move_forward( frontier, iter, k );
          // Extend the patch.
          patch += frontier;
          while( !iter[ frontier.get_nodes() ] ) {
            move_forward( frontier, iter, k );
            add_node( patch, frontier.get_nodes().back() );
          }
          clear( frontier );
        }
      }
      catch( const std::range_error& ) {
        if ( length( patch ) > 0 ) {
          paths.add_path( patch );
          --iter;  // save the traversed path and reset the Haplotyper iterator.
        }
      }
      iter.raise_on_end = false;
    }

  template< typename TPathSet >
      inline bool
    get_uniq_patched_haplotype( TPathSet& paths,
        typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
        unsigned int context_len )
    {
      if ( level( iter ) == 0 ) {
        get_uniq_full_haplotype( paths, iter );
        return true;
      }
      unsigned int paths_no = length( paths.string_set );
      get_uniq_patches( paths, iter, context_len );
      if ( paths_no == length( paths.string_set ) ) return false;
      return true;
    }

  /* END OF Haplotyper iterator interface functions  ----------------------------- */
}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */
