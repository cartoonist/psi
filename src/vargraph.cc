/*
 * =====================================================================================
 *
 * Filename: vargraph.cc
 *
 * Created: Fri Nov 11, 2016  23:12
 *
 * Description: VarGraph class implementation.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#include <stdexcept>
#include <functional>
#include <iostream>
#include <ios>
#include <exception>

#include "stream/stream.h"

#include "vargraph.h"
#include "logger.h"

namespace grem
{
  /* VarGraph iterators template specialization  --------------------------------- */

  template< typename TSpec >
    using VarGraphIterTraits = typename GraphIter< VarGraph, TSpec >::TTraits;

  /* BFS template specialization  ------------------------------------------------ */

  /* Interface functions specialization. */

  template < >
    bool at_end ( GraphIter < VarGraph, BFS > &it )
    {
      return it.visiting_buffer.empty();
    }

  template < >
    GraphIter< VarGraph, BFS >
    begin ( const VarGraph &g, VarGraphIterTraits< BFS >::Value start )
    {
      typedef VarGraphIterTraits< BFS > TTraits;

      GraphIter < VarGraph, BFS > begin_itr;
      TTraits::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.rank_to_id( 1 );
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.visiting_buffer.push_back(std::make_pair(start_node_id, 0));
      begin_itr.visited.insert(std::make_pair(start_node_id, 0));
      begin_itr.itr_value = begin_itr.visiting_buffer.front().first;

      return begin_itr;
    }

  template < >
    void go_begin ( GraphIter < VarGraph, BFS > &it,
        VarGraphIterTraits< BFS >::Value start )
    {
      typedef VarGraphIterTraits< BFS > TTraits;

      TTraits::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.vargraph_ptr->rank_to_id( 1 );
      }

      it.visiting_buffer.clear();
      it.visiting_buffer.push_back ( std::make_pair ( start_node_id, 0 ) );
      it.visited.clear();
      it.visited.insert ( std::make_pair ( start_node_id, 0 ) );
      it.itr_value = it.visiting_buffer.front().first;
    }  /* -----  end of template function go_begin  ----- */

  template < >
    VarGraphIterTraits< BFS >::Level level( GraphIter < VarGraph, BFS > & it )
    {
      if ( !it.visiting_buffer.empty() ) {
        return it.visiting_buffer.front().second;
      }

      throw std::runtime_error( "invalid level query on the end of iterator." );
    }

  /* Member functions specialization. */

  template < >
    GraphIter<VarGraph, BFS> &
    GraphIter<VarGraph, BFS>::operator++ ( )
    {
      typedef VarGraphIterTraits< BFS > TTraits;

      TTraits::Level plevel = level(*this);
      if (this->vargraph_ptr->has_edges_from(this->itr_value))
      {
        auto edges = this->vargraph_ptr->edges_from(this->itr_value);
        for (auto it = edges.begin(); it != edges.end(); ++it)
        {
          TTraits::Value adj_node = (*it).to();
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

  /* Backtracker template specialization  ---------------------------------------- */

  /* Interface functions specialization. */
  template < >
    bool at_end ( GraphIter < VarGraph, Backtracker > &it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template < >
    GraphIter < VarGraph, Backtracker >
    begin ( const VarGraph &g, VarGraphIterTraits< Backtracker >::Value start )
    {
      typedef VarGraphIterTraits< Backtracker > TTraits;

      GraphIter < VarGraph, Backtracker > begin_itr;
      TTraits::Value start_node_id;

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

  template < >
    void go_begin ( GraphIter < VarGraph, Backtracker > &it,
        VarGraphIterTraits< Backtracker >::Value start )
    {
      typedef VarGraphIterTraits< Backtracker > TTraits;

      TTraits::Value start_node_id;

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

  template < >
    GraphIter < VarGraph, Backtracker > &
    GraphIter < VarGraph, Backtracker >::operator++ ( )
    {
      typedef VarGraphIterTraits< Backtracker > TTraits;

      if ( this->state.buffer != 0 ) {                             // Any node buffered?
        this->itr_value = this->state.buffer;                      // Use it.
        this->state.buffer = 0;                                    // Clear up buffer.
      }
      else {                                                  // else
        TTraits::Value cnode_id = this->itr_value;
        if ( this->vargraph_ptr->has_edges_from(cnode_id) ) {  // Any forward edge?
          // Go forward.
          this->itr_value = this->vargraph_ptr->edges_from(cnode_id).at(0).to();
          // On each branch nodes enqueue other branches for traversing later.
          auto edges = this->vargraph_ptr->edges_from(cnode_id);
          for ( int i = edges.size() - 1; i >= 1; --i ) {
            this->visiting_buffer.push_back ( std::make_pair(cnode_id, edges[i].to()) );
          }
        }
        else {
          this->state.end = true;
        }
      }

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Backtracker <> >::operator++  ----- */

  template < >
    GraphIter < VarGraph, Backtracker > &
    GraphIter < VarGraph, Backtracker >::operator-- ( )
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
    }  /* -----  end of method GraphIter < VarGraph, Backtracker <> >::operator--  ----- */

  /* END OF Backtracker template specialization  --------------------------------- */

  /* Haplotyper template specialization  ----------------------------------------- */

  /* Interface functions specialization. */
  template < >
    bool at_end ( GraphIter < VarGraph, Haplotyper > &it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper >
    begin ( const VarGraph &g, VarGraphIterTraits< Haplotyper >::Value start )
    {
      typedef VarGraphIterTraits< Haplotyper > TTraits;

      GraphIter < VarGraph, Haplotyper > begin_itr;
      TTraits::Value start_node_id;

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

  template < >
    void go_begin ( GraphIter < VarGraph, Haplotyper > &it,
        VarGraphIterTraits< Haplotyper >::Value start )
    {
      typedef VarGraphIterTraits< Haplotyper > TTraits;

      TTraits::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.state.start;  // Re-use start node.
      }

      it.itr_value = start_node_id;
      it.visiting_buffer.clear();
      it.state.end = false;  // Re-set at-end flag.
      it.visited.clear();
      clear( it.state.current_path );
      add_node( it.state.current_path, it.itr_value );
      it.state.setback = 0;
    }  /* -----  end of template function go_begin  ----- */

  /* Member functions specialization. */

  template < >
    void
    GraphIter < VarGraph, Haplotyper >::set_setback ( )
    {
      this->state.setback = (( this->visited.size() == 0 /* first path */||
                               this->visited.size() % 2 /* odd */) ?
                             this->visited.size() : this->visited.size() + 1 );
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper<> >::set_setback  ----- */

  /**
   *  A setback path is a sequence of last 's' nodes of currently generating haplotype.
   *  An unvisited setback path is a path that does not occur as a subset of any
   *  previously generated haplotypes. This function search in adjacent nodes set for
   *  a node that together with 's-1' previously selected nodes forms an unvisited
   *  setback in order to cover more k-mers from all paths in the graph; i.e. generating
   *  more diverse haplotypes.
   */
  template < >
    GraphIter < VarGraph, Haplotyper > &
    GraphIter < VarGraph, Haplotyper >::operator++ ( )
    {
      typedef VarGraphIterTraits< Haplotyper > TTraits;

      if ( !this->vargraph_ptr->has_edges_from ( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      if ( this->state.setback != 0 &&
          this->visiting_buffer.size() >= this->state.setback ) {
        this->visiting_buffer.pop_front();
      }

      TTraits::Value next_candidate = 0;
      auto fwd_edges = this->vargraph_ptr->edges_from ( this->itr_value );
      if ( this->state.setback == 0 || fwd_edges.size() == 1 ) {
        next_candidate = fwd_edges[0].to();
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        for ( auto e : fwd_edges ) {
          this->visiting_buffer.push_back ( e.to() );
          if ( covered_by ( this->visiting_buffer, this->visited ) ) {  // Visited?
            this->visiting_buffer.pop_back();
            continue;                             // Next edge.
          }
          this->visiting_buffer.pop_back();       // No change to the iterator state.
          next_candidate = e.to();                // Found!
        }
      }
      // If no unvisited setback path found, use a node with least path coverage.
      if ( next_candidate == 0 ) {
        next_candidate = least_covered_adjacent ( *this->vargraph_ptr,
            this->itr_value, this->visited );
      }
      // If all forward edges are visited, pick one randomly with uniform distribution.
      if ( next_candidate == 0 ) {
        next_candidate =
          get_random_adjacent ( ( *this->vargraph_ptr ),  this->itr_value );
      }

      this->itr_value = next_candidate;
      if ( this->state.setback != 0 ) {
        this->visiting_buffer.push_back ( this->itr_value );
      }
      add_node( this->state.current_path, this->itr_value );

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper <> >::operator++  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper > &
    GraphIter < VarGraph, Haplotyper >::operator-- ( int )
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
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper<> >::operator--  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper > &
    GraphIter < VarGraph, Haplotyper >::operator-- ( )
    {
      this->visited.push_back( this->state.current_path );
      this->set_setback();
      (*this)--;
      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper <> >::operator--  ----- */

  /**
   *  @brief  Check if the given path is present in paths generated so far.
   *
   *  @param  path A container of node IDs indicating nodes in a path.
   *  @return `true` if the path is present; `false` otherwise.
   *
   *  Check whether the given path is generated before or not.
   */
  template < >
  template < typename TContainer >
    bool
    GraphIter < VarGraph, Haplotyper >::operator[] ( const TContainer &path )
    {
      return covered_by ( path, this->visited );
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper<> >::operator[]  ----- */

  /* END OF Haplotyper template specialization  ---------------------------------- */

  /* Haplotyper iterator interface function  ------------------------------------- */

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
    void
  get_uniq_haplotype( Path<>& haplotype,
      typename seqan::Iterator< VarGraph, Haplotyper >::Type& iter,
      int tries )
  {
    do {
      clear( haplotype );
      while ( !at_end ( iter ) ) {
        add_node( haplotype, *iter );
        ++iter;
      }
      if ( tries-- && iter [ haplotype.get_nodes() ] ) {
        iter--;  // discard the traversed path and reset the Haplotyper iterator.
      }
      else {
        --iter;  // save the traversed path and reset the Haplotyper iterator.
        break;
      }
      /* trying again */
    } while (true);
  }

  /* END OF Haplotyper iterator interface function  ------------------------------ */
}
