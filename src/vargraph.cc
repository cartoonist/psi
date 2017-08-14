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

#include <functional>
#include <iostream>
#include <ios>
#include <exception>

#include <stream/src/stream.hpp>

#include "vargraph.h"
#include "logger.h"

namespace grem
{
  void
    VarGraph::extend(vg::Graph &vg_graph)
  {
    for (int i = 0; i < vg_graph.node_size(); ++i)
    {
      vg::Node *node = vg_graph.mutable_node(i);

      try
      {
        this->add_node(node);
      }
      catch(std::runtime_error &e)
      {
        LOG(WARNING) << "handling a std::runtime_error "
                     << "while adding a node: " << e.what();
      }
    }

    for (int i = 0; i < vg_graph.edge_size(); ++i)
    {
      vg::Edge *edge = vg_graph.mutable_edge(i);

      try
      {
        this->add_edge(edge);
      }
      catch(std::runtime_error &e)
      {
        LOG(WARNING) << "handling a std::runtime_error "
                     << "while adding an edge: " << e.what();
      }
    }

    // TODO: add paths.
  }

  void
    VarGraph::extend_from_file(std::ifstream &ifs)
  {
    TIMED_SCOPE(loadGraphTimer, "load-graph");

    if (!ifs.is_open())
    {
      throw std::ios_base::failure("could not open the file.");
    }

    // Extend graph callback function.
    std::function< void(vg::Graph&) > extend_graph = [this] (vg::Graph& g)
    {
      this->extend(g);
    };
    // Handle count callback function.
    std::function< void(uint64_t) > handle_count = [](uint64_t count) {
      LOG(DEBUG) << "Loading " << count << " graph(s)...";
    };

    stream::for_each(ifs, extend_graph, handle_count);
  }

  void
    VarGraph::extend_from_file(const std::string &filename)
  {
    std::ifstream ifs(filename, std::ifstream::in | std::ifstream::binary);
    this->extend_from_file(ifs);
  }

  void
    VarGraph::extend_from_file(const char *filename)
  {
    std::string fname(filename);
    this->extend_from_file(fname);
  }

  bool
    VarGraph::has_node(const vg::Node* node) const
  {
    return this->has_node(node->id());
  }

  bool
    VarGraph::has_node(NodeID node_id) const
  {
    auto got = this->nodes_by_id.find(node_id);

    if (got == this->nodes_by_id.end()) return false;

    return true;
  }

  bool
    VarGraph::has_fwd_edge(vg::Node *node) const
  {
    return this->has_fwd_edge(node->id());
  }

  bool
    VarGraph::has_fwd_edge(NodeID node_id) const
  {
    auto got = this->edges_by_id.find(node_id);

    if (got == this->edges_by_id.end()) return false;

    return true;
  }

  bool
    VarGraph::is_branch ( vg::Node *node ) const
    {
      return this->is_branch(node->id());
    }  /* -----  end of method VarGraph::is_branch  ----- */

  bool
    VarGraph::is_branch ( NodeID node_id ) const
    {
      if ( this->has_fwd_edge(node_id) ) {
        return this->fwd_edges(node_id).size() > 1;
      }
      else {
        return false;
      }
    }  /* -----  end of method VarGraph::is_branch  ----- */


  bool
    VarGraph::has_bwd_edge(vg::Node *node) const
  {
    return this->has_bwd_edge(node->id());
  }

  bool
    VarGraph::has_bwd_edge(NodeID node_id) const
  {
    auto got = this->redges_by_id.find(node_id);

    if (got == this->redges_by_id.end()) return false;

    return true;
  }

  bool
    VarGraph::is_merge ( vg::Node *node ) const
    {
      return this->is_merge(node->id());
    }  /* -----  end of method VarGraph::is_merge  ----- */

  bool
    VarGraph::is_merge ( NodeID node_id ) const
    {
      if ( this->has_bwd_edge(node_id) ) {
        return this->bwd_edges(node_id).size() > 1;
      }
      else {
        return false;
      }
    }  /* -----  end of method VarGraph::is_merge  ----- */

  /**
   *  @brief  Get string value of the given path.
   *
   *  @param  path The given path as a set of node IDs.
   *  @return string representation of the path in the variation graph.
   *
   *  Get string representation of a path in the variation graph.
   */
  std::string
    VarGraph::get_string ( std::vector < VarGraph::NodeID > &path ) const
    {
      std::string repr_str;
      for ( auto it = path.begin(); it != path.end(); ++it ) {
        repr_str += this->node_by(*it).sequence();
      }
      return repr_str;
    }  /* -----  end of method VarGraph::get_string  ----- */

  void
    VarGraph::add_node(vg::Node *node)
  {
    if (node->id() == 0)
    {
      throw std::runtime_error("node ID 0 is not allowed in 'vg'. Skipping.");
    }
    else if (this->has_node(node))
    {
      std::string msg = "node ID " + std::to_string(node->id()) +
                        " appears multiple times. Skipping.";
      throw std::runtime_error(msg);
    }

    vg::Node *new_node = this->vg_graph.add_node();
    *new_node = *node;
    this->nodes_by_id[new_node->id()] = new_node;
  }

  bool
    VarGraph::has_edge(vg::Edge *edge) const
  {
    // TODO: check for duplication. It needs map<pair<NodeID,NodeID>,vg::Edge*>.
    return false;
  }

  void
    VarGraph::add_edge(vg::Edge *edge)
  {
    if (this->has_edge(edge))
    {
      std::string msg = "edge " + std::to_string(edge->from()) +
                        (edge->from_start() ? " start" : " end") +
                        " <-> " + std::to_string(edge->to()) +
                        (edge->to_end() ? " end" : " start") +
                        " appears multiple times. Skipping.";
      throw std::runtime_error(msg);
    }

    vg::Edge *new_edge = this->vg_graph.add_edge();
    *new_edge = *edge;
    this->edges_by_id[new_edge->from()].push_back(new_edge);
    this->redges_by_id[new_edge->to()].push_back(new_edge);
  }

  /* BFS template specialization  ------------------------------------------------ */

  /* Interface functions specialization. */

  template < >
    bool at_end ( GraphIter < VarGraph, BFS<> > &it )
    {
      return it.visiting_buffer.empty();
    }

  template < >
    GraphIter< VarGraph, BFS<> >
    begin ( const VarGraph &g, BFS<>::Value start )
    {
      GraphIter < VarGraph, BFS<> > begin_itr;
      BFS<>::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.node_at(0).id();
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.visiting_buffer.push_back(std::make_pair(start_node_id, 0));
      begin_itr.visited.insert(std::make_pair(start_node_id, 0));
      begin_itr.itr_value = begin_itr.visiting_buffer.front().first;

      return begin_itr;
    }

  template < >
    void go_begin ( GraphIter < VarGraph, BFS<> > &it,
        BFS<>::Value start )
    {
      BFS<>::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = it.vargraph_ptr->node_at(0).id();
      }

      it.visiting_buffer.clear();
      it.visiting_buffer.push_back ( std::make_pair ( start_node_id, 0 ) );
      it.visited.clear();
      it.visited.insert ( std::make_pair ( start_node_id, 0 ) );
      it.itr_value = it.visiting_buffer.front().first;
    }  /* -----  end of template function go_begin  ----- */

  template < >
    BFS <>::Level level( GraphIter < VarGraph, BFS <> > & it )
    {
      if ( !it.visiting_buffer.empty() ) {
        return it.visiting_buffer.front().second;
      }

      return -1;
    }

  /* Member functions specialization. */

  template < >
    GraphIter<VarGraph, BFS<>> &
    GraphIter<VarGraph, BFS<>>::operator++ ( )
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

  /* Backtracker template specialization  ---------------------------------------- */

  /* Interface functions specialization. */
  template < >
    bool at_end ( GraphIter < VarGraph, Backtracker<> > &it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template < >
    GraphIter < VarGraph, Backtracker <> >
    begin ( const VarGraph &g, Backtracker<>::Value start )
    {
      GraphIter < VarGraph, Backtracker <> > begin_itr;
      Backtracker<>::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.node_at(0).id();
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.itr_value = start_node_id;
      begin_itr.state.buffer = 0;
      begin_itr.state.end = false;
      begin_itr.state.start = start_node_id;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  template < >
    void go_begin ( GraphIter < VarGraph, Backtracker<> > &it,
        Backtracker<>::Value start )
    {
      Backtracker<>::Value start_node_id;

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
    GraphIter < VarGraph, Backtracker <> > &
    GraphIter < VarGraph, Backtracker <> >::operator++ ( )
    {
      if ( this->state.buffer != 0 ) {                             // Any node buffered?
        this->itr_value = this->state.buffer;                      // Use it.
        this->state.buffer = 0;                                    // Clear up buffer.
      }
      else {                                                  // else
        Backtracker<>::Value cnode_id = this->itr_value;
        if ( this->vargraph_ptr->has_fwd_edge(cnode_id) ) {  // Any forward edge?
          // Go forward.
          this->itr_value = this->vargraph_ptr->fwd_edges(cnode_id).at(0)->to();
          // On each branch nodes enqueue other branches for traversing later.
          auto edges = this->vargraph_ptr->fwd_edges(cnode_id);
          for ( int i = edges.size() - 1; i >= 1; --i ) {
            this->visiting_buffer.push_back ( std::make_pair(cnode_id, edges[i]->to()) );
          }
        }
        else {
          this->state.end = true;
        }
      }

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Backtracker <> >::operator++  ----- */

  template < >
    GraphIter < VarGraph, Backtracker <> > &
    GraphIter < VarGraph, Backtracker <> >::operator-- ( )
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
    bool at_end ( GraphIter < VarGraph, Haplotyper<> > &it )
    {
      return it.state.end;
    }  /* -----  end of template function at_end  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper <> >
    begin ( const VarGraph &g, Haplotyper<>::Value start )
    {
      GraphIter < VarGraph, Haplotyper <> > begin_itr;
      Haplotyper<>::Value start_node_id;

      if ( start != 0 ) {
        start_node_id = start;
      }
      else {
        start_node_id = g.node_at(0).id();
      }

      begin_itr.vargraph_ptr = &g;
      begin_itr.itr_value = start_node_id;
      begin_itr.state.start = start_node_id;
      begin_itr.state.end = false;
      begin_itr.state.current_path.insert( begin_itr.itr_value );
      begin_itr.state.setback = 0;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  template < >
    void go_begin ( GraphIter < VarGraph, Haplotyper<> > &it,
        Haplotyper<>::Value start )
    {
      Haplotyper<>::Value start_node_id;

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
      it.state.current_path.clear();
      it.state.current_path.insert( it.itr_value );
      it.state.setback = 0;
    }  /* -----  end of template function go_begin  ----- */

  /* Member functions specialization. */

  template < >
    void
    GraphIter < VarGraph, Haplotyper<> >::set_setback ( )
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
    GraphIter < VarGraph, Haplotyper <> > &
    GraphIter < VarGraph, Haplotyper <> >::operator++ ( )
    {
      if ( !this->vargraph_ptr->has_fwd_edge ( this->itr_value ) ) {
        this->state.end = true;
        return *this;
      }

      if ( this->state.setback != 0 &&
          this->visiting_buffer.size() >= this->state.setback ) {
        this->visiting_buffer.pop_front();
      }

      Haplotyper<>::Value next_candidate = 0;
      auto fwd_edges = this->vargraph_ptr->fwd_edges ( this->itr_value );
      if ( this->state.setback == 0 || fwd_edges.size() == 1 ) {
        next_candidate = fwd_edges[0]->to();
      }
      else {
        // Search for a forward node such that the setback path is not in previous paths.
        for ( auto e : fwd_edges ) {
          this->visiting_buffer.push_back ( e->to() );
          if ( covered_by ( this->visiting_buffer, this->visited ) ) {  // Visited?
            this->visiting_buffer.pop_back();
            continue;                             // Next edge.
          }
          this->visiting_buffer.pop_back();       // No change to the iterator state.
          next_candidate = e->to();               // Found!
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
      this->state.current_path.insert( this->itr_value );

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper <> >::operator++  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper <> > &
    GraphIter < VarGraph, Haplotyper<> >::operator-- ( int )
    {
      this->itr_value = this->state.start;    // Reset the iterator to the start node.
      this->visiting_buffer.clear();
      if ( this->state.setback != 0 ) {
        this->visiting_buffer.push_back( this->itr_value );
      }
      this->state.end = false;                // Reset at-end flag.
      this->state.current_path.clear();
      this->state.current_path.insert( this->itr_value );
      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper<> >::operator--  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper <> > &
    GraphIter < VarGraph, Haplotyper <> >::operator-- ( )
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
    GraphIter < VarGraph, Haplotyper<> >::operator[] ( const TContainer &path )
    {
      return covered_by ( path, this->visited );
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper<> >::operator[]  ----- */

  /* END OF Haplotyper template specialization  ---------------------------------- */

  /* Haplotyper iterator interface functions  ------------------------------------ */

  /**
   *  @brief  Simulate a unique haplotype.
   *
   *  @param[out]  haplotype The simulated haplotype as a list of node IDs.
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
    get_uniq_haplotype ( std::vector < VarGraph::NodeID > &haplotype,
        typename seqan::Iterator < VarGraph, Haplotyper<> >::Type &iter,
        int tries )
    {
      do {
        haplotype.clear();
        while ( !at_end ( iter ) ) {
          haplotype.push_back ( *iter );
          ++iter;
        }
        if ( tries-- && iter [ haplotype ] ) {
          iter--;  // discard the traversed path and reset the Haplotyper iterator.
        }
        else {
          --iter;  // save the traversed path and reset the Haplotyper iterator.
          break;
        }
        /* trying again */
      } while (true);
    }

  /* END OF Haplotyper iterator interface functions  ----------------------------- */
}
