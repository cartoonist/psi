/*
 * =====================================================================================
 *
 * Filename: vargraph.cc
 *
 * Created: Fri Nov 11, 2016  23:12
 * Last modified: Mon Mar 06, 2017  11:16
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
#include <random>
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
    VarGraph::extend_from_file(std::string &filename)
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

  /* Meta-functions specialization. */

  template < >
    bool at_end ( GraphIter < VarGraph, BFS<> > &it )
    {
      return it.visiting_buffer.empty();
    }

  template < >
    GraphIter< VarGraph, BFS<> > begin ( const VarGraph &g, BFS<>::Value start )
    {
      GraphIter < VarGraph, BFS<> > begin_itr;
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

  /* Meta-functions specialization. */
  template < >
    bool at_end ( GraphIter < VarGraph, Backtracker<> > &it )
    {
      return it.visiting_buffer.empty() &&
        !it.vargraph_ptr->has_fwd_edge(*it);
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
      begin_itr.visited = 0;  // Next node ID from current node. 0 = nothing buffered.

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  /* Member functions specialization. */

  template < >
    GraphIter < VarGraph, Backtracker <> > &
    GraphIter < VarGraph, Backtracker <> >::operator++ ( )
    {
      if ( this->visited != 0 ) {                             // Any node buffered?
        this->itr_value = this->visited;                      // Use it.
        this->visited = 0;                                    // Clear up buffer.
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
      }

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Backtracker <> >::operator++  ----- */

  template < >
    GraphIter < VarGraph, Backtracker <> > &
    GraphIter < VarGraph, Backtracker <> >::operator-- ( )
    {
      if ( this->visited != 0 ) {                             // Any node buffered?
        while (                // Remove all buffered branches of the current node.
            !this->visiting_buffer.empty() &&
            this->visiting_buffer.back().first == this->itr_value ) {
          this->visiting_buffer.pop_back();
        }
      }

      if ( !this->visiting_buffer.empty() ) {                 // Go back in buffer.
        this->itr_value = this->visiting_buffer.back().first;
        this->visited = this->visiting_buffer.back().second;
        this->visiting_buffer.pop_back();
      }

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Backtracker <> >::operator--  ----- */

  /* END OF Backtracker template specialization  --------------------------------- */

  /* Haplotyper template specialization  ----------------------------------------- */

  /* Meta-functions specialization. */
  template < >
    bool at_end ( GraphIter < VarGraph, Haplotyper<> > &it )
    {
      return !it.vargraph_ptr->has_fwd_edge(*it);
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
      begin_itr.visiting_buffer = start_node_id;

      return begin_itr;
    }  /* -----  end of template function begin  ----- */

  /* Member functions specialization. */

  template < >
    GraphIter < VarGraph, Haplotyper <> > &
    GraphIter < VarGraph, Haplotyper <> >::operator++ ( )
    {
      Haplotyper<>::Value cnode_id = this->itr_value;
      if ( !this->vargraph_ptr->has_fwd_edge ( cnode_id ) ) {    // No forward edges?
        return *this;                                            // Return.
      }

      auto fwd_edges = this->vargraph_ptr->fwd_edges(cnode_id);  // Forward edges.
      // Search for a forward node that is not in visited branches.
      for ( auto e_itr = fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
        const Haplotyper<>::Value &next_node = ( *e_itr )->to();
        if ( this->visited.find( next_node )                     // Visited?
            != this->visited.end() ) {
          continue;                                              // Next edge.
        }

        this->itr_value = next_node;                             // Not visited? Use it.
        // Only nodes whose parent is a branch node are added to the visited node set.
        if ( this->vargraph_ptr->is_branch ( cnode_id ) ) {      // Parent is branch?
          this->visited.insert ( next_node );                    // Add to visited.
        }
        return *this;                                            // Found! Return.
      }

      // If all forward edges are visited, pick one randomly with uniform distribution.
      std::random_device rd;  // Will be used to obtain a seed for the random no. engine
      std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
      std::uniform_int_distribution<> dis(0, fwd_edges.size() - 1);
      this->itr_value = dis(gen);

      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper <> >::operator++  ----- */

  template < >
    GraphIter < VarGraph, Haplotyper <> > &
    GraphIter < VarGraph, Haplotyper <> >::operator-- ( )
    {
      this->itr_value = this->visiting_buffer;  // Reset the iterator to the start node.
      return *this;
    }  /* -----  end of method GraphIter < VarGraph, Haplotyper <> >::operator--  ----- */

  /* END OF Haplotyper template specialization  ---------------------------------- */

  /* Haplotyper iterator meta-functions  ----------------------------------------- */

  /**
   *  @brief  Simulate a unique haplotype.
   *
   *  @param[in,out]  iter Haplotyper graph iterator.
   *  @return Node IDs list of the simulated haplotype.
   *
   *  This function gets a Haplotyper graph iterator and generate a unique haplotype
   *  if available. The input Haplotyper iterator stores required information of the
   *  previous simulated haplotypes for which the iterator is used. So, in order to
   *  simulate multiple unique haplotypes use the same iterator as the input.
   */
  inline std::vector < VarGraph::NodeID >
    get_uniq_haplotype ( typename seqan::Iterator < VarGraph, Haplotyper<> >::Type &iter )
    {
      std::vector < VarGraph::NodeID > haplotype;
      --iter;                                 // reset the Haplotyper iterator.
      while ( !at_end ( iter ) ) {
        haplotype.push_back ( *iter );
        ++iter;
      }
      return haplotype;
    }

  /* END OF Haplotyper iterator meta-functions  ---------------------------------- */
}
