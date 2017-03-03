/*
 * =====================================================================================
 *
 * Filename: vargraph.cc
 *
 * Created: Fri Nov 11, 2016  23:12
 * Last modified: Thu Mar 02, 2017  23:01
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
}
