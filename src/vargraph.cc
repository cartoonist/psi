/*
 * =====================================================================================
 *
 * Filename: vargraph.cc
 *
 * Created: Fri Nov 11, 2016  23:12
 * Last modified: Wed Nov 16, 2016  09:31
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
#include <sstream>
#include <exception>

#include <stream/src/stream.hpp>

#include "vargraph.h"
#include "release.h"

namespace grem
{
  void
    VarGraph::extend(vg::Graph &vg_graph)
  {
    for (id_t i = 0; i < this->vg_graph.node_size(); ++i)
    {
      vg::Node *node = this->vg_graph.mutable_node(i);

      try
      {
        this->add_node(node);
      }
      catch(std::runtime_error &e)
      {
        std::cerr << "[" << PACKAGE << "] Warning: "
                  << "handling a std::runtime_error." << std::endl
                  << "While adding a node: " << e.what() << std::endl;
      }
    }

    for (id_t i = 0; i < this->vg_graph.edge_size(); ++i)
    {
      vg::Edge *edge = this->vg_graph.mutable_edge(i);

      try
      {
        this->add_edge(edge);
      }
      catch(std::runtime_error &e)
      {
        std::cerr << "[" << PACKAGE << "] Warning: "
                  << "handling a std::runtime_error." << std::endl
                  << "While adding an edge: " << e.what() << std::endl;
      }
    }

    // TODO: add paths.
  }

  void
    VarGraph::load_file(std::ifstream &ifs)
  {
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
      std::cout << "Count: " << count << std::endl;
    };

    stream::for_each(ifs, extend_graph, handle_count);
  }

  void
    VarGraph::load_file(std::string &filename)
  {
    std::ifstream ifs(filename, std::ifstream::in | std::ifstream::binary);
    this->load_file(ifs);
  }

  void
    VarGraph::load_file(const char *filename)
  {
    std::string fname(filename);
    this->load_file(fname);
  }

  bool
    VarGraph::has_node(vg::Node *node)
  {
    auto got = this->nodes_by_id.find(node->id());

    if (got == this->nodes_by_id.end()) return false;

    return true;
  }

  void
    VarGraph::add_node(vg::Node *node)
  {
    if (node->id() == 0)
    {
      throw std::runtime_error("node ID 0 is not allowed in 'vg'. Skipping.");
    }
    else if (this->has_node(node))
    {
      std::stringstream ss;
      ss << "node ID " << node->id() << " appears multiple times. Skipping.";
      throw std::runtime_error(ss.str());
    }

    vg::Node *new_node = this->vg_graph.add_node();
    *new_node = *node;
    this->nodes_by_id[new_node->id()] = new_node;
  }

  bool
    VarGraph::has_edge(vg::Edge *edge)
  {
    // TODO: check for duplication. It needs map<pair<id_t,id_t>,vg::Edge*>.
    return false;
  }

  void
    VarGraph::add_edge(vg::Edge *edge)
  {
    if (this->has_edge(edge))
    {
      std::stringstream ss;
      ss << "edge " << edge->from() << (edge->from_start() ? " start" : " end")
        << " <-> " << edge->to() << (edge->to_end() ? " end" : " start")
        << " appears multiple times. Skipping.";
      throw std::runtime_error(ss.str());
    }

    vg::Edge *new_edge = this->vg_graph.add_edge();
    *new_edge = *edge;
    this->edges_by_id[new_edge->from()].push_back(new_edge);
    this->redges_by_id[new_edge->to()].push_back(new_edge);
  }
}
