/*
 * =====================================================================================
 *
 * Filename: vargraph.h
 *
 * Created: Fri Nov 11, 2016  01:08
 * Last modified: Tue Feb 14, 2017  01:57
 *
 * Description: VarGraph class definition.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef VARGRAPH_H__
#define VARGRAPH_H__

#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>

#include "types.h"
#include "vg.pb.h"


namespace grem
{
  class VarGraph
  {
    public:
      // typedefs
      typedef vg::Node Node;

      // Constructors
      VarGraph(std::ifstream &ifs, std::string &name_) : name(name_)
      { this->extend_from_file(ifs); }

      VarGraph(std::ifstream &ifs) : name("")
      { this->extend_from_file(ifs); }

      VarGraph(std::string &filename, std::string &name_) : name(name_)
      { this->extend_from_file(filename); }

      VarGraph(const char *filename, std::string &name_) : name(name_)
      { this->extend_from_file(filename); }

      VarGraph(std::string &filename) : name("")
      { this->extend_from_file(filename); }

      VarGraph(const char *filename) : name("")
      { this->extend_from_file(filename); }

      VarGraph(vg::Graph &vg_graph, std::string &name_) : name(name_)
      { this->extend(vg_graph); }

      VarGraph() : name("") {}

      // TODO: Move/copy constructors.
      // TODO: Move/copy assignment operators.
      // TODO: Destructor.
      // Public methods
      void                                   extend(vg::Graph &vg_graph);
      void                                   extend_from_file(std::ifstream &ifs);
      void                                   extend_from_file(std::string &filename);
      void                                   extend_from_file(const char *filename);
      inline unsigned int                    nodes_size() const
      { return this->vg_graph.node_size(); }

      bool                                   has_node(const vg::Node* node) const;
      bool                                   has_node(id_t node_id) const;

      inline const vg::Node&                 node_at(unsigned int idx) const
      { return this->vg_graph.node(idx); }

      inline vg::Node*                       mutable_node_at(unsigned int idx)
      { return this->vg_graph.mutable_node(idx); }

      inline const vg::Node&                 node_by(id_t node_id) const
      { return *(this->nodes_by_id.at(node_id)); }

      inline vg::Node*                       mutable_node_by(id_t node_id)
      { return this->nodes_by_id.at(node_id); }

      inline unsigned int                    edges_size() const
      { return this->vg_graph.edge_size(); }

      bool                                   has_edge(vg::Edge *edge) const;

      inline const vg::Edge&                 edge_at(unsigned int idx) const
      { return this->vg_graph.edge(idx); }

      inline vg::Edge*                       mutable_edge_at(unsigned int idx)
      { return this->vg_graph.mutable_edge(idx); }

      bool                                   has_fwd_edge(vg::Node *node) const;
      bool                                   has_fwd_edge(id_t node_id) const;

      inline const std::vector< vg::Edge* >& fwd_edges(id_t node_id) const
      { return this->edges_by_id.at(node_id); }

      inline std::vector< vg::Edge* >&       mutable_fwd_edges(id_t node_id)
      { return this->edges_by_id.at(node_id); }

      bool                                   has_bwd_edge(vg::Node *node) const;
      bool                                   has_bwd_edge(id_t node_id) const;

      inline const std::vector< vg::Edge* >& bwd_edges(id_t node_id) const
      { return this->redges_by_id.at(node_id); }

      inline std::vector< vg::Edge* >&       mutable_bwd_edges(id_t node_id)
      { return this->redges_by_id.at(node_id); }

      // TODO: incomplete methods for accessing paths in the graph.
      inline unsigned int                    paths_size() const
      { return this->vg_graph.path_size(); }

      inline const vg::Path&                 path_at(unsigned int idx) const
      { return this->vg_graph.path(idx); }

      inline vg::Path*                       mutable_path_at(unsigned int idx)
      { return this->vg_graph.mutable_path(idx); }

      // Attributes getters and setters
      inline const std::string&              get_name() const
      { return this->name; }

      inline vg::Graph&                      get_vg_graph()
      { return this->vg_graph; }
    private:
      // Attributes
      std::string                                              name;
      vg::Graph                                                vg_graph;
      std::unordered_map< id_t, vg::Node* >                    nodes_by_id;
      std::unordered_map< id_t, std::vector< vg::Edge* >>      edges_by_id;
      std::unordered_map< id_t, std::vector< vg::Edge* >>      redges_by_id;

      // internal methods
      void add_node(vg::Node *node);
      void add_edge(vg::Edge *edge);
  };
}

#endif  // VARGRAPH_H__
