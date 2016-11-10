/*
 * =====================================================================================
 *
 * Filename: vargraph.h
 *
 * Created: Fri Nov 11, 2016  01:08
 * Last modified: Mon Nov 14, 2016  00:17
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
      // Constructors
      VarGraph(std::ifstream &ifs, std::string &name_);
      VarGraph(std::ifstream &ifs);
      VarGraph(std::string &filename, std::string &name_);
      VarGraph(const char *filename, std::string &name_);
      VarGraph(std::string &filename);
      VarGraph(const char *filename);
      VarGraph(vg::Graph &vg_graph, std::string &name_);
      // TODO: Default constructor.
      // TODO: Move/copy constructors.
      // TODO: Move/copy assignment operators.
      // Public methods
      inline void                            extend(vg::Graph &vg_graph);
      inline unsigned int                    nodes_size() const;
      inline const vg::Node&                 node_at(unsigned int idx) const;
      inline vg::Node*                       mutable_node_at(unsigned int idx);
      inline const vg::Node&                 node_by(id_t node_id) const;
      inline vg::Node*                       mutable_node_by(id_t node_id);
      inline unsigned int                    edges_size() const;
      inline const vg::Edge&                 edge_at(unsigned int idx) const;
      inline vg::Edge*                       mutable_edge_at(unsigned int idx);
      inline const std::vector< vg::Edge* >& fwd_edges(id_t node_id) const;
      inline std::vector< vg::Edge* >&       mutable_fwd_edges(id_t node_id);
      inline const std::vector< vg::Edge* >& bwd_edges(id_t node_id) const;
      inline std::vector< vg::Edge* >&       mutable_bwd_edges(id_t node_id);
      // TODO: incomplete methods for accessing paths in the graph.
      inline unsigned int                    paths_size() const;
      inline const vg::Path&                 path_at(unsigned int idx) const;
      inline vg::Path*                       mutable_path_at(unsigned int idx);
      // Attributes getters and setters
      inline const std::string&              get_name() const;
      inline vg::Graph&                      get_vg_graph();
    private:
      // Attributes
      std::string                                              name;
      vg::Graph                                                vg_graph;
      std::unordered_map< id_t, vg::Node* >                    nodes_by_id;
      std::unordered_map< id_t, std::vector< vg::Edge* >>      edges_by_id;
      std::unordered_map< id_t, std::vector< vg::Edge* >>      redges_by_id;

      // internal methods
      inline void load_file(std::ifstream &ifs);
      inline void load_file(std::string &filename);
      inline void load_file(const char *filename);
      inline bool has_node(vg::Node *node);
      inline void add_node(vg::Node *node);
      inline bool has_edge(vg::Edge *edge);
      inline void add_edge(vg::Edge *edge);
  };
}

#endif  // VARGRAPH_H__
