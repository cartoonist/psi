/*
 * =====================================================================================
 *
 * Filename: vargraph.h
 *
 * Created: Fri Nov 11, 2016  01:08
 * Last modified: Fri Mar 17, 2017  10:27
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
#include <unordered_set>
#include <deque>
#include <utility>
#include <cstdint>

#include <seqan/index.h>

#include "graph_iter.h"
#include "vg.pb.h"


namespace grem
{
  class VarGraph
  {
    public:
      // typedefs
      typedef vg::Node Node;
      typedef uint64_t NodeID;   /**< @brief Node ID type. */

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
      bool                                   has_node(NodeID node_id) const;

      inline const vg::Node&                 node_at(unsigned int idx) const
      { return this->vg_graph.node(idx); }

      inline vg::Node*                       mutable_node_at(unsigned int idx)
      { return this->vg_graph.mutable_node(idx); }

      inline const vg::Node&                 node_by(NodeID node_id) const
      { return *(this->nodes_by_id.at(node_id)); }

      inline vg::Node*                       mutable_node_by(NodeID node_id)
      { return this->nodes_by_id.at(node_id); }

      inline unsigned int                    edges_size() const
      { return this->vg_graph.edge_size(); }

      bool                                   has_edge(vg::Edge *edge) const;

      inline const vg::Edge&                 edge_at(unsigned int idx) const
      { return this->vg_graph.edge(idx); }

      inline vg::Edge*                       mutable_edge_at(unsigned int idx)
      { return this->vg_graph.mutable_edge(idx); }

      bool                                   has_fwd_edge(vg::Node *node) const;
      bool                                   has_fwd_edge(NodeID node_id) const;
      bool                                   is_branch(vg::Node *node) const;
      bool                                   is_branch(NodeID node_id) const;

      inline const std::vector< vg::Edge* >& fwd_edges(NodeID node_id) const
      { return this->edges_by_id.at(node_id); }

      inline std::vector< vg::Edge* >&       mutable_fwd_edges(NodeID node_id)
      { return this->edges_by_id.at(node_id); }

      bool                                   has_bwd_edge(vg::Node *node) const;
      bool                                   has_bwd_edge(NodeID node_id) const;
      bool                                   is_merge(vg::Node *node) const;
      bool                                   is_merge(NodeID node_id) const;

      inline const std::vector< vg::Edge* >& bwd_edges(NodeID node_id) const
      { return this->redges_by_id.at(node_id); }

      inline std::vector< vg::Edge* >&       mutable_bwd_edges(NodeID node_id)
      { return this->redges_by_id.at(node_id); }

      // TODO: incomplete methods for accessing paths in the graph.
      inline unsigned int                    paths_size() const
      { return this->vg_graph.path_size(); }

      inline const vg::Path&                 path_at(unsigned int idx) const
      { return this->vg_graph.path(idx); }

      inline vg::Path*                       mutable_path_at(unsigned int idx)
      { return this->vg_graph.mutable_path(idx); }

      // Helper functions.
      std::string get_string ( std::vector < VarGraph::NodeID > &path ) const;

      // Attributes getters and setters
      inline const std::string&              get_name() const
      { return this->name; }

      inline vg::Graph&                      get_vg_graph()
      { return this->vg_graph; }
    private:
      // Attributes
      std::string                                              name;
      vg::Graph                                                vg_graph;
      std::unordered_map< NodeID, vg::Node* >                    nodes_by_id;
      std::unordered_map< NodeID, std::vector< vg::Edge* >>      edges_by_id;
      std::unordered_map< NodeID, std::vector< vg::Edge* >>      redges_by_id;

      // internal methods
      void add_node(vg::Node *node);
      void add_edge(vg::Edge *edge);
  };

  /* GRAPH ITERATORS  ============================================================ */

  /* Tags template specialization  --------------------------------------------- */

  /**
   *  @brief  Breadth-first search graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BFSIter for VarGraph.
   */
  template < typename TSpec >
    struct BFSIter < VarGraph, TSpec >
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

  template < typename TSpec = void >
    using BFS = BFSIter < VarGraph, TSpec >;

  /**
   *  @brief  Backtracker graph iterator tag.
   *
   *  Specialization of generic graph iterator tag BacktrackerIter for VarGraph.
   */
  template < typename TSpec >
    struct BacktrackerIter < VarGraph, TSpec > {
      typedef VarGraph::NodeID Value;
      typedef Value Level;
      typedef std::deque< std::pair< Value, Value > > TContainer;
      typedef std::vector< Value > TSet;        /**< @brief Used as a buffer. */
    };  /* ----------  end of struct Backtracker  ---------- */

  template < typename TSpec = void >
    using Backtracker = BacktrackerIter < VarGraph, TSpec >;

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag HaplotyperIter for VarGraph.
   */
  template < typename TSpec >
    struct HaplotyperIter < VarGraph, TSpec > {
      typedef VarGraph::NodeID Value;
      typedef Value Level;
      typedef std::vector< Value > TContainer;  /**< @brief Used to store start node. */
      typedef std::unordered_set < Value > TSet;
    };  /* ----------  end of struct HaplotyperIter  ---------- */

  template < typename TSpec = void >
    using Haplotyper = HaplotyperIter < VarGraph, TSpec >;

  /* END OF tags template specialization  -------------------------------------- */

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Iterator template class specialization for VarGraph.
   *
   *  This class extends existing Iterator class in seqan namespace.
   */
  template < typename TSpec >
    class Iterator < grem::VarGraph, TSpec >
    {
      public:
        typedef grem::GraphIter < grem::VarGraph, TSpec > Type;
        /* ====================  TYPEDEFS      ======================================= */
    };  /* ----------  end of template class Iterator  ---------- */
}  /* -----  end of namespace seqan  ----- */

namespace grem {
  /* Haplotyper iterator meta-function declarations  ----------------------------- */

  void
    get_uniq_haplotype ( std::vector < VarGraph::NodeID > &haplotype,
        typename seqan::Iterator < VarGraph, Haplotyper<> >::Type &iter );

  /* END OF Haplotyper iterator meta-function declarations  ---------------------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */
