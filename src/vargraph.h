/*
 * =====================================================================================
 *
 * Filename: vargraph.h
 *
 * Created: Fri Nov 11, 2016  01:08
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

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <utility>
#include <cstdint>
#include <random>

#include <seqan/basic.h>

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
      typedef std::unordered_set< NodeID > NodeCoverage;

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
      void                                   extend_from_file(const std::string &filename);
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

  /* Graph interface functions  ------------------------------------------------ */

  /**
   *  @brief  Check whether the given node is on the provided path.
   *
   *  @param  node_id The ID of the node.
   *  @param  path The path as a container of `VarGraph::NodeID`s.
   *  @return `true` if the node ID is in the path.
   *
   *  Check if the given node ID presents in the provided path. This function assumes
   *  that the `TContainer` type has `find` and `end` methods behaving just like
   *  ordinary containers in the C++ standard library.
   */
  template < class TContainer >
    inline bool
    on_path ( VarGraph::NodeID node_id, const TContainer &path )
    {
      return path.find ( node_id ) != path.end();
    }

  /**
   *  @brief  Check whether the given node is on the provided path.
   *
   *  @param  node The node.
   *  @param  path The path as a container of `VarGraph::NodeID`s.
   *  @return `true` if the node is in the path.
   *
   *  Check if the given node presents in the provided path. This function assumes
   *  that the `TContainer` type has `find` and `end` methods behaving just like
   *  ordinary containers in the C++ standard library.
   */
  template < class TContainer >
    inline bool
    on_path ( VarGraph::Node node, const TContainer &path )
    {
      return path.find ( static_cast<VarGraph::NodeID>( node.id() ) ) != path.end();
    }

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  path_begin The begin iterator of the path as set of `VarGraph::NodeID`.
   *  @param  path_end The end iterator of the path as set of `VarGraph::NodeID`.
   *  @param  path_set_begin The begin iterator of the paths set.
   *  @param  path_set_end The end iterator of the paths set.
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  This function checks if all nodes of the given path is covered by at least one
   *  path in the given set of paths. The path set should be a container of a type that
   *  has `find` and `end` methods behaving just like ordinary containers in the C++
   *  standard library.
   */
  template < class Iter1, class Iter2 >
    inline bool
    covered_by ( Iter1 path_begin, Iter1 path_end,
        Iter2 paths_set_begin, Iter2 paths_set_end )
    {
      if ( path_end - path_begin > 0 ) {
        for ( auto path_itr = paths_set_begin; path_itr != paths_set_end; ++path_itr ) {
          const VarGraph::NodeCoverage &coverage = *path_itr;
          auto &&on_path_wrapper =
            [&coverage]( VarGraph::NodeID i ) { return on_path ( i, coverage ); };

          if ( std::all_of ( path_begin, path_end, on_path_wrapper ) ) return true;
        }
      }

      return false;
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  path the given path to check as a container of node IDs
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return true if the given path is a subset of a path in the paths set; otherwise
   *          false -- including the case that the path is empty.
   *
   *  This function checks if all nodes of the given path is covered by at least one
   *  path in the given set of paths.
   */
  template < class TContainer1, class TContainer2 >
    inline bool
    covered_by ( const TContainer1 &path, const TContainer2 &paths_coverage )
    {
      return covered_by ( path.begin(), path.end(),
          paths_coverage.begin(), paths_coverage.end() );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a node is covered by a set of paths.
   *
   *  @param  node_id the ID of the given node to check
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return true if the node is on one of the path in the paths set, otherwise false.
   *
   *  This function simply checks if a node is on any of the path in the given set of
   *  paths.
   */
  template < >
    inline bool
    covered_by ( const VarGraph::NodeID & node_id,
        const std::vector< VarGraph::NodeCoverage > &paths_coverage )
    {
      for ( const auto & coverage : paths_coverage ) {
        if ( on_path( node_id, coverage ) ) return true;
      }

      return false;
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Check whether a node is covered by a set of paths.
   *
   *  @param  node the given node to check
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return true if the node is on one of the path in the paths set, otherwise false.
   *
   *  This function simply checks if a node is on any of the path in the given set of
   *  paths.
   */
  template < >
    inline bool
    covered_by ( const VarGraph::Node & node,
        const std::vector< VarGraph::NodeCoverage > &paths_coverage )
    {
      return covered_by ( static_cast<VarGraph::NodeID>( node.id() ), paths_coverage );
    }  /* -----  end of template function covered_by  ----- */

  /**
   *  @brief  Return the path coverage of a given node by its ID.
   *
   *  @param  node_id the node ID of the node to check
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return the number of paths that cover the given node.
   *
   *  This function get a node ID and return the number of paths that cover the node.
   */
    inline unsigned int
  get_path_coverage ( const VarGraph::NodeID & node_id,
      const std::vector< VarGraph::NodeCoverage > &paths_coverage )
  {
    unsigned int coverage = 0;
    for ( const auto & path : paths_coverage ) {
      if ( path.find( node_id ) != path.end() ) {
        ++coverage;
      }
    }
    return coverage;
  }  /* -----  end of function get_path_coverage  ----- */

  /**
   *  @brief  Return the path coverage of a given node.
   *
   *  @param  node the given node to check
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return the number of paths that cover the given node.
   *
   *  This function get a node and return the number of paths that cover the node.
   */
    inline unsigned int
  get_path_coverage ( const VarGraph::Node & node,
      const std::vector< VarGraph::NodeCoverage > &paths_coverage )
  {
    return
      get_path_coverage ( static_cast<VarGraph::NodeID>( node.id() ), paths_coverage );
  }  /* -----  end of function get_path_coverage  ----- */

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
    inline VarGraph::NodeID
  get_random_adjacent ( const VarGraph &vargraph, VarGraph::NodeID node_id )
  {
    if ( !vargraph.has_fwd_edge ( node_id ) ) {
      return 0;
    }

    auto fwd_edges = vargraph.fwd_edges(node_id);

    std::random_device rd;  // Will be used to obtain a seed for the random no. engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, fwd_edges.size() - 1);

    return fwd_edges[dis(gen)]->to();
  }  /* -----  end of function get_random_adjacent  ----- */

  /**
   *  @brief  Get the ID of an adjacent node with least coverage.
   *
   *  @param  vargraph The variation graph.
   *  @param  node_id The ID ot the node whose an adjacent node should be returned.
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return the ID of the one of adjacent nodes with least coverage. If there are
   *          multiple nodes with minimum coverage, it returns one of them. If all nodes
   *          are covered equally or no forward edge exists, it returns zero.
   *
   *  Calculate coverage for all adjacent nodes and find the smallest one.
   */
    inline VarGraph::NodeID
  least_covered_adjacent ( const VarGraph &vargraph, VarGraph::NodeID node_id,
      const std::vector< VarGraph::NodeCoverage > &paths_coverage )
  {
    if ( vargraph.has_fwd_edge ( node_id ) ) {
      auto fwd_edges = vargraph.fwd_edges( node_id );

      VarGraph::NodeID first_adj_id = ( *fwd_edges.begin() )->to();

      VarGraph::NodeID lc_node_id = first_adj_id;
      unsigned int lc = get_path_coverage ( first_adj_id, paths_coverage );
      bool equally_covered = true;

      for ( auto e_itr = ++fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
        const VarGraph::NodeID &next_node = ( *e_itr )->to();
        unsigned int next_node_cov = get_path_coverage ( next_node, paths_coverage );

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
  }  /* -----  end of function least_covered_adjacent  ----- */

  /* END OF graph interface functions  ----------------------------------------- */

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
      typedef void* TState;
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
      typedef void* TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        Value buffer;                           /**< @brief Buffer node ID. 0=nothing */
        bool end;                               /**< @brief End flag. */
      } TState;
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
      typedef std::deque< Value > TContainer;
      /**< @brief Set of visited paths. */
      typedef std::vector< VarGraph::NodeCoverage > TSet;
      typedef struct {
        Value start;                            /**< @brief Start node ID. */
        bool end;                               /**< @brief End flag. */
        VarGraph::NodeCoverage current_path;
        unsigned int setback;
      } TState;
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
  /* Haplotyper iterator interface function declarations  ------------------------ */

  void
    get_uniq_haplotype ( std::vector < VarGraph::NodeID > &haplotype,
        typename seqan::Iterator < VarGraph, Haplotyper<> >::Type &iter,
        int tries=0 );

  /* END OF Haplotyper iterator interface function declarations  ----------------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */
