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

#include <seqan/basic.h>

#include "graph_iter.h"


namespace grem
{
  class VarGraph : public xg::XG
  {
    public:
      // typedefs
      typedef vg::Node node_type;               /**< @brief Node type. */
      typedef size_t nodeid_type;               /**< @brief Node ID type. */
      typedef size_t rank_type;                 /**< @brief Node ID type. */
      typedef std::unordered_set< nodeid_type > NodeCoverage;

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

      // Helper functions.
      std::string get_string ( std::vector < nodeid_type > &path ) const;
  };

  /* Graph interface functions  ------------------------------------------------ */

  /**
   *  @brief  Check whether the given node is on the provided path.
   *
   *  @param  node_id The ID of the node.
   *  @param  path The path as a container of `VarGraph::nodeid_type`s.
   *  @return `true` if the node ID is in the path.
   *
   *  Check if the given node ID presents in the provided path. This function assumes
   *  that the `TContainer` type has `find` and `end` methods behaving just like
   *  ordinary containers in the C++ standard library.
   */
  template < class TContainer >
    inline bool
    on_path ( VarGraph::nodeid_type node_id, const TContainer &path )
    {
      return path.find ( node_id ) != path.end();
    }

  /**
   *  @brief  Check whether a path is covered by a set of paths.
   *
   *  @param  path_begin The begin iterator of the path as set of `VarGraph::nodeid_type`.
   *  @param  path_end The end iterator of the path as set of `VarGraph::nodeid_type`.
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
            [&coverage]( VarGraph::nodeid_type i ) { return on_path ( i, coverage ); };

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
    covered_by ( const VarGraph::nodeid_type & node_id,
        const std::vector< VarGraph::NodeCoverage > &paths_coverage )
    {
      for ( const auto & coverage : paths_coverage ) {
        if ( on_path( node_id, coverage ) ) return true;
      }

      return false;
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
  get_path_coverage ( const VarGraph::nodeid_type & node_id,
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
   *  @param  paths_coverage a set of paths as a vector of `VarGraph::NodeCoverage`
   *  @return the ID of the one of adjacent nodes with least coverage. If there are
   *          multiple nodes with minimum coverage, it returns one of them. If all nodes
   *          are covered equally or no forward edge exists, it returns zero.
   *
   *  Calculate coverage for all adjacent nodes and find the smallest one.
   */
    inline VarGraph::nodeid_type
  least_covered_adjacent ( const VarGraph &vargraph, VarGraph::nodeid_type node_id,
      const std::vector< VarGraph::NodeCoverage > &paths_coverage )
  {
    if ( vargraph.has_edges_from ( node_id ) ) {
      auto fwd_edges = vargraph.edges_from( node_id );

      VarGraph::nodeid_type first_adj_id = ( *fwd_edges.begin() ).to();

      VarGraph::nodeid_type lc_node_id = first_adj_id;
      unsigned int lc = get_path_coverage ( first_adj_id, paths_coverage );
      bool equally_covered = true;

      for ( auto e_itr = ++fwd_edges.begin(); e_itr != fwd_edges.end(); ++e_itr ) {
        const VarGraph::nodeid_type &next_node = ( *e_itr ).to();
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

  template < typename TSpec = void >
    using Backtracker = BacktrackerIter < VarGraph, TSpec >;

  /**
   *  @brief  Haplotyper graph iterator tag.
   *
   *  Specialization of generic graph iterator tag HaplotyperIter for VarGraph.
   */
  template < typename TSpec >
    struct HaplotyperIter < VarGraph, TSpec > {
      typedef VarGraph::nodeid_type Value;
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
    get_uniq_haplotype ( std::vector < VarGraph::nodeid_type > &haplotype,
        typename seqan::Iterator < VarGraph, Haplotyper<> >::Type &iter,
        int tries=0 );

  /* END OF Haplotyper iterator interface function declarations  ----------------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef VARGRAPH_H__  ----- */
