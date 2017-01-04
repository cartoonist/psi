/*
 * =====================================================================================
 *
 * Filename: traverser.cc
 *
 * Created: Mon Nov 14, 2016  01:13
 * Last modified: Mon Dec 19, 2016  16:38
 *
 * Description: Traversers class implementations.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#include "traverser.h"

namespace grem
{

  /** Traverse interface functions **
   *  NOTE: These methods should be specialized for any other PathTraverser
   *        classes.
   **/
  void
    move_forward(PathTraverser &ptrav, std::vector< PathTraverser > &new_ptravs)
  {
    if (ptrav.finished)
      throw std::runtime_error("cannot move forward on a finalized path.");

    ptrav.one_node_forward();

    id_t c_node_id = ptrav.c_locus.node_id();

    if (ptrav.is_seed_hit() ||
        ptrav.iters_state.empty() ||
        !ptrav.vargraph->has_fwd_edge(c_node_id))
    {
      ptrav.finished = true;
    }

    if (ptrav.finished) return;
    else
    {
      auto edges = ptrav.vargraph->fwd_edges(c_node_id);
      auto it = edges.begin();
      ptrav.c_locus.set_node_id((*it)->to());
      ptrav.c_locus.set_offset(0);
      ++it;

      for (; it != edges.end(); ++it)
      {
        vg::Position new_pos;
        new_pos.set_node_id((*it)->to());
        new_pos.set_offset(0);
        new_ptravs.push_back(PathTraverser(ptrav, std::move(new_pos)));
      }
    }
  }

  bool
    is_finished(PathTraverser &ptrav)
  { return ptrav.finished; }

  bool
    is_valid(PathTraverser &ptrav)
  {
    return ptrav.is_seed_hit();
  }

  void
    get_results(PathTraverser &ptrav, std::vector< PathTraverser::Output > &results)
  {
    if (is_valid(ptrav))
      ptrav.get_results(results);
  }

  /**  PathTraverser  **/

  // explicit instantiation of GraphTraverser for PathTraverser.
  template class GraphTraverser< PathTraverser >;

  PathTraverser::PathTraverser(const VarGraph *graph,
                               PathTraverser::Param *trav_params,
                               vg::Position start) :
    vargraph(graph), parameters(trav_params), s_locus(start), c_locus(start),
    path_length(0), finished(false)
  {
    this->iters_state.push_back(
        IterState({
          DnaSSIndexIter(this->parameters->reads_index),
          0})
        );
  }

  PathTraverser::PathTraverser(const VarGraph &graph, PathTraverser::Param &trav_params,
                               vg::Position start) :
    PathTraverser(&graph, &trav_params, start)
  {}

  PathTraverser::PathTraverser(const PathTraverser & other)
  {
    this->vargraph = other.vargraph;
    this->parameters = other.parameters;
    this->s_locus = other.s_locus;
    this->c_locus = other.c_locus;
    this->iters_state = other.iters_state;
    this->path_length = other.path_length;
    this->finished = other.finished;
  }

  PathTraverser::PathTraverser(PathTraverser && other) noexcept
  {
    this->vargraph = other.vargraph;
    this->parameters = other.parameters;
    this->s_locus = std::move(other.s_locus);
    this->c_locus = std::move(other.c_locus);
    this->iters_state = std::move(other.iters_state);
    this->path_length = other.path_length;
    this->finished = other.finished;
  }

  PathTraverser & PathTraverser::operator=(const PathTraverser & other)
  {
    PathTraverser tmp(other);
    *this = std::move(tmp);
    return *this;
  }

  PathTraverser & PathTraverser::operator=(PathTraverser && other) noexcept
  {
    this->vargraph = other.vargraph;
    this->parameters = other.parameters;
    this->s_locus = std::move(other.s_locus);
    this->c_locus = std::move(other.c_locus);
    this->iters_state = std::move(other.iters_state);
    this->path_length = other.path_length;
    this->finished = other.finished;

    return *this;
  }

  PathTraverser::PathTraverser(const PathTraverser &other, vg::Position new_locus) :
    PathTraverser(other)
  {
    this->c_locus = new_locus;
  }

  bool
    PathTraverser::is_seed_hit()
  {
    return (this->path_length == this->parameters->seed_len);
  }

  bool
    PathTraverser::go_down(IterState &its, seqan::Value<DnaSeq>::Type c)
  {
    if (its.boffset == 0) {
      if (!seqan::goDown(its.iter, c)) return false;

      its.boffset = parentEdgeLength(its.iter) - 1;
    } else if (c == parentEdgeLabel(its.iter)[ parentEdgeLength(its.iter) - its.boffset ]) {
      --its.boffset;
    } else {
      return false;
    }

    return true;
  }

  void
    PathTraverser::go_down_all(seqan::Value<DnaSeq>::Type c)
  {
    static std::vector<int> to_be_deleted;
    for (unsigned int i = 0; i < this->iters_state.size(); ++i)
    {
      if(!this->go_down(this->iters_state[i], c)) to_be_deleted.push_back(i);
    }

    if (to_be_deleted.size() < this->iters_state.size()) ++this->path_length;

    for (auto idx : to_be_deleted)
    {
      this->iters_state.erase(this->iters_state.begin()+idx);
    }
    to_be_deleted.clear();
  }

  void
    PathTraverser::one_node_forward()
  {
    id_t c_node_id = this->c_locus.node_id();
    const vg::Node &c_node = this->vargraph->node_by(c_node_id);
    DnaSeq partseq = c_node.sequence().substr(this->c_locus.offset());

    long unsigned int i;
    for (i = 0;
        i < seqan::length(partseq) &&
        !this->iters_state.empty() &&
        this->path_length < this->parameters->get_seed_len();
        ++i)
    {
      this->go_down_all(partseq[i]);
    }
  }

  void
    PathTraverser::get_results(std::vector< PathTraverser::Output > &results)
  {
    for (auto its : this->iters_state)
    {
      typedef seqan::SAValue<DnaSeqSetIndex>::Type TSAValue;
      seqan::String<TSAValue> saPositions = getOccurrences(its.iter);
      for (unsigned i = 0; i < length(saPositions); ++i)
      {
        PathTraverser::SeedHit seed_hit;
        seed_hit.seed_locus = this->s_locus;
        seed_hit.read_id = this->parameters->reads.ids[saPositions[i].i1];
        seed_hit.read_pos = saPositions[i].i2;

        results.push_back(std::move(seed_hit));
      }
    }
  }

  /**  GraphTraverser **/
  template <class TPathTraverser>
    GraphTraverser<TPathTraverser>::GraphTraverser(const VarGraph *graph,
        const std::vector< vg::Position > *start_loci) :
      vargraph(graph)
  {
    if (start_loci != nullptr) this->starting_points = *start_loci;
  }

  template <class TPathTraverser>
    GraphTraverser<TPathTraverser>::GraphTraverser(const VarGraph &graph,
        const std::vector< vg::Position > *start_loci) :
      GraphTraverser(&graph, start_loci)
  {}

  template <class TPathTraverser>
    GraphTraverser<TPathTraverser>::GraphTraverser(const VarGraph *graph) :
      GraphTraverser(graph, nullptr)
  {}

  template <class TPathTraverser >
    GraphTraverser<TPathTraverser>::GraphTraverser(const VarGraph &graph) :
      GraphTraverser(&graph)
  {}

  template <class TPathTraverser>
  void
    GraphTraverser<TPathTraverser>::traverse(
        typename TPathTraverser::Param trav_params,
        std::function< void(typename TPathTraverser::Output &) > callback)
  {
    TIMED_FUNC(traverseTimer);
    unsigned int locus_counter = 0;

    for (auto locus : this->starting_points)
    {
      this->traverse_from_locus(trav_params, callback, locus);

      ++locus_counter;
      if (locus_counter % TRAVERSE_CHECKPOINT_LOCI_NO == 0)
      {
        locus_counter = 0;
        PERFORMANCE_CHECKPOINT(traverseTimer);
      }
    }
  }

  template <class TPathTraverser>
  void
    GraphTraverser<TPathTraverser>::add_all_loci()
  {
    TIMED_FUNC(addAllLociTimer);

    for (unsigned int i = 0; i < this->vargraph->nodes_size(); ++i)
    {
      const vg::Node &node = this->vargraph->node_at(i);
      for (unsigned int j = 0; j < node.sequence().length(); ++j)
      {
        vg::Position s_point;
        s_point.set_node_id(node.id());
        s_point.set_offset(j);

        this->add_start(s_point);
      }
    }
  }

  template < class TPathTraverser >
  void
    GraphTraverser<TPathTraverser>::traverse_from_locus(
        typename TPathTraverser::Param & trav_params,
        std::function< void(typename TPathTraverser::Output &) > & callback,
        const vg::Position & locus)
  {
    static std::vector< TPathTraverser > path_traversers;
    static std::vector< int > deleted_paths_idx;
    static std::vector< typename TPathTraverser::Output > seeds;
    static std::vector< TPathTraverser > new_ptravs;

    static long unsigned int total_nof_ptravs = 0;
    static double      avg_path_lengths = 0;

    path_traversers.push_back(TPathTraverser(this->vargraph, &trav_params, locus));
    while (!path_traversers.empty())
    {
      for (unsigned int i = 0; i < path_traversers.size(); ++i)
      {
        TPathTraverser &ptrav = path_traversers[i];
        if (is_finished(ptrav))
        {
#ifndef NDEBUG
          // XXX: compute average go downs.
          avg_path_lengths = avg_path_lengths * total_nof_ptravs /
            (total_nof_ptravs + 1) + ptrav.get_path_length() / (total_nof_ptravs + 1);
          ++total_nof_ptravs;
          if (total_nof_ptravs % AVG_GODOWNS_SAMPLES == 0)
          {
            LOG(DEBUG) << "Average number of go downs (" << AVG_GODOWNS_SAMPLES
                       << " samples): " << avg_path_lengths;
          }
#endif

          get_results(ptrav, seeds);
          for (auto s : seeds) callback(s);
          deleted_paths_idx.push_back(i);
          seeds.clear();
        }
        else
        {
          move_forward(ptrav, new_ptravs);
          path_traversers.reserve(path_traversers.size() + new_ptravs.size());
          std::move(std::begin(new_ptravs), std::end(new_ptravs),
                    std::back_inserter(path_traversers));
          new_ptravs.clear();
        }
      }

      for (auto idx : deleted_paths_idx)
      {
        path_traversers.erase(path_traversers.begin()+idx);
      }

      deleted_paths_idx.clear();
    }
  }
}
