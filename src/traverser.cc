/*
 * =====================================================================================
 *
 * Filename: traverser.cc
 *
 * Created: Mon Nov 14, 2016  01:13
 * Last modified: Mon Nov 21, 2016  00:13
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
    auto edges = ptrav.vargraph->fwd_edges(c_node_id);

    if (ptrav.path.length() == ptrav.parameters->get_seed_len() ||
        ptrav.iters_state.empty() ||
        edges.empty())
    {
      ptrav.finished = true;
    }

    if (ptrav.finished) return;
    else
    {
      auto it = edges.begin();
      ptrav.c_locus.set_node_id((*it)->to());
      ptrav.c_locus.set_offset(0);
      ++it;

      for (; it != edges.end(); ++it)
      {
        vg::Position new_pos;
        new_pos.set_node_id((*it)->to());
        new_pos.set_offset(0);
        new_ptravs.push_back(PathTraverser(ptrav, new_pos));
      }
    }

    return;
  }

  bool
    is_finished(PathTraverser &ptrav)
  { return ptrav.finished; }

  void
    get_results(PathTraverser &ptrav, std::vector< PathTraverser::Output > &results)
  {
    ptrav.get_results(results);
  }

  /**  PathTraverser  **/

  // explicit instantiation of GraphTraverser for PathTraverser.
  template class GraphTraverser< PathTraverser >;

  PathTraverser::PathTraverser(const VarGraph *graph, PathTraverser::Param *trav_params,
                               vg::Position start) :
    vargraph(graph), parameters(trav_params), s_locus(start), c_locus(start),
    finished(false)
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

  PathTraverser::PathTraverser(const PathTraverser &other, vg::Position new_locus) :
    vargraph(other.vargraph), parameters(other.parameters), s_locus(other.s_locus),
    c_locus(new_locus), iters_state(other.iters_state), path(other.path),
    finished(false)
  {}

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

    this->path.set_length(this->path.length() + 1);

    return true;
  }

  void
    PathTraverser::go_down_all(seqan::Value<DnaSeq>::Type c)
  {
    std::vector<int> to_be_deleted;
    for (unsigned int i = 0; i < this->iters_state.size(); ++i)
    {
      if(!this->go_down(this->iters_state[i], c)) to_be_deleted.push_back(i);
    }

    for (auto it = to_be_deleted.begin(); it != to_be_deleted.end(); ++it)
    {
      this->iters_state.erase(this->iters_state.begin()+(*it));
    }
  }

  void
    PathTraverser::extend_path(unsigned int visit_len)
  {
    if (visit_len == 0)
      throw std::runtime_error("cannot extend a path by the length of zero.");

    auto new_mapping = this->path.add_mapping();
    new_mapping->set_allocated_position(new vg::Position(this->c_locus));
    auto new_edit = new_mapping->add_edit();
    new_edit->set_from_length(visit_len);
    new_mapping->set_rank(this->path.mapping_size());
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
        this->path.length() < this->parameters->get_seed_len();
        ++i)
    {
      this->go_down_all(partseq[i]);
    }

    this->extend_path(i);
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
        // Create a new alignment.
        vg::Alignment new_aln;
        // Set alignment name to the read ID.
        seqan::CharString rid = this->parameters->reads.ids[saPositions[i].i1];
        new_aln.set_name(seqan::toCString(rid));
        // Set alignment's 'sequence' field.
        seqan::String<char, seqan::CStyle> seq_cstr = seqan::representative(its.iter);
        new_aln.set_sequence(seq_cstr);
        // Set alignment's 'path' field.
        auto new_path = new vg::Path(this->path);
        new_path->set_name(std::to_string(saPositions[i].i2));
        new_aln.set_allocated_path(new_path);
        // Process the new alignment.
        results.push_back(std::move(new_aln));
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
    GraphTraverser<TPathTraverser>::traverse(typename TPathTraverser::Param trav_params,
        std::function< void(typename TPathTraverser::Output &) > callback)
  {
    for (auto locus : this->starting_points)
    {
      std::vector< TPathTraverser > path_traversers;
      path_traversers.push_back(TPathTraverser(this->vargraph, &trav_params, locus));
      while (!path_traversers.empty())
      {
        for (auto it = path_traversers.begin(); it != path_traversers.end(); ++it)
        {
          if (is_finished(*it))
          {
            std::vector< vg::Alignment > seeds;
            get_results(*it, seeds);
            for (auto s : seeds) callback(s);
            it = --(path_traversers.erase(it));
          }
          else
          {
            std::vector< PathTraverser > new_ptravs;
            move_forward(*it, new_ptravs);
            for (auto npt : new_ptravs)
            {
              path_traversers.push_back(npt);
            }
          }
        }
      }
    }
  }
}
