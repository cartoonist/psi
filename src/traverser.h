/*
 * =====================================================================================
 *
 * Filename: traverser.h
 *
 * Created: Mon Nov 14, 2016  01:11
 * Last modified: Sun Jan 29, 2017  03:12
 *
 * Description: Traversers template class.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef TRAVERSER_H__
#define TRAVERSER_H__

#include <vector>
#include <functional>

#include "vg.pb.h"
#include "vargraph.h"
#include "vargraph_iter.h"
#include "types.h"

#include <easyloggingpp/src/easylogging++.h>

// TODO: refactor: types (const, * and &).

namespace grem
{
  /* Forwards */
  template< typename TIndex, typename TIterSpec >
    class PathTraverser;

  /** Traverse interface functions **
   *  NOTE: These methods should be specialized for any other PathTraverser
   *        classes.
   **/
  template< typename TIndex, typename TIterSpec >
  void
    move_forward(PathTraverser< TIndex, TIterSpec > &ptrav,
                 std::vector< PathTraverser< TIndex, TIterSpec >> &new_ptravs)
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
        new_ptravs.push_back(PathTraverser< TIndex, TIterSpec >(ptrav, std::move(new_pos)));
      }
    }
  }

  template< typename TIndex, typename TIterSpec >
  bool
    is_finished(PathTraverser< TIndex, TIterSpec > & ptrav)
  {
    return ptrav.finished;
  }

  template< typename TIndex, typename TIterSpec >
  bool
    is_valid(PathTraverser< TIndex, TIterSpec > & ptrav)
  {
    return ptrav.is_seed_hit();
  }

  template< typename TIndex, typename TIterSpec >
  void
    get_results(PathTraverser< TIndex, TIterSpec > &ptrav,
                std::vector< typename PathTraverser< TIndex, TIterSpec >::Output > &results)
  {
    if (is_valid(ptrav))
      ptrav.get_results(results);
  }

  template< typename TIndex, typename TIterSpec >
    class PathTraverser
    {
      public:
        // Member typedefs and classes
        typedef struct
        {
          vg::Position      seed_locus;
          seqan::CharString read_id;
          unsigned int      read_pos;
        } SeedHit;

        typedef SeedHit Output;
        // Traverse parameters
        class Param
        {
          friend class PathTraverser;

          public:
          // Constructors
          Param(const ReadsChunk &reads_, unsigned int seed_len_) :
            reads(reads_), reads_index(reads.seqs), seed_len(seed_len_)
          {
            TIMED_FUNC(readsIndexTimer);
          }

          // Attributes getters and setters
          inline const ReadsChunk     &get_reads()
          { return this->reads; }

          inline const DnaSeqSetIndex< TIndex > &get_reads_index()
          { return this->reads_index; }

          inline unsigned int          get_seed_len()
          { return this->seed_len; }

          private:
          ReadsChunk     reads;
          DnaSeqSetIndex< TIndex > reads_index;
          unsigned int   seed_len;
        };

        // Constructors
        PathTraverser(const VarGraph *graph,
                      PathTraverser::Param *trav_params,
                      vg::Position start) :
          vargraph(graph), parameters(trav_params), s_locus(start),
          c_locus(start), path_length(0), finished(false)
        {
          this->iters_state.push_back(
              IterState({
                DnaSSIndexIter< TIndex, TIterSpec >(this->parameters->reads_index),
                0})
              );
        }

        PathTraverser(const VarGraph &graph,
                      PathTraverser::Param &trav_params,
                      vg::Position start) :
          PathTraverser(&graph, &trav_params, start)
        {}

        PathTraverser(const PathTraverser & other)
        {
          this->vargraph = other.vargraph;
          this->parameters = other.parameters;
          this->s_locus = other.s_locus;
          this->c_locus = other.c_locus;
          this->iters_state = other.iters_state;
          this->path_length = other.path_length;
          this->finished = other.finished;
        }

        PathTraverser(PathTraverser && other) noexcept
        {
          this->vargraph = other.vargraph;
          this->parameters = other.parameters;
          this->s_locus = std::move(other.s_locus);
          this->c_locus = std::move(other.c_locus);
          this->iters_state = std::move(other.iters_state);
          this->path_length = other.path_length;
          this->finished = other.finished;
        }

        PathTraverser & operator=(const PathTraverser & other)
        {
          PathTraverser tmp(other);
          *this = std::move(tmp);
          return *this;
        }

        PathTraverser & operator=(PathTraverser && other) noexcept
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

        ~PathTraverser() noexcept {}

        PathTraverser(const PathTraverser & other, vg::Position new_locus) :
          PathTraverser(other)
        {
          this->c_locus = new_locus;
        }

        // Traverse interface functions (are friends!)
        friend void move_forward< TIndex, TIterSpec >(PathTraverser &ptrav,
                                                      std::vector< PathTraverser > &new_ptravs);
        friend bool is_finished< TIndex, TIterSpec >(PathTraverser &ptrav);
        friend bool is_valid< TIndex, TIterSpec >(PathTraverser &ptrav);
        friend void get_results< TIndex, TIterSpec >(PathTraverser &ptrav,
                                                     std::vector< PathTraverser::Output > &results);

        // Attributes getters and setters
        inline const VarGraph *              get_vargraph()
        { return this->vargraph; }

        inline const PathTraverser< TIndex, TIterSpec >::Param *  get_paramters()
        { return this->parameters; }

        inline vg::Position                  get_s_locus()
        { return this->s_locus; }

        inline vg::Position                  get_c_locus()
        { return this->c_locus; }

        inline unsigned int                  get_path_length()
        { return this->path_length; }

#ifndef NDEBUG
        static inline unsigned long long int inc_total_go_down(unsigned int by=0)
        {
          static unsigned long long int total_go_down = 0;
          if (by != 0) total_go_down += by;
          return total_go_down;
        }
#endif

      private:
        // Internal typedefs and classes
        typedef struct {
          DnaSSIndexIter< TIndex, TIterSpec > iter;
          unsigned int   boffset;
        } IterState;

        // Attributes
        const VarGraph *         vargraph;        // pointer to variation graph.
        PathTraverser< TIndex, TIterSpec >::Param *   parameters;      // pointer to params (shared between traversers).
        vg::Position             s_locus;         // starting locus
        vg::Position             c_locus;         // current locus
        std::vector< IterState > iters_state;
        unsigned int             path_length;
        bool                     finished;
#ifndef NDEBUG
#endif

        // Internal methods
        inline bool is_seed_hit()
        {
          return (this->path_length == this->parameters->seed_len);
        }

        inline bool go_down(IterState &its, seqan::Value<DnaSeq>::Type c)
        {
#ifndef NDEBUG
          PathTraverser::inc_total_go_down(1);
#endif
          // XXX: assume "N" as a mismatch.
          if (c == 'N' || c == 'n') return false;

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

        inline void go_down_all(seqan::Value<DnaSeq>::Type c)
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

        inline void one_node_forward()
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

        inline void get_results(std::vector< PathTraverser< TIndex, TIterSpec >::Output > &results)
        {
          for (auto its : this->iters_state)
          {
            using TSAValue = typename seqan::SAValue< DnaSeqSetIndex< TIndex >>::Type;
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
    };

  template <class TPathTraverser>
    class GraphTraverser
    {
      public:
        // Constructors
        GraphTraverser(const VarGraph *graph,
                       const std::vector< vg::Position > *start_loci) :
          vargraph(graph)
        {
          if (start_loci != nullptr) this->starting_points = *start_loci;
        }

        GraphTraverser(const VarGraph &graph,
                       const std::vector< vg::Position > *start_loci) :
          GraphTraverser(&graph, start_loci)
        {}

        GraphTraverser(const VarGraph *graph) : GraphTraverser(graph, nullptr)
        {}

        GraphTraverser(const VarGraph &graph) : GraphTraverser(&graph)
        {}

        GraphTraverser(const GraphTraverser & other)
        {
          this->vargraph = other.vargraph;
          this->starting_points = other.starting_points;
        }

        GraphTraverser(const GraphTraverser && other) noexcept
        {
          this->vargraph = other.vargraph;
          this->starting_points = std::move(other.starting_points);
        }

        GraphTraverser & operator=(const GraphTraverser & other)
        {
          GraphTraverser tmp(other);
          *this = std::move(tmp);
          return *this;
        }

        GraphTraverser & operator=(GraphTraverser && other) noexcept
        {
          this->vargraph = other.vargraph;
          this->starting_points = std::move(other.starting_points);
        }

        ~GraphTraverser() noexcept {}

        // Public methods
        inline void add_start(vg::Position &locus)
        {
          this->starting_points.push_back(locus);
        }

        inline void traverse(typename TPathTraverser::Param trav_params,
                      std::function< void(typename TPathTraverser::Output &) > callback)
        {
          TIMED_FUNC(traverseTimer);
          unsigned int locus_counter = 0;
          unsigned int nof_reports = 0;

          for (auto locus : this->starting_points)
          {
            this->traverse_from_locus(trav_params, callback, locus);

            ++locus_counter;
            if (locus_counter % TRAVERSE_CHECKPOINT_LOCI_NO == 0)
            {
              locus_counter = 0;
              ++nof_reports;
              CLOG(INFO, "performance") << "Traversed "
                                        << nof_reports * TRAVERSE_CHECKPOINT_LOCI_NO
                                        << " starting points:";
              PERFORMANCE_CHECKPOINT(traverseTimer);
            }
          }
        }

        inline void add_all_loci(unsigned int step=1)
        {
          // TODO: mention in the documentation that the `step` is approximately preserved in
          //       the whole graph. This means that for example add_all_loci(2) would add
          //       the first loci in each node and then add other loci with distance 2 (every
          //       other loci) within the node. So at the end, it won't necessarily preserve
          //       this distance between inter-node loci.
          // **UPDATE** This algorithm use better approximation.
          // TODO: Add documentation.
          TIMED_FUNC(addAllLociTimer);

          // TODO: Old method -- remove
          /*
          for (unsigned int i = 0; i < this->vargraph->nodes_size(); ++i)
          {
            const vg::Node &node = this->vargraph->node_at(i);
            for (unsigned int j = 0; j < node.sequence().length(); j += step)
            {
              vg::Position s_point;
              s_point.set_node_id(node.id());
              s_point.set_offset(j);

              this->add_start(s_point);
            }
          }
          */

          // New algorithm
          Iterator<VarGraph, BfsIterator<>> itr(this->vargraph);

          unsigned long int prenode_remain = 0;
          unsigned long int remain_estimate = 0;
          id_t prenode_level = 0;
          std::string seq;
          while (!at_end(itr))
          {
            if (prenode_level != level(itr))
            {
              prenode_remain = remain_estimate;
              remain_estimate = 0;
              prenode_level = level(itr);
            }

            seq = this->vargraph->node_by(*itr).sequence();

            unsigned long int cursor = (step - prenode_remain) % step;
            while (cursor < seq.length())
            {
              vg::Position s_point;
              s_point.set_node_id(*itr);
              s_point.set_offset(cursor);
              this->add_start(s_point);

              cursor += step;
            }

            unsigned long int new_remain;
            if (step - prenode_remain > seq.length())
            {
              new_remain = prenode_remain + seq.length();
            }
            else
            {
              new_remain = (seq.length() - step + prenode_remain) % step;
            }

            if (remain_estimate < new_remain)
            {
              remain_estimate = new_remain;
            }

            ++itr;
          }
        }

        inline std::vector< vg::Position > const & get_starting_points()
        {
          return this->starting_points;
        }

      private:
        // Attributes
        const VarGraph *             vargraph;
        std::vector< vg::Position >  starting_points;

        // internal methods
        inline void traverse_from_locus(typename TPathTraverser::Param & trav_params,
            std::function< void(typename TPathTraverser::Output &) > & callback,
            const vg::Position & locus)
        {
          static std::vector< TPathTraverser > path_traversers;
          static std::vector< int > deleted_paths_idx;
          static std::vector< typename TPathTraverser::Output > seeds;
          static std::vector< TPathTraverser > new_ptravs;

#ifndef NDEBUG
          static long unsigned int total_nof_ptravs = 0;
          static double      avg_path_lengths = 0;
#endif

          path_traversers.push_back(TPathTraverser(this->vargraph, &trav_params, locus));
          while (!path_traversers.empty())
          {
            for (unsigned int i = 0; i < path_traversers.size(); ++i)
            {
              TPathTraverser &ptrav = path_traversers[i];
              if (is_finished(ptrav))
              {
#ifndef NDEBUG
                // XXX: compute average path length.
                avg_path_lengths = avg_path_lengths * (total_nof_ptravs / (total_nof_ptravs+1.0))
                                   + ptrav.get_path_length() / (total_nof_ptravs+1.0);
                ++total_nof_ptravs;
                if (total_nof_ptravs % AVG_GODOWNS_SAMPLES == 0)
                {
                  static int nof_reports = 1;
                  LOG(DEBUG) << "Average traversed path length (from "
                             << nof_reports * AVG_GODOWNS_SAMPLES
                             << " samples): " << avg_path_lengths;
                  ++nof_reports;
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
    };
}

#endif
