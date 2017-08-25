/**
 *    @file  traverser.h
 *   @brief  Traverser template class.
 *
 *  Traverser template class definition and its template specialisations.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Nov 14, 2016  01:11
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2016, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef TRAVERSER_H__
#define TRAVERSER_H__

// :TODO:Fri Aug 25 13:15:\@cartoonist: remove unused headers.
// :TODO:Fri Aug 25 13:15:\@cartoonist: add required headers explicitly.
#include <vector>

#include <seqan/seeds.h>

#include "vargraph.h"
#include "logger.h"
#include "sequence.h"
#include "index.h"
#include "index_iter.h"

// TODO: refactor: types (const, * and &).

namespace grem
{
  /* Forwards */
  template< typename TIndexSpec >
    class Traverser;

  /** Traverse interface functions **
   *  @note These methods should be specialized for any other Traverser
   *        classes.
   **/
  template< typename TIndexSpec >
  void
    move_forward(Traverser< TIndexSpec > &ptrav,
                 std::vector< Traverser< TIndexSpec >> &new_ptravs)
  {
    if (ptrav.finished)
      throw std::runtime_error("cannot move forward on a finalized path.");

    ptrav.one_node_forward();

    VarGraph::nodeid_type c_node_id = ptrav.c_locus.node_id();

    if (ptrav.is_seed_hit() ||
        ptrav.iters_state.empty() ||
        !ptrav.vargraph->has_edges_from(c_node_id))
    {
      ptrav.finished = true;
    }

    if (ptrav.finished) return;
    else
    {
      auto edges = ptrav.vargraph->edges_from(c_node_id);
      auto it = edges.begin();
      ptrav.c_locus.set_node_id((*it).to());
      ptrav.c_locus.set_offset(0);
      ++it;

      for (; it != edges.end(); ++it)
      {
        vg::Position new_pos;
        new_pos.set_node_id((*it).to());
        new_pos.set_offset(0);
        new_ptravs.push_back(Traverser< TIndexSpec >(ptrav, std::move(new_pos)));
      }
    }
  }

  template< typename TIndexSpec >
  bool
    is_finished(Traverser< TIndexSpec > & ptrav)
  {
    return ptrav.finished;
  }

  template< typename TIndexSpec >
  bool
    is_valid(Traverser< TIndexSpec > & ptrav)
  {
    return ptrav.is_seed_hit();
  }

  template< typename TIndexSpec >
  void
    get_results(Traverser< TIndexSpec > &ptrav,
                std::vector< typename Traverser< TIndexSpec >::Output > &results)
  {
    if (is_valid(ptrav))
      ptrav.get_results(results);
  }

  template< typename TIndexSpec >
    class Traverser
    {
      public:
        // Member typedefs and classes
        typedef struct
        {
          vg::Position      seed_locus;
          seqan::CharString read_id;
          unsigned int      read_pos;
        } SeedHit;

        // defined types
        typedef seqan::Seed < seqan::Simple > Output;
        typedef TIndexSpec IndexType;
        typedef seqan::TopDown< seqan::ParentLinks<> > IterType;
        // Traverse parameters
        class Param
        {
          friend class Traverser;

          public:
          // Constructors
          Param(const Dna5QRecords &reads_, unsigned int seed_len_)
          {
            this->reads = reads_;
            TIMED_BLOCK(readsIndexTimer, "index-reads") {
              this->reads_index = Dna5QStringSetIndex < TIndexSpec >(this->reads.str);
            }
            this->seed_len = seed_len_;
            TIMED_BLOCK(seedingTimer, "seeding") {
              this->seeds = seeding ( this->reads.str, this->seed_len,
                                      FixedLengthNonOverlapping() );
            }
            TIMED_BLOCK(seedsIndexTimer, "index-seeds") {
              this->seeds_index = Dna5QStringSetIndex < TIndexSpec >(this->seeds);
            }
          }

          // Attributes getters and setters
          inline const Dna5QRecords     &get_reads()
          { return this->reads; }

          inline const Dna5QStringSetIndex < TIndexSpec > &get_reads_index()
          { return this->reads_index; }

          inline Dna5QStringSetIndex < TIndexSpec > &mutable_get_reads_index()
          { return this->reads_index; }

          inline const Dna5QStringSet &get_seeds()
          { return this->seeds; }

          inline Dna5QStringSet &mutable_get_seeds()
          { return this->seeds; }

          inline const Dna5QStringSetIndex < TIndexSpec > &get_seeds_index()
          { return this->seeds_index; }

          inline Dna5QStringSetIndex < TIndexSpec > &mutable_get_seeds_index()
          { return this->seeds_index; }

          inline unsigned int          get_seed_len()
          { return this->seed_len; }

          private:
          Dna5QRecords     reads;
          Dna5QStringSetIndex < TIndexSpec > reads_index;
          Dna5QStringSet seeds;
          Dna5QStringSetIndex < TIndexSpec > seeds_index;
          unsigned int   seed_len;
        };

        // Constructors
        Traverser(const VarGraph *graph,
                      Traverser::Param *trav_params,
                      vg::Position start) :
          vargraph(graph), parameters(trav_params), s_locus(start),
          c_locus(start), path_length(0), finished(false)
        {
          this->iters_state.push_back(
              IterState({
                TIndexIter< Dna5QStringSetIndex < TIndexSpec >, Traverser::IterType >(this->parameters->reads_index),
                0})
              );
        }

        Traverser(const VarGraph &graph,
                      Traverser::Param &trav_params,
                      vg::Position start) :
          Traverser(&graph, &trav_params, start)
        {}

        Traverser(const Traverser & other)
        {
          this->vargraph = other.vargraph;
          this->parameters = other.parameters;
          this->s_locus = other.s_locus;
          this->c_locus = other.c_locus;
          this->iters_state = other.iters_state;
          this->path_length = other.path_length;
          this->finished = other.finished;
        }

        Traverser(Traverser && other) noexcept
        {
          this->vargraph = other.vargraph;
          this->parameters = other.parameters;
          this->s_locus = std::move(other.s_locus);
          this->c_locus = std::move(other.c_locus);
          this->iters_state = std::move(other.iters_state);
          this->path_length = other.path_length;
          this->finished = other.finished;
        }

        Traverser & operator=(const Traverser & other)
        {
          Traverser tmp(other);
          *this = std::move(tmp);
          return *this;
        }

        Traverser & operator=(Traverser && other) noexcept
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

        ~Traverser() noexcept {}

        Traverser(const Traverser & other, vg::Position new_locus) :
          Traverser(other)
        {
          this->c_locus = new_locus;
        }

        // Traverse interface functions (are friends!)
        friend void move_forward< TIndexSpec >(Traverser &ptrav,
                                                      std::vector< Traverser > &new_ptravs);
        friend bool is_finished< TIndexSpec >(Traverser &ptrav);
        friend bool is_valid< TIndexSpec >(Traverser &ptrav);
        friend void get_results< TIndexSpec >(Traverser &ptrav,
                                                     std::vector< Traverser::Output > &results);

        // Attributes getters and setters
        inline const VarGraph *              get_vargraph()
        { return this->vargraph; }

        inline const Traverser< TIndexSpec >::Param *  get_paramters()
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
          TIndexIter< Dna5QStringSetIndex < TIndexSpec >, Traverser::IterType > iter;
          unsigned int   boffset;
        } IterState;

        // Attributes
        const VarGraph *         vargraph;        // pointer to variation graph.
        Traverser< TIndexSpec >::Param *   parameters;      // pointer to params (shared between traversers).
        vg::Position             s_locus;         // starting locus
        vg::Position             c_locus;         // current locus
        std::vector< IterState > iters_state;
        unsigned int             path_length;
        bool                     finished;

        // Internal methods
        inline bool is_seed_hit()
        {
          return (this->path_length == this->parameters->seed_len);
        }

        inline bool go_down(IterState &its, seqan::Value< seqan::Dna5QString >::Type c)
        {
#ifndef NDEBUG
          Traverser::inc_total_go_down(1);
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

        inline void go_down_all(seqan::Value< seqan::Dna5QString >::Type c)
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
          VarGraph::nodeid_type c_node_id = this->c_locus.node_id();
          const VarGraph::node_type &c_node = this->vargraph->node(c_node_id);
          seqan::Dna5QString partseq = c_node.sequence().substr(this->c_locus.offset());

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

        inline void get_results(std::vector< Traverser< TIndexSpec >::Output > &results)
        {
          for (auto its : this->iters_state)
          {
            using TSAValue = typename seqan::SAValue< Dna5QStringSetIndex < TIndexSpec >>::Type;
            seqan::String<TSAValue> saPositions = getOccurrences(its.iter);
            for (unsigned i = 0; i < length(saPositions); ++i)
            {
              Traverser::Output hit;
              seqan::setBeginPositionH ( hit, this->s_locus.node_id());
              seqan::setEndPositionH ( hit, this->s_locus.offset());
              seqan::setBeginPositionV ( hit, saPositions[i].i1);  // Read ID.
              seqan::setEndPositionV ( hit, saPositions[i].i2);    // Position in the read.

              results.push_back(std::move(hit));
            }
          }
        }
    };
}

#endif  // end of TRAVERSER_H__
