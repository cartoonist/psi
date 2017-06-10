/*
 * =====================================================================================
 *
 * Filename: traverser.h
 *
 * Created: Mon Nov 14, 2016  01:11
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

#include <algorithm>
#include <iterator>
#include <vector>
#include <functional>
#include <unordered_set>

#include <seqan/seeds.h>

#include "vg.pb.h"
#include "vargraph.h"
#include "logger.h"
#include "sequence.h"
#include "index.h"
#include "index_iter.h"

// TODO: refactor: types (const, * and &).

namespace grem
{
  /* Forwards */
  template< typename TIndexSpec, typename TIterSpec >
    class PathTraverser;

  /** Traverse interface functions **
   *  @note These methods should be specialized for any other PathTraverser
   *        classes.
   **/
  template< typename TIndexSpec, typename TIterSpec >
  void
    move_forward(PathTraverser< TIndexSpec, TIterSpec > &ptrav,
                 std::vector< PathTraverser< TIndexSpec, TIterSpec >> &new_ptravs)
  {
    if (ptrav.finished)
      throw std::runtime_error("cannot move forward on a finalized path.");

    ptrav.one_node_forward();

    VarGraph::NodeID c_node_id = ptrav.c_locus.node_id();

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
        new_ptravs.push_back(PathTraverser< TIndexSpec, TIterSpec >(ptrav, std::move(new_pos)));
      }
    }
  }

  template< typename TIndexSpec, typename TIterSpec >
  bool
    is_finished(PathTraverser< TIndexSpec, TIterSpec > & ptrav)
  {
    return ptrav.finished;
  }

  template< typename TIndexSpec, typename TIterSpec >
  bool
    is_valid(PathTraverser< TIndexSpec, TIterSpec > & ptrav)
  {
    return ptrav.is_seed_hit();
  }

  template< typename TIndexSpec, typename TIterSpec >
  void
    get_results(PathTraverser< TIndexSpec, TIterSpec > &ptrav,
                std::vector< typename PathTraverser< TIndexSpec, TIterSpec >::Output > &results)
  {
    if (is_valid(ptrav))
      ptrav.get_results(results);
  }

  template< typename TIndexSpec, typename TIterSpec >
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

        // defined types
        typedef seqan::Seed < seqan::Simple > Output;
        typedef TIndexSpec IndexType;
        typedef TIterSpec IterType;
        // Traverse parameters
        class Param
        {
          friend class PathTraverser;

          public:
          // Constructors
          Param(const Dna5QRecords &reads_, unsigned int seed_len_)
          {
            TIMED_SCOPE(readsIndexTimer, "index-read");
            this->reads = reads_;
            this->reads_index = Dna5QStringSetIndex < TIndexSpec >(reads.str);
            this->seed_len = seed_len_;
          }

          // Attributes getters and setters
          inline const Dna5QRecords     &get_reads()
          { return this->reads; }

          inline const Dna5QStringSetIndex < TIndexSpec > &get_reads_index()
          { return this->reads_index; }

          inline Dna5QStringSetIndex < TIndexSpec > &mutable_get_reads_index()
          { return this->reads_index; }

          inline unsigned int          get_seed_len()
          { return this->seed_len; }

          private:
          Dna5QRecords     reads;
          Dna5QStringSetIndex < TIndexSpec > reads_index;
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
                TIndexIterator< Dna5QStringSetIndex < TIndexSpec >, TIterSpec >(this->parameters->reads_index),
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
        friend void move_forward< TIndexSpec, TIterSpec >(PathTraverser &ptrav,
                                                      std::vector< PathTraverser > &new_ptravs);
        friend bool is_finished< TIndexSpec, TIterSpec >(PathTraverser &ptrav);
        friend bool is_valid< TIndexSpec, TIterSpec >(PathTraverser &ptrav);
        friend void get_results< TIndexSpec, TIterSpec >(PathTraverser &ptrav,
                                                     std::vector< PathTraverser::Output > &results);

        // Attributes getters and setters
        inline const VarGraph *              get_vargraph()
        { return this->vargraph; }

        inline const PathTraverser< TIndexSpec, TIterSpec >::Param *  get_paramters()
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
          TIndexIterator< Dna5QStringSetIndex < TIndexSpec >, TIterSpec > iter;
          unsigned int   boffset;
        } IterState;

        // Attributes
        const VarGraph *         vargraph;        // pointer to variation graph.
        PathTraverser< TIndexSpec, TIterSpec >::Param *   parameters;      // pointer to params (shared between traversers).
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
          VarGraph::NodeID c_node_id = this->c_locus.node_id();
          const vg::Node &c_node = this->vargraph->node_by(c_node_id);
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

        inline void get_results(std::vector< PathTraverser< TIndexSpec, TIterSpec >::Output > &results)
        {
          for (auto its : this->iters_state)
          {
            using TSAValue = typename seqan::SAValue< Dna5QStringSetIndex < TIndexSpec >>::Type;
            seqan::String<TSAValue> saPositions = getOccurrences(its.iter);
            for (unsigned i = 0; i < length(saPositions); ++i)
            {
              PathTraverser::Output hit;
              seqan::setBeginPositionH ( hit, this->s_locus.node_id());
              seqan::setEndPositionH ( hit, this->s_locus.offset());
              seqan::setBeginPositionV ( hit, saPositions[i].i1);  // Read ID.
              seqan::setEndPositionV ( hit, saPositions[i].i2);    // Position in the read.

              results.push_back(std::move(hit));
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

        /**
         *  @brief  Pick n paths from the variation graph.
         *
         *  @param[out]  paths The set of generated paths are added to this string set.
         *  @param[out]  covered_nodes The set of covered nodes by the generated paths.
         *  @param[in]  n Number of paths.
         *
         *  This method generates a set of (probably) unique whole-genome path from the
         *  variation graph.
         */
        void
          pick_paths ( Dna5QStringSet &paths,
              std::unordered_set < VarGraph::NodeID > &covered_nodes, int n )
          {
            if ( n == 0 ) return;

            seqan::Iterator < VarGraph, Haplotyper<> >::Type hap_itr ( this->vargraph );
            seqan::Dna5QString new_path;
            std::vector < VarGraph::NodeID > new_hap;

            new_hap.reserve ( this->vargraph->nodes_size() );
            covered_nodes.reserve ( covered_nodes.size() + this->vargraph->nodes_size() );

            for ( int i = 0; i < n; ++i ) {
              get_uniq_haplotype ( new_hap, hap_itr );

              std::copy ( new_hap.begin(), new_hap.end(),
                  std::inserter ( covered_nodes, covered_nodes.end() ) );

              new_path = this->vargraph->get_string ( new_hap );

              // :TODO:Mon Mar 06 13:00:\@cartoonist: faked quality score.
              char fake_qual = 'I';
              assignQualities ( new_path, std::string ( length(new_path), fake_qual ) );

              appendValue ( paths, new_path );

              new_hap.clear();
            }

            LOG(INFO) << "Picked " << n << " path(s).";
          }  /* -----  end of function pick_paths  ----- */

        /**
         *  @brief  Find seeds on a set of whole-genome paths for the input reads chunk.
         *
         *  @param[in]  paths_index The index of the set of paths used for finding the seeds.
         *  @param[in]  trav_params Traverse parameters including reads chunk and its index.
         *  @param[in]  callback The call back function applied on the found seeds.
         *
         *  This function uses a set of paths from variation graph to find seeds of the
         *  input set of reads on these paths by traversing the virtual suffix tree of
         *  both indexes of reads chunk and whole-genome paths.
         */
        // :TODO:Mon Mar 06 11:56:\@cartoonist: Function intention and naming is vague.
        inline void
          seeds_on_paths ( Dna5QStringSetIndex < seqan::IndexEsa<> > &paths_index,
              typename TPathTraverser::Param trav_params,
              std::function< void(typename TPathTraverser::Output const &) > callback )
          {
            if ( length ( indexText ( paths_index ) ) == 0 ) return;

            LOG(INFO) << "Finding seeds on paths...";

            // :TODO:Mon Mar 06 13:00:\@cartoonist: IndexEsa<> -> IndexFM<>
            typedef Dna5QStringSetIndex < seqan::IndexEsa<> > TPathIndex;
            typedef typename TPathTraverser::IndexType TReadsIndexSpec;
            typedef Dna5QStringSetIndex < TReadsIndexSpec > TReadsIndex;
            typedef seqan::Seed < seqan::Simple > TSimpleSeed;
            typedef seqan::SeedSet < TSimpleSeed > TSimpleSeedSet;
            typedef seqan::Iterator < TSimpleSeedSet >::Type TSeedIterator;

            TFineIterator < TPathIndex, seqan::ParentLinks<> > paths_itr (paths_index);
            TFineIterator < TReadsIndex, seqan::ParentLinks<> > reads_itr ( trav_params.mutable_get_reads_index() );

            TSimpleSeedSet seeds_set;
            kmer_exact_matches < TPathIndex, TReadsIndex > ( seeds_set, paths_itr, reads_itr,
                trav_params.get_seed_len() );

            // :TODO:Tue Mar 21 10:30:\@cartoonist: Remove this log message.
            LOG(INFO) << "Number of seeds found on paths: " << length ( seeds_set );

            for ( TSeedIterator it = begin ( seeds_set, seqan::Standard() );
                it != end ( seeds_set, seqan::Standard() );
                ++it )
            {
              callback ( *it );
            }

            LOG(INFO) << "Finding seeds on paths: Done.";

          }  /* -----  end of method GraphTraverser::seeds_on_paths  ----- */

        inline void traverse(typename TPathTraverser::Param trav_params,
                      std::function< void(typename TPathTraverser::Output const &) > callback)
        {
          TIMED_SCOPE(traverseTimer, "traverse");
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
              PERFORMANCE_CHECKPOINT_WITH_ID(traverseTimer,
                  std::to_string(nof_reports * TRAVERSE_CHECKPOINT_LOCI_NO));
            }
          }
        }

        inline void add_all_loci(unsigned int step=1,
            std::unordered_set < VarGraph::NodeID > *exclude_nodes=nullptr)
        {
          // TODO: mention in the documentation that the `step` is approximately preserved in
          //       the whole graph. This means that for example add_all_loci(2) would add
          //       the first loci in each node and then add other loci with distance 2 (every
          //       other loci) within the node. So at the end, it won't necessarily preserve
          //       this distance between inter-node loci.
          // **UPDATE** This algorithm use better approximation.
          // TODO: Add documentation.
          TIMED_SCOPE(addAllLociTimer, "add-starts");

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
          seqan::Iterator<VarGraph, BFS<>>::Type itr(this->vargraph);

          unsigned long int prenode_remain = 0;
          unsigned long int remain_estimate = 0;
          VarGraph::NodeID prenode_level = 0;
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
            if ( exclude_nodes == nullptr ||
                exclude_nodes->find(*itr) == exclude_nodes->end() ) {
              bool set = false;
              while (cursor < seq.length())
              {
                vg::Position s_point;
                s_point.set_node_id(*itr);
                s_point.set_offset(cursor);
                this->add_start(s_point);
                set = true;

                cursor += step;
              }

              if ( exclude_nodes != nullptr && exclude_nodes->size() != 0 && !set ) {
                vg::Position s_point;
                s_point.set_node_id(*itr);
                s_point.set_offset(0);
                this->add_start(s_point);
                set = true;
              }
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

          LOG(INFO) << "Number of starting points selected (from "
                    << this->vargraph->nodes_size() << "): "
                    << this->starting_points.size();
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
            std::function< void(typename TPathTraverser::Output const &) > & callback,
            const vg::Position & locus)
        {
          // TODO: Thread unsafe!
          //       Possible solution: non-static variable passed by caller.
          static std::vector< TPathTraverser > path_traversers;
          static std::vector< int > deleted_paths_idx;
          static std::vector< typename TPathTraverser::Output > seeds;
          static std::vector< TPathTraverser > new_ptravs;

#ifndef NDEBUG
          // TODO: Use a structure with atomic increase operation to collect statistics.
          static unsigned long int ptrav_counter = 0;
          static unsigned long int ptrav_len_sum = 0;
          static double      pre_avg_paths_len = 0;
          static double      avg_paths_len = 0;
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
                ptrav_len_sum += ptrav.get_path_length();
                ++ptrav_counter;
                if (ptrav_counter == NOF_PATHLEN_SAMPLES)
                {
                  static int nof_reports = 0;
                  avg_paths_len = ptrav_len_sum / static_cast<float>(NOF_PATHLEN_SAMPLES);
                  if (pre_avg_paths_len != 0)
                  {
                    avg_paths_len = (avg_paths_len + pre_avg_paths_len) / 2.0;
                  }
                  pre_avg_paths_len = avg_paths_len;

                  LOG(INFO) << "Average traversed path length (from "
                            << nof_reports * NOF_PATHLEN_SAMPLES + ptrav_counter
                            << " samples): " << avg_paths_len;
                  ++nof_reports;

                  ptrav_counter = 0;
                  ptrav_len_sum = 0;
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
