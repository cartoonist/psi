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
#include <fstream>
#include <memory>

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

  template <class TPathTraverser>
    class Mapper
    {
      public:
        // Constructors
        Mapper(const VarGraph *graph,
            const std::vector< vg::Position > *start_loci) :
          vargraph(graph)
        {
          if (start_loci != nullptr) this->starting_points = *start_loci;
        }

        Mapper(const VarGraph &graph,
            const std::vector< vg::Position > *start_loci) :
          Mapper(&graph, start_loci)
        {}

        Mapper(const VarGraph &graph,
            const std::vector< vg::Position > &start_loci) :
          Mapper(&graph, &start_loci)
        {}

        Mapper(const VarGraph *graph,
            const std::vector< vg::Position > &start_loci) :
          Mapper(graph, &start_loci)
        {}

        Mapper(const VarGraph *graph) : Mapper(graph, nullptr)
        {}

        Mapper(const VarGraph &graph) : Mapper(&graph)
        {}

        // Public methods
        inline void add_start(const vg::Position &locus)
        {
          this->starting_points.push_back(locus);
        }

        inline void add_start ( const VarGraph::NodeID & node_id,
            const unsigned long int & offset )
        {
          vg::Position locus;
          locus.set_node_id(node_id);
          locus.set_offset(offset);
          this->add_start ( locus );
        }

        /**
         *  @brief  Pick n paths from the variation graph.
         *
         *  @param[out]  paths The set of generated paths are added to this string set.
         *  @param[out]  paths_covered_nodes The vector of `NodeCoverage` of each path.
         *  @param[in]  n Number of paths.
         *
         *  This method generates a set of (probably) unique whole-genome path from the
         *  variation graph.
         */
        void
          pick_paths ( Dna5QStringSet &paths,
              std::vector< VarGraph::NodeCoverage > &paths_covered_nodes, int n )
          {
            if ( n == 0 ) return;

            TIMED_SCOPE(pickPathsTimer, "pick-paths");

            LOG(INFO) << "Picking " << n << " different path(s) on the graph...";

            seqan::Iterator < VarGraph, Haplotyper<> >::Type hap_itr ( this->vargraph );
            std::vector < VarGraph::NodeID > new_path;
            seqan::Dna5QString new_path_str;
            VarGraph::NodeCoverage covered_nodes;

            new_path.reserve ( this->vargraph->nodes_size() );
            covered_nodes.reserve ( this->vargraph->nodes_size() );
            paths_covered_nodes.reserve ( n );

            for ( int i = 0; i < n; ++i ) {
              get_uniq_haplotype ( new_path, hap_itr );

              std::copy ( new_path.begin(), new_path.end(),
                  std::inserter ( covered_nodes, covered_nodes.end() ) );
              paths_covered_nodes.push_back ( covered_nodes );

              new_path_str = this->vargraph->get_string ( new_path );

              // :TODO:Mon Mar 06 13:00:\@cartoonist: faked quality score.
              char fake_qual = 'I';
              assignQualities ( new_path_str, std::string ( length(new_path_str), fake_qual ) );

              appendValue ( paths, new_path_str );

              new_path.clear();
              covered_nodes.clear();
            }
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

            TIMED_SCOPE(pathsSeedFindTimer, "paths-seed-find");

            LOG(INFO) << "Finding seeds on paths...";

            // :TODO:Mon Mar 06 13:00:\@cartoonist: IndexEsa<> -> IndexFM<>
            //typedef Dna5QStringSetIndex < seqan::IndexEsa<> > TPathIndex;
            //typedef Dna5QStringSetIndex < typename TPathTraverser::IndexType > TReadsIndex;

            //TFineIndexIter < TPathIndex, seqan::ParentLinks<> > paths_itr (paths_index);
            //TFineIndexIter < TReadsIndex, seqan::ParentLinks<> > reads_itr ( trav_params.mutable_get_reads_index() );
            //seqan::Iterator < TPathIndex, seqan::TopDown<seqan::ParentLinks<>> >::Type paths_itr (paths_index);
            //typename seqan::Iterator < TReadsIndex, seqan::TopDown<seqan::ParentLinks<>> >::Type reads_itr ( trav_params.mutable_get_reads_index() );

            //std::vector < seqan::Seed < seqan::Simple > > seeds_set;
            //kmer_exact_matches < TReadsIndex, TPathIndex > ( seeds_set, reads_itr, paths_itr, trav_params.get_seed_len() );
            kmer_exact_matches ( paths_index, trav_params.mutable_get_seeds_index(),
               trav_params.get_seed_len(), callback );

            //std::for_each ( seeds_set.begin(), seeds_set.end(), callback );
          }  /* -----  end of method Mapper::seeds_on_paths  ----- */

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

        inline void add_all_loci(std::vector < VarGraph::NodeCoverage > &paths_coverage,
            unsigned int k, unsigned int step=1)
        {
          if ( paths_coverage.size() == 0 ) return this->add_all_loci(step);

          seqan::Iterator< VarGraph, Backtracker<> >::Type bt_itr ( this->vargraph );
          std::vector< VarGraph::NodeID > trav_path;
          unsigned int trav_len = 0;

          // :TODO:Sun Jun 11 21:36:\@cartoonist: traverse the graph using BFS instead
          //   of iterating over node list would be more cache oblivious.
          for ( unsigned long int idx = 0; idx < this->vargraph->nodes_size(); ++idx ) {
            const VarGraph::Node &start_node = this->vargraph->node_at ( idx );
            VarGraph::NodeID start_node_id = start_node.id();
            unsigned int label_len = start_node.sequence().length();

            bool set = false;
            unsigned int init_offset = ( label_len < k - 1 ) ? 0 : label_len - k + 1;
            for ( unsigned int offset = init_offset; offset < label_len; offset += step ) {
              // :TODO:Mon May 22 14:40:\@cartoonist: missed some locations when the
              //     the length of branch node's label is less than k.
              if ( ! this->vargraph->is_branch ( start_node_id ) &&
                  covered_by ( start_node_id, paths_coverage ) &&
                  this->vargraph->has_fwd_edge ( start_node_id ) &&
                    this->vargraph->node_by (
                      this->vargraph->fwd_edges ( start_node_id ).at(0)->to() )
                      .sequence().length() > k ) {
                  continue;
              }

              if ( set ) {
                this->add_start( start_node_id, offset );
                continue;
              }

              go_begin ( bt_itr, start_node_id );

              while ( !at_end( bt_itr ) ) {
                while ( !at_end( bt_itr ) ) {
                  trav_path.push_back ( *bt_itr );
                  if ( *bt_itr != start_node_id ) {
                    trav_len += this->vargraph->node_by ( *bt_itr ).sequence().length();
                  }
                  else {
                    trav_len = label_len - offset;
                  }

                  if ( trav_len < k ) ++bt_itr;
                  else break;
                }

                if ( ! covered_by ( trav_path, paths_coverage ) ) {
                  this->add_start( start_node_id, offset );
                  set = true;
                  break;
                }

                --bt_itr;

                VarGraph::NodeID poped_id = 0;
                while ( !trav_path.empty() && poped_id != *bt_itr ) {
                  poped_id = trav_path.back();
                  trav_len -= this->vargraph->node_by ( poped_id ).sequence().length();
                  trav_path.pop_back();
                }
              }
            }
          }

          LOG(INFO) << "Number of starting points selected (from "
                    << this->vargraph->nodes_size() << "): "
                    << this->starting_points.size();
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
          TIMED_SCOPE(addAllLociTimer, "add-starts");

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
            while (cursor < seq.length())
            {
              this->add_start(*itr, cursor);
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

  /* interface functions */

  /**
   *  @brief  Save coverage (nodes covered by each path) of the picked paths.
   *
   *  @param[in]  paths_covered_nodes The vector of `NodeCoverage` of each path.
   *  @param[in]  paths_index_file The file path prefix of the paths index.
   *
   *  It saves the paths coverage; i.e. the list of node IDs covered by each path,
   *  separatedly in files. These files in conjuction with the paths index can be
   *  considered as an offline index for seed finding.
   */
  void
    save_paths_coverage ( std::vector< VarGraph::NodeCoverage > &paths_covered_nodes,
        const std::string &path_prefix )
    {
      std::ofstream file_stream;
      for ( unsigned int i = 0; i < paths_covered_nodes.size(); ++i ) {
        const auto &covered_nodes = paths_covered_nodes[i];
        file_stream.open ( path_prefix + "_path_" + std::to_string ( i ),
            std::ofstream::out | std::ofstream::binary );

        uint64_t set_size = covered_nodes.size();
        file_stream.write
          ( reinterpret_cast<char*>( &set_size ), sizeof ( uint64_t ) );
        for ( const auto &node_id : covered_nodes ) {
          file_stream.write
            ( reinterpret_cast<const char*>( &node_id ), sizeof ( VarGraph::NodeID ) );
        }

        file_stream.close();
      }
    }  /* -----  end of function save_paths_coverage  ----- */

  /**
   *  @brief  Load coverage (node covered by each path) from file.
   *
   *  @param[out]  paths_covered_nodes The vector of `NodeCoverage` to be filled.
   *  @param[in]   path_prefix The file path prefix of the saved paths.
   *  @return true if it successfully load the node coverage of the paths from file;
   *          otherwise false.
   *
   *  It reads the saved files and insert each node ID to the unordered set. The file is
   *  prefixed by the number of node IDs stored in the file.
   */
  bool
    load_paths_coverage ( std::vector< VarGraph::NodeCoverage > &paths_covered_nodes,
        const std::string &path_prefix, unsigned int path_num )
    {
      std::ifstream file_stream;
      VarGraph::NodeCoverage covered_nodes;
      VarGraph::NodeID node_id;

      paths_covered_nodes.reserve ( path_num );

      for ( unsigned int i = 0; i < path_num; ++i ) {
        file_stream.open ( path_prefix + "_path_" + std::to_string ( i ),
            std::ifstream::in | std::ifstream::binary );
        if ( !file_stream ) {
          paths_covered_nodes.clear();
          return false;
        }

        uint64_t set_size;
        file_stream.read ( reinterpret_cast<char*>( &set_size ), sizeof ( uint64_t ) );
        covered_nodes.reserve ( set_size );
        for ( unsigned int j = 0; j < set_size; ++j ) {
          file_stream.read
            ( reinterpret_cast<char*>( &node_id ), sizeof ( VarGraph::NodeID ) );
          covered_nodes.insert ( node_id );
        }
        paths_covered_nodes.push_back ( covered_nodes );

        file_stream.close();
        covered_nodes.clear();
      }

      return true;
    }  /* -----  end of function load_paths_coverage  ----- */
}

#endif
