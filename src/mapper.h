/**
 *    @file  mapper.h
 *   @brief  Mapper template class.
 *
 *  Mapper template class definition and its template specialisations.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Aug 23, 2017  12:44
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef MAPPER_H__
#define MAPPER_H__

#include <fstream>
#include <vector>
#include <iterator>
#include <functional>
#include <algorithm>

#include "vargraph.h"
#include "sequence.h"
#include "index.h"
#include "index_iter.h"
#include "utils.h"
#include "logger.h"

namespace grem
{
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

        inline void add_start ( const VarGraph::nodeid_type & node_id,
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
            std::vector < VarGraph::nodeid_type > new_path;
            seqan::Dna5QString new_path_str;
            VarGraph::NodeCoverage covered_nodes;

            new_path.reserve ( this->vargraph->node_count );
            covered_nodes.reserve ( this->vargraph->node_count );
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
          std::vector< VarGraph::nodeid_type > trav_path;
          unsigned int trav_len = 0;

          // :TODO:Sun Jun 11 21:36:\@cartoonist: traverse the graph using BFS instead
          //   of iterating over node list would be more cache oblivious.
          for ( VarGraph::rank_type rank = 1; rank <= this->vargraph->max_node_rank(); ++rank ) {
            VarGraph::nodeid_type start_node_id = this->vargraph->rank_to_id( rank );
            unsigned int label_len = this->vargraph->node_sequence( start_node_id ).length();

            bool set = false;
            unsigned int init_offset = ( label_len < k - 1 ) ? 0 : label_len - k + 1;
            for ( unsigned int offset = init_offset; offset < label_len; offset += step ) {
              // :TODO:Mon May 22 14:40:\@cartoonist: missed some locations when the
              //     the length of branch node's label is less than k.
              if ( ! this->vargraph->is_branch ( start_node_id ) &&
                  covered_by ( start_node_id, paths_coverage ) &&
                  this->vargraph->has_edges_from ( start_node_id ) &&
                    this->vargraph->node_sequence(
                      this->vargraph->edges_from ( start_node_id ).at(0).to() )
                      .length() > k ) {
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
                    trav_len += this->vargraph->node_sequence( *bt_itr ).length();
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

                VarGraph::nodeid_type poped_id = 0;
                while ( !trav_path.empty() && poped_id != *bt_itr ) {
                  poped_id = trav_path.back();
                  trav_len -= this->vargraph->node_sequence( poped_id ).length();
                  trav_path.pop_back();
                }
              }
            }
          }

          LOG(INFO) << "Number of starting points selected (from "
                    << this->vargraph->node_count << "): "
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
          VarGraph::nodeid_type prenode_level = 0;
          std::string seq;
          while (!at_end(itr))
          {
            if (prenode_level != level(itr))
            {
              prenode_remain = remain_estimate;
              remain_estimate = 0;
              prenode_level = level(itr);
            }

            seq = this->vargraph->node_sequence(*itr);

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
                    << this->vargraph->node_count << "): "
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
    save_paths_coverage ( std::vector< VarGraph::NodeCoverage >& paths_covered_nodes,
        const std::string path_prefix )
    {
      std::ofstream file_stream;
      unsigned int i = 0;
      for ( const auto& covered_nodes : paths_covered_nodes ) {
        std::string file_path = path_prefix + std::to_string ( i );
        file_stream.open( file_path, std::ofstream::out | std::ofstream::binary );
        serialize( file_stream, covered_nodes,
            covered_nodes.begin(), covered_nodes.end() );
        file_stream.close();
        ++i;
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
    load_paths_coverage ( std::vector< VarGraph::NodeCoverage >& paths_covered_nodes,
        const std::string path_prefix, unsigned int path_num )
    {
      std::ifstream file_stream;

      VarGraph::NodeCoverage covered_nodes;
      paths_covered_nodes.reserve ( path_num );

      for ( unsigned int i = 0; i < path_num; ++i ) {
        std::string file_path = path_prefix + std::to_string ( i );
        file_stream.open( file_path, std::ifstream::in | std::ifstream::binary );

        if ( !file_stream ) {
          paths_covered_nodes.clear();
          return false;
        }

        deserialize( file_stream, covered_nodes,
            std::inserter( covered_nodes, covered_nodes.end() ));
        paths_covered_nodes.push_back ( covered_nodes );
        file_stream.close();
        covered_nodes.clear();
      }

      return true;
    }  /* -----  end of function load_paths_coverage  ----- */
}

#endif  // end of MAPPER_H__
