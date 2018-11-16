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
#include "pathset.h"
#include "utils.h"
#include "logger.h"

namespace grem
{
  template <class TTraverser>
    class Mapper
    {
      public:
        /* ====================  LIFECYCLE      ====================================== */
        Mapper( const VarGraph *graph,
            Dna5QRecords&& reads,
            unsigned int len,
            unsigned char mismatches = 0 )
          : vargraph( graph ), reads( std::move( reads ) ), seed_len( len ),
          seed_mismatches( mismatches )
        {
          if ( length( reads.str ) != 0 ) {
            this->index_reads();
          }
        }

        Mapper( const VarGraph *graph,
            const Dna5QRecords& reads,
            unsigned int len,
            unsigned char mismatches = 0 )
          : Mapper( graph , Dna5QRecords( reads ), len, mismatches )
        { }

        Mapper( const VarGraph *graph,
            unsigned int len,
            unsigned char mismatches = 0 )
          : Mapper( graph , Dna5QRecords( ), len, mismatches )
        { }
        /* ====================  ACCESSORS      ====================================== */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const VarGraph*
        get_vargraph ( ) const
        {
          return this->vargraph;
        }  /* -----  end of method get_vargraph  ----- */

        /**
         *  @brief  getter function for starting_loci.
         */
          inline const std::vector< vg::Position >&
        get_starting_loci( ) const
        {
          return this->starting_loci;
        }  /* -----  end of method get_starting_loci  ----- */

        /**
         *  @brief  getter function for seed_len.
         */
          inline unsigned int
        get_seed_len ( ) const
        {
          return this->seed_len;
        }  /* -----  end of method get_seed_len  ----- */

        /**
         *  @brief  getter function for seed_mismatches.
         */
          inline unsigned char
        get_seed_mismatches ( ) const
        {
          return this->seed_mismatches;
        }  /* -----  end of method get_seed_mismatches  ----- */

        /**
         *  @brief  getter function for reads.
         */
          inline const Dna5QRecords&
        get_reads ( ) const
        {
          return this->reads;
        }  /* -----  end of method get_reads  ----- */
        /* ====================  MUTATORS       ====================================== */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph ( const VarGraph* value )
        {
          this->vargraph = value;
        }  /* -----  end of method set_vargraph  ----- */

        /**
         *  @brief  setter function for starting_loci.
         *
         *  Copy assignment.
         */
          inline void
        set_starting_loci( const std::vector< vg::Position >& loci )
        {
          this->starting_loci = loci;
        }  /* -----  end of method set_starting_loci  ----- */

        /**
         *  @brief  setter function for starting_loci.
         *
         *  Move assignment.
         */
          inline void
        set_starting_loci( std::vector< vg::Position >&& loci )
        {
          this->starting_loci = std::move( loci );
        }  /* -----  end of method set_starting_loci  ----- */

        /**
         *  @brief  setter function for seed_len.
         */
          inline void
        set_seed_len ( unsigned int value )
        {
          this->seed_len = value;
        }  /* -----  end of method set_seed_len  ----- */

        /**
         *  @brief  setter function for seed_mismatches.
         */
          inline void
        set_seed_mismatches ( unsigned char value )
        {
          this->seed_mismatches = value;
        }  /* -----  end of method set_seed_mismatches  ----- */

        /**
         *  @brief  setter function for reads.
         *
         *  Move assignment.
         */
          inline void
        set_reads ( Dna5QRecords&& value )
        {
          this->reads = std::move( value );
          this->index_reads();
        }  /* -----  end of method set_reads  ----- */

        /**
         *  @brief  setter function for reads.
         *
         *  Copy assignment.
         */
          inline void
        set_reads ( const Dna5QRecords& value )
        {
          this->set_reads( Dna5QRecords( value ) );
        }  /* -----  end of method set_reads  ----- */

          inline void
        add_start( const vg::Position &locus )
        {
          this->starting_loci.push_back( locus );
        }

          inline void
        add_start( VarGraph::nodeid_type node_id, VarGraph::offset_type offset )
        {
          vg::Position locus;
          locus.set_node_id( node_id );
          locus.set_offset( offset );
          this->add_start ( locus );
        }
        /* ====================  METHODS        ====================================== */
        /**
         *  @brief  Pick n paths from the variation graph.
         *
         *  @param[out]  paths Set of paths to be generated.
         *  @param[in]  n Number of paths.
         *
         *  This method generates a set of (probably) unique whole-genome paths from the
         *  variation graph.
         */
        template< typename TIndexSpec >
            void
          pick_paths( PathSet< TIndexSpec >& paths, int n )
          {
            if ( n == 0 ) return;

            LOG(INFO) << "Picking " << n << " different path(s) on the graph...";
            TIMED_SCOPE(pickPathsTimer, "pick-paths");

            seqan::Iterator< VarGraph, Haplotyper >::Type hap_itr( this->vargraph );

            paths.reserve( n );
            for ( int i = 0; i < n; ++i ) {
              Path<> new_path( this->vargraph );
              reserve( new_path, this->vargraph->node_count );
              get_uniq_haplotype( new_path, hap_itr );
              paths.add_path( std::move( new_path ) );
            }
          }  /* -----  end of template function pick_paths  ----- */

        /**
         *  @brief  Find seeds on a set of whole-genome paths for the input reads chunk.
         *
         *  @param[in]  paths The set of paths used for finding the seeds.
         *  @param[in]  callback The call back function applied on the found seeds.
         *
         *  This function uses a set of paths from variation graph to find seeds of the
         *  input set of reads on these paths by traversing the virtual suffix tree of
         *  both indexes of reads chunk and whole-genome paths.
         */
        // :TODO:Mon Mar 06 11:56:\@cartoonist: Function intention and naming is vague.
        template< typename TIndexSpec >
            inline void
          seeds_on_paths( PathSet< TIndexSpec >& paths,
              std::function< void(typename TTraverser::output_type const &) >& callback )
          {
            if ( length( paths.string_set ) == 0 ) return;

            LOG(INFO) << "Finding seeds on paths...";
            TIMED_SCOPE(pathsSeedFindTimer, "paths-seed-find");

            // :TODO:Tue Aug 29 14:48:\@cartoonist: there is a newer `kmer_exact_matches` function!
            //                                      Check `index_iter.h`.
            kmer_exact_matches( paths.index, this->seeds_index, this->seed_len, callback );
          }  /* -----  end of method template Mapper::seeds_on_paths  ----- */

        template< typename TIndexSpec >
            inline void
          add_all_loci( PathSet< TIndexSpec >& paths, unsigned int seed_len,
              unsigned int step=1)
          {
            if ( paths.size() == 0 ) return this->add_all_loci( step );

            seqan::Iterator< VarGraph, Backtracker >::Type bt_itr ( this->vargraph );
            std::vector< VarGraph::nodeid_type > trav_path;
            unsigned int trav_len = 0;

            // :TODO:Sun Jun 11 21:36:\@cartoonist: traverse the graph using BFS instead
            //   of iterating over node list would be more cache oblivious.
            for ( VarGraph::rank_type rank = 1; rank <= this->vargraph->max_node_rank(); ++rank ) {
              VarGraph::nodeid_type start_node_id = this->vargraph->rank_to_id( rank );
              unsigned int label_len = this->vargraph->node_sequence( start_node_id ).length();

              bool set = false;
              unsigned int init_offset = ( label_len < seed_len - 1 ) ? 0 : label_len - seed_len + 1;
              for ( unsigned int offset = init_offset; offset < label_len; offset += step ) {
                // :TODO:Mon May 22 14:40:\@cartoonist: missed some locations when the
                //     the length of branch node's label is less than seed_len.
                if ( ! this->vargraph->is_branch ( start_node_id ) &&
                    covered_by ( start_node_id, paths.paths_set ) &&
                    this->vargraph->has_edges_from ( start_node_id ) &&
                      this->vargraph->node_sequence(
                        this->vargraph->edges_from ( start_node_id ).at(0).to() )
                        .length() > seed_len ) {
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

                    if ( trav_len < seed_len ) ++bt_itr;
                    else break;
                  }

                  if ( ! covered_by ( trav_path, paths.paths_set ) ) {
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
                      << this->starting_loci.size();
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

          seqan::Iterator<VarGraph, BFS>::Type itr(this->vargraph);

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
                    << this->starting_loci.size();
        }
        /* ====================  METHODS       ======================================= */
          inline void
        traverse ( std::function< void( typename TTraverser::output_type const& ) >& callback )
        {
          TIMED_SCOPE(traverseTimer, "traverse");
#ifndef NDEBUG
          unsigned int locus_counter = 0;
          unsigned int nof_reports = 0;
#endif
          TTraverser traverser( this->vargraph, &(this->reads_index), this->seed_len );
          for (auto locus : this->starting_loci)
          {
            traverser.set_start_locus( locus );
            traverser.run( callback );
#ifndef NDEBUG
            ++locus_counter;
            if (locus_counter % TRAVERSE_CHECKPOINT_LOCI_NO == 0)
            {
              locus_counter = 0;
              ++nof_reports;
              PERFORMANCE_CHECKPOINT_WITH_ID(traverseTimer,
                  std::to_string(nof_reports * TRAVERSE_CHECKPOINT_LOCI_NO));
            }
#endif
          }
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;
        std::vector< vg::Position > starting_loci;
        Dna5QRecords reads;
        unsigned int seed_len;
        unsigned char seed_mismatches;  /**< @brief Allowed mismatches in a seed hit. */
        typename TTraverser::index_type reads_index;
        Dna5QStringSet seeds;
        typename TTraverser::index_type seeds_index;
        /* ====================  METHODS       ======================================= */
          inline void
        index_reads( )
        {
          TIMED_BLOCK(readsIndexTimer, "index-reads")
          {
            this->reads_index = typename TTraverser::index_type( this->reads.str );
          }
          TIMED_BLOCK(seedingTimer, "seeding") {
            this->seeds =
              seeding ( this->reads.str, this->seed_len, FixedLengthNonOverlapping() );
          }
          TIMED_BLOCK(seedsIndexTimer, "index-seeds") {
            this->seeds_index = typename TTraverser::index_type( this->seeds );
          }
        }
    };
}

#endif  // end of MAPPER_H__
