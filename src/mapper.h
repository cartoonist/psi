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
#include <type_traits>
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
#include "stat.h"


namespace grem
{
  // :TODO:Mon Sep 25 00:59:\@cartoonist: Fix MapperStat to be thread-safe like TraverserStat.
  /**
   *  @brief  MapperStat template class.
   *
   *  Collect statistics from a `Mapper` class instance(s) in running time.
   */
  template< typename TSpec = void >
    class MapperStat : public Timer
    {
      public:
        /**
         *  @brief  MapperStat constructor.
         *
         *  Start the timer.
         */
        MapperStat( const std::string& name )
          : Timer( name )
        { }

        /**
         *  @brief  Get the current processing locus.
         *
         *  @return The reference to static variable `current_locus`.
         *
         *  It is set to the current processing locus when the mapper is in the middle
         *  of the "traversal" phase.
         */
          static inline vg::Position&
        get_current_locus( )
        {
          static vg::Position current_locus;
          return current_locus;
        }  /* -----  end of method get_current_locus  ----- */

        /**
         *  @brief  Set the current processing locus.
         *
         *  @param  value The value to be set as current locus.
         *
         *  It sets the current processing locus when the mapper is in the middle of
         *  the "traversal" phase.
         */
          static inline void
        set_current_locus( const vg::Position& value )
        {
          get_current_locus() = value;
        }  /* -----  end of method set_current_locus  ----- */

        /**
         *  @brief  Get the index of the current processing locus in the starting loci.
         *
         *  @return The reference to static variable `current_locus_idx`.
         *
         *  It is set to the index of the current processing locus in the starting loci
         *  vector when the mapper is in the middle of the "traversal" phase.
         */
          static inline std::size_t&
        get_current_locus_idx( )
        {
          static std::size_t current_locus_idx;
          return current_locus_idx;
        }  /* -----  end of method get_current_locus_idx  ----- */

        /**
         *  @brief  Set the index of the current processing locus in the starting loci.
         *
         *  @param  value The value to be set as current locus index.
         *
         *  It sets the index of the current processing locus in the starting loci
         *  vector when the mapper is in the middle of the "traversal" phase.
         */
          static inline void
        set_current_locus_idx( const std::size_t& value )
        {
          get_current_locus_idx() = value;
        }  /* -----  end of method set_current_locus_idx  ----- */

        /**
         *  @brief  Get the total number of loci in the starting loci vector.
         *
         *  @return The reference to static variable `total_nof_loci`.
         */
          static inline std::size_t&
        get_total_nof_loci( )
        {
          static std::size_t total_nof_loci;
          return total_nof_loci;
        }  /* -----  end of method get_total_nof_loci  ----- */

        /**
         *  @brief  Set the total number of loci in the starting loci vector.
         *
         *  @param  value The value to be set as total number of loci.
         */
          static inline void
        set_total_nof_loci( const std::size_t& value )
        {
          get_total_nof_loci() = value;
        }  /* -----  end of method set_total_nof_loci  ----- */
    };  /* ----------  end of template class MapperStat  ---------- */

  /**
   *  @brief  MapperStat template specialization for no-stat.
   *
   *  Do nothing.
   */
  template< >
    class MapperStat< NoStat >
    {
      public:
        /* ====================  METHODS       ======================================= */
        MapperStat( const std::string& ) { }
        ~MapperStat() { }
        static inline std::chrono::microseconds get_duration( const std::string& )
        {
          return std::chrono::duration_cast< std::chrono::microseconds >(
              Timer::clock_type::duration::zero() );
        }
        static inline std::chrono::microseconds get_lap( const std::string& )
        {
          return std::chrono::duration_cast< std::chrono::microseconds >(
              Timer::clock_type::duration::zero() );
        }
          static inline vg::Position&
        get_current_locus( )
        {
          static vg::Position current_locus;
          return current_locus;
        }
        static inline void set_current_locus( const vg::Position& value ) { }
          static inline std::size_t&
        get_current_locus_idx( )
        {
          static std::size_t current_locus_idx;
          return current_locus_idx;
        }
        static inline void set_current_locus_idx( const std::size_t& value ) { }
          static inline std::size_t&
        get_total_nof_loci( )
        {
          static std::size_t total_nof_loci;
          return total_nof_loci;
        }
        static inline void set_total_nof_loci( const std::size_t& value ) { }
    };  /* ----------  end of template class MapperStat  ---------- */

  template< class TTraverser, typename TStatSpec = void >
    class Mapper
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename Stat< Mapper >::Type stats_type;
        typedef TTraverser traverser_type;
        /* ====================  LIFECYCLE      ====================================== */
        Mapper( const VarGraph *graph,
            Dna5QRecords&& r,
            unsigned int len,
            unsigned char mismatches = 0 )
          : vargraph( graph ), reads( std::move( r ) ), seed_len( len ),
          seed_mismatches( mismatches )
        {
          if ( length( this->reads.str ) != 0 ) {
            this->index_reads();
          }
        }

        Mapper( const VarGraph *graph,
            const Dna5QRecords& r,
            unsigned int len,
            unsigned char mismatches = 0 )
          : Mapper( graph , Dna5QRecords( r ), len, mismatches )
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
         *
         *  XXX: We assume that each connect component in the graph has one and only one
         *  path indicating a sample haplotype in that region.
         */
        template< typename TIndexSpec >
            void
          pick_paths( PathSet< TIndexSpec >& paths, int n )
          {
            if ( n == 0 ) return;
            auto timer = stats_type( "pick-paths" );

            paths.reserve( n * this->vargraph->path_count );
            for ( std::size_t rank = 1; rank <= this->vargraph->max_path_rank(); ++rank ) {
              const auto& path_name = this->vargraph->path_name( rank );
              auto s = this->vargraph->node_at_path_position( path_name, 0 );
              seqan::Iterator< VarGraph, Haplotyper >::Type hap_itr( this->vargraph, s );
              for ( int i = 0; i < n; ++i ) {
                Path<> new_path( this->vargraph );
                get_uniq_haplotype( new_path, hap_itr );
                paths.add_path( std::move( new_path ) );
              }
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

            auto timer = stats_type( "paths-seed-find" );

            // :TODO:Tue Aug 29 14:48:\@cartoonist: there is a newer `kmer_exact_matches` function!
            //                                      Check `index_iter.h`.
            kmer_exact_matches( paths.index, this->seeds_index, this->seed_len, callback );
          }  /* -----  end of method template Mapper::seeds_on_paths  ----- */

        template< typename TIndexSpec >
            inline void
          add_all_loci( PathSet< TIndexSpec >& paths, unsigned int k,
              unsigned int step=1)
          {
            if ( paths.size() == 0 ) return this->add_all_loci( step );
            auto timer = stats_type( "add-starts" );

            seqan::Iterator< VarGraph, Backtracker >::Type bt_itr ( this->vargraph );
            Path< Full > trav_path( this->vargraph );

            for ( VarGraph::rank_type rank = 1; rank <= this->vargraph->max_node_rank(); ++rank ) {
              VarGraph::nodeid_type id = this->vargraph->rank_to_id( rank );
              auto label_len = this->vargraph->node_length( id );
              std::make_unsigned< VarGraph::offset_type >::type offset = label_len;

              go_begin ( bt_itr, id );
              while ( !at_end( bt_itr ) && offset != 0 ) {
                while ( !at_end( bt_itr ) ) {
                  add_node( trav_path, *bt_itr );
                  // :TODO:Fri Sep 29 00:45:\@cartoonist: in some cases when this loop
                  //     breaks when `bt_itr` is at end but below condition is not met,
                  //     a few unnecessary loci would be added. Bring `add_start` in the
                  //     outer loop and add neccessary loci in the loop; not outside it.
                  if ( trav_path.get_sequence().length() < offset - 1 + k ) ++bt_itr;
                  else break;
                }

                Path< Full > current_path = trav_path;
                while ( !covered_by( current_path.get_nodes(), paths.paths_set ) ) {
                  auto trimmed_len = current_path.get_sequence().length()
                    - this->vargraph->node_length( current_path.get_nodes().back() );
                  if ( trimmed_len <= k - 1 ) {
                    offset = 0;
                    break;
                  }
                  offset = trimmed_len - k + 1;
                  trim( current_path );
                }
                --bt_itr;
                trim( trav_path, *bt_itr );
              }

              for ( auto f = offset; f < label_len; f += step ) {
                this->add_start( id, f );
              }

              clear( trav_path );
            }
          }

        inline void add_all_loci(unsigned int step=1)
        {
          // TODO: Add documentation.
          // TODO: mention in the documentation that the `step` is approximately preserved in
          //       the whole graph.
          auto timer = stats_type( "add-starts" );

          seqan::Iterator<VarGraph, BFS>::Type itr(this->vargraph);

          unsigned long int prenode_remain = 0;
          unsigned long int remain_estimate = 0;
          VarGraph::nodeid_type prenode_level = 0;
          while (!at_end(itr)) {
            if (prenode_level != level(itr)) {
              prenode_remain = remain_estimate;
              remain_estimate = 0;
              prenode_level = level(itr);
            }

            auto seq_len = this->vargraph->node_length(*itr);
            unsigned long int cursor = (step - prenode_remain) % step;
            while (cursor < seq_len) {
              this->add_start(*itr, cursor);
              cursor += step;
            }

            unsigned long int new_remain;
            if (step - prenode_remain > seq_len) {
              new_remain = prenode_remain + seq_len;
            }
            else {
              new_remain = (seq_len - step + prenode_remain) % step;
            }

            if (remain_estimate < new_remain) {
              remain_estimate = new_remain;
            }

            ++itr;
          }
        }
        /* ====================  METHODS       ======================================= */
          inline void
        traverse( std::function< void( typename TTraverser::output_type const& ) >& callback )
        {
          auto timer = stats_type( "traverse" );
          stats_type::set_total_nof_loci( this->starting_loci.size() );

          TTraverser traverser( this->vargraph, &(this->reads_index), this->seed_len );
          for ( std::size_t idx = 0; idx < this->starting_loci.size(); ++idx )
          {
            const auto& locus = this->starting_loci[ idx ];
            stats_type::set_current_locus( locus );
            stats_type::set_current_locus_idx( idx );

            traverser.set_start_locus( locus );
            traverser.run( callback );
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
          {
            auto timer = stats_type( "index-reads" );
            this->reads_index = typename TTraverser::index_type( this->reads.str );
          }
          {
            auto timer = stats_type( "seeding" );
            this->seeds =
              seeding ( this->reads.str, this->seed_len, FixedLengthNonOverlapping() );
          }
          {
            auto timer = stats_type( "index-seeds" );
            this->seeds_index = typename TTraverser::index_type( this->seeds );
          }
        }
    };

  /**
   *  @brief  Stat template class specialization for `MapperStat`.
   */
  template< class TTraverser, typename TSpec >
    class Stat< Mapper< TTraverser, TSpec > >
    {
      public:
        typedef MapperStat< TSpec > Type;
    };  /* ----------  end of template class Stat  ---------- */
}

#endif  // end of MAPPER_H__
