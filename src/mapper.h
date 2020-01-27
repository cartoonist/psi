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
#include <atomic>
#include <vector>
#include <iterator>
#include <functional>
#include <algorithm>

#include <sdsl/bit_vectors.hpp>

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
  /**
   *  @brief  MapperStat template class.
   *
   *  Collect statistics from a `Mapper` class instance(s) in running time.
   */
  template< typename TSpec = void >
    class MapperStat : public Timer
    {
      public:
        /* ====================  MEMBER TYPES  ======================================= */
        struct Coordinates {
          VarGraph::nodeid_type node_id;
          VarGraph::offset_type offset;
        };
        /* ====================  LIFECYCLE     ======================================= */
        /**
         *  @brief  MapperStat constructor.
         *
         *  Start the timer.
         */
        MapperStat( const std::string& name )
          : Timer( name )
        { }
        /* ====================  ACCESSORS     ======================================= */
        /**
         *  @brief  Get the last processing locus.
         *
         *  @return The reference to static variable `lastproc_locus`.
         *
         *  It is set to the last processing locus when the mapper is in the middle
         *  of the "traversal" phase.
         */
          static inline std::atomic< Coordinates >&
        get_lastproc_locus( )
        {
          static std::atomic< Coordinates > lastproc_locus( { 0, 0 } );
          return lastproc_locus;
        }  /* -----  end of method get_lastproc_locus  ----- */

        /**
         *  @brief  Get the index of the last done processing locus in the starting loci.
         *
         *  @return The reference to static variable `lastdone_locus_idx`.
         *
         *  It is set to the index of the last done processing locus in the starting loci
         *  vector when the mapper is in the middle of the "traversal" phase.
         */
          static inline std::atomic< std::size_t >&
        get_lastdone_locus_idx( )
        {
          static std::atomic< std::size_t > lastdone_locus_idx( 0 );
          return lastdone_locus_idx;
        }  /* -----  end of method get_lastdone_locus_idx  ----- */

        /**
         *  @brief  Get the total number of loci in the starting loci vector.
         *
         *  @return The reference to static variable `total_nof_loci`.
         */
          static inline std::size_t&
        get_total_nof_loci( )
        {
          static std::size_t total_nof_loci = 0;
          return total_nof_loci;
        }  /* -----  end of method get_total_nof_loci  ----- */
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  Set the last processing locus.
         *
         *  @param  value The value to be set as last processing locus.
         *
         *  It sets the last processing locus when the mapper is in the middle of
         *  the "traversal" phase.
         */
          static inline void
        set_lastproc_locus( const vg::Position& value )
        {
          get_lastproc_locus().store( { value.node_id(), value.offset() } );
        }  /* -----  end of method set_lastproc_locus  ----- */

        /**
         *  @brief  Set the index of the last done processing locus in the starting loci.
         *
         *  @param  value The value to be set as last done locus index.
         *
         *  It sets the index of the last done processing locus in the starting loci
         *  vector when the mapper is in the middle of the "traversal" phase.
         */
          static inline void
        set_lastdone_locus_idx( const std::size_t& value )
        {
          get_lastdone_locus_idx().store( value );
        }  /* -----  end of method set_lastdone_locus_idx  ----- */

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
        /* ====================  MEMBER TYPES  ======================================= */
        struct Coordinates {
          VarGraph::nodeid_type node_id;
          VarGraph::offset_type offset;
        };
        /* ====================  LIFECYCLE     ======================================= */
        MapperStat( const std::string& ) { }
        ~MapperStat() { }
        /* ====================  METHODS       ======================================= */
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
        /* ====================  ACCESSORS     ======================================= */
          static inline std::atomic< Coordinates >&
        get_lastproc_locus( )
        {
          static std::atomic< Coordinates > lastproc_locus( { 0, 0 } );
          return lastproc_locus;
        }
          static inline std::atomic< std::size_t >&
        get_lastdone_locus_idx( )
        {
          static std::atomic< std::size_t > lastdone_locus_idx( 0 );
          return lastdone_locus_idx;
        }
          static inline std::size_t&
        get_total_nof_loci( )
        {
          static std::size_t total_nof_loci = 0;
          return total_nof_loci;
        }
        /* ====================  METHODS       ======================================= */
        static inline void set_lastproc_locus( const vg::Position& value ) { }
        static inline void set_lastdone_locus_idx( const std::size_t& value ) { }
        static inline void set_total_nof_loci( const std::size_t& value ) { }
    };  /* ----------  end of template class MapperStat  ---------- */

  template< class TTraverser, typename TStatSpec = void >
    class Mapper
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef typename Stat< Mapper >::Type stats_type;
        typedef TTraverser traverser_type;
        typedef Records< typename TTraverser::stringset_type > readsrecord_type;
        typedef typename TTraverser::index_type readsindex_type;
        /* ====================  LIFECYCLE      ====================================== */
        Mapper( const VarGraph* graph,
            readsrecord_type&& r,
            unsigned int len,
            unsigned char mismatches = 0 )
          : vargraph( graph ), reads( std::move( r ) ), seed_len( len ),
          seed_mismatches( mismatches )
        {
          if ( length( this->reads.str ) != 0 ) {
            this->index_reads();
          }
        }

        Mapper( const VarGraph* graph,
            const readsrecord_type& r,
            unsigned int len,
            unsigned char mismatches = 0 )
          : Mapper( graph , readsrecord_type( r ), len, mismatches )
        { }

        Mapper( const VarGraph* graph,
            unsigned int len,
            unsigned char mismatches = 0 )
          : Mapper( graph , readsrecord_type( ), len, mismatches )
        { }
        /* ====================  ACCESSORS      ====================================== */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const VarGraph*
        get_vargraph( ) const
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
        get_seed_len( ) const
        {
          return this->seed_len;
        }  /* -----  end of method get_seed_len  ----- */

        /**
         *  @brief  getter function for seed_mismatches.
         */
          inline unsigned char
        get_seed_mismatches( ) const
        {
          return this->seed_mismatches;
        }  /* -----  end of method get_seed_mismatches  ----- */

        /**
         *  @brief  getter function for reads.
         */
          inline const readsrecord_type&
        get_reads( ) const
        {
          return this->reads;
        }  /* -----  end of method get_reads  ----- */
        /* ====================  MUTATORS       ====================================== */
        /**
         *  @brief  setter function for vargraph.
         */
          inline void
        set_vargraph( const VarGraph* value )
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
        set_seed_len( unsigned int value )
        {
          this->seed_len = value;
        }  /* -----  end of method set_seed_len  ----- */

        /**
         *  @brief  setter function for seed_mismatches.
         */
          inline void
        set_seed_mismatches( unsigned char value )
        {
          this->seed_mismatches = value;
        }  /* -----  end of method set_seed_mismatches  ----- */

        /**
         *  @brief  setter function for reads.
         *
         *  Move assignment.
         */
          inline void
        set_reads( readsrecord_type&& value )
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
        set_reads( const readsrecord_type& value )
        {
          this->set_reads( readsrecord_type( value ) );
        }  /* -----  end of method set_reads  ----- */

          inline void
        add_start( const vg::Position& locus )
        {
          this->starting_loci.push_back( locus );
        }

          inline void
        add_start( VarGraph::nodeid_type node_id, VarGraph::offset_type offset )
        {
          vg::Position locus;
          locus.set_node_id( node_id );
          locus.set_offset( offset );
          this->add_start( locus );
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
        template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection >
            void
          pick_paths( PathSet< TGraph, TText, TIndexSpec, TSequenceDirection >& paths,
              int n, bool patched=true )
          {
            if ( n == 0 ) return;
            auto timer = stats_type( "pick-paths" );

            paths.reserve( n * this->vargraph->path_count );
            for ( std::size_t rank = 1; rank <= this->vargraph->max_path_rank(); ++rank ) {
              const auto& path_name = this->vargraph->path_name( rank );
              auto s = this->vargraph->node_at_path_position( path_name, 0 );
              seqan::Iterator< VarGraph, Haplotyper >::Type hap_itr( this->vargraph, s );
              for ( int i = 0; i < n; ++i ) {
                if ( patched ) get_uniq_patched_haplotype( paths, hap_itr, this->seed_len );
                else get_uniq_full_haplotype( paths, hap_itr );
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
        template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection >
            inline void
          seeds_on_paths( PathSet< TGraph, TText, TIndexSpec, TSequenceDirection >& paths,
              std::function< void(typename TTraverser::output_type const &) >& callback )
          {
            if ( length( indexText( paths.index ) ) == 0 ) return;

            auto timer = stats_type( "paths-seed-find" );

            kmer_exact_matches( paths.index, &paths, &(this->reads), this->seed_len, GreedyOverlapping(), callback );
          }  /* -----  end of method template Mapper::seeds_on_paths  ----- */

        template< typename TGraph, typename TText, typename TIndexSpec, typename TSequenceDirection >
            inline void
          add_all_loci( PathSet< TGraph, TText, TIndexSpec, TSequenceDirection >& paths, unsigned int k,
              unsigned int step=1)
          {
            if ( paths.size() == 0 ) return this->add_all_loci( step );
            auto timer = stats_type( "add-starts" );

            seqan::Iterator< VarGraph, Backtracker >::Type bt_itr( this->vargraph );
            Path< VarGraph > trav_path( this->vargraph );
            Path< VarGraph > current_path( this->vargraph );
            sdsl::bit_vector bv_starts( this->vargraph->get_max_node_len(), 0 );

            for ( VarGraph::rank_type rank = 1; rank <= this->vargraph->max_node_rank(); ++rank ) {
              VarGraph::nodeid_type id = this->vargraph->rank_to_id( rank );
              auto label_len = this->vargraph->node_length( id );
              std::make_unsigned< VarGraph::offset_type >::type offset = label_len;

              go_begin( bt_itr, id );
              while ( !at_end( bt_itr ) && offset != 0 ) {
                extend_to_k( trav_path, bt_itr, offset - 1 + k );
                if ( trav_path.get_sequence_len() >= k ) current_path = trav_path;
                while ( current_path.get_sequence_len() != 0 &&
                    !covered_by( current_path, paths, Unordered() ) ) {
                  auto trimmed_len = current_path.get_sequence_len()
                    - this->vargraph->node_length( current_path.get_nodes().back() );
                  if ( trimmed_len <= k - 1 ) {
                    offset = 0;
                    break;
                  }
                  offset = trimmed_len - k + 1;
                  trim_back( current_path );
                }
                for ( auto f = offset;
                    f < label_len && f + k < trav_path.get_sequence_len() + 1;
                    f += step ) {
                  bv_starts[f] = 1;
                }

                --bt_itr;
                trim_back( trav_path, *bt_itr );
                clear( current_path );
              }

              for ( std::size_t f = 0; f < label_len; ++f ) {
                if ( bv_starts[ f ] == 1 ) {
                  bv_starts[ f ] = 0;
                  this->add_start( id, f );
                }
              }

              clear( trav_path );
            }
          }

        inline void add_all_loci( unsigned int step=1 )
        {
          // TODO: Add documentation.
          // TODO: mention in the documentation that the `step` is approximately preserved in
          //       the whole graph.
          auto timer = stats_type( "add-starts" );

          seqan::Iterator<VarGraph, BFS>::Type itr( this->vargraph );

          unsigned long int prenode_remain = 0;
          unsigned long int remain_estimate = 0;
          VarGraph::nodeid_type prenode_level = 0;
          while ( !at_end( itr ) ) {
            if ( prenode_level != level( itr ) ) {
              prenode_remain = remain_estimate;
              remain_estimate = 0;
              prenode_level = level( itr );
            }

            auto seq_len = this->vargraph->node_length( *itr );
            unsigned long int cursor = ( step - prenode_remain ) % step;
            while ( cursor < seq_len ) {
              this->add_start( *itr, cursor );
              cursor += step;
            }

            unsigned long int new_remain;
            if ( step - prenode_remain > seq_len ) {
              new_remain = prenode_remain + seq_len;
            }
            else {
              new_remain = ( seq_len - step + prenode_remain ) % step;
            }

            if ( remain_estimate < new_remain ) {
              remain_estimate = new_remain;
            }

            ++itr;
          }
        }

          inline void
        traverse( std::function< void( typename TTraverser::output_type const& ) >& callback )
        {
          auto timer = stats_type( "traverse" );
          stats_type::set_total_nof_loci( this->starting_loci.size() );

          TTraverser traverser( this->vargraph, &(this->reads), &(this->reads_index), this->seed_len );
          for ( std::size_t idx = 0; idx < this->starting_loci.size(); ++idx )
          {
            const auto& locus = this->starting_loci[ idx ];
            stats_type::set_lastproc_locus( locus );

            traverser.set_start_locus( locus );
            traverser.run( callback );

            stats_type::set_lastdone_locus_idx( idx );
          }
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;
        std::vector< vg::Position > starting_loci;
        readsrecord_type reads;
        unsigned int seed_len;
        unsigned char seed_mismatches;  /**< @brief Allowed mismatches in a seed hit. */
        readsindex_type reads_index;
        /* ====================  METHODS       ======================================= */
          inline void
        index_reads( )
        {
          auto timer = stats_type( "index-reads" );
          this->reads_index = readsindex_type( this->reads.str );
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
