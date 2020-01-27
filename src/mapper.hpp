/**
 *    @file  mapper.hpp
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

#ifndef PSI_MAPPER_HPP__
#define PSI_MAPPER_HPP__

#include <fstream>
#include <type_traits>
#include <atomic>
#include <vector>
#include <unordered_set>
#include <iterator>
#include <functional>
#include <algorithm>

#include <sdsl/bit_vectors.hpp>
#include "stream.hpp"

#include "graph.hpp"
#include "sequence.hpp"
#include "index.hpp"
#include "index_iter.hpp"
#include "pathindex.hpp"
#include "utils.hpp"
#include "logger.hpp"
#include "stat.hpp"


namespace psi {
  /**
   *  @brief  MapperStat template class.
   *
   *  Collect statistics from a `Mapper` class instance(s) in running time.
   */
  template< typename TSpec = void >
    class MapperStat : public Timer<>
    {
      public:
        /* ====================  MEMBER TYPES  ======================================= */
        struct Coordinates {
          VarGraph::nodeid_type node_id;
          VarGraph::offset_type offset;
          Coordinates( VarGraph::nodeid_type nid, VarGraph::offset_type noff )
            : node_id( nid ), offset( noff ) { }
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
          get_lastproc_locus().store( Coordinates( value.node_id(), value.offset() ) );
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
        typedef Timer<> timer_type;
        struct Coordinates {
          VarGraph::nodeid_type node_id;
          VarGraph::offset_type offset;
        };
        /* ====================  LIFECYCLE     ======================================= */
        MapperStat( const std::string& ) { }
        ~MapperStat() { }
        /* ====================  METHODS       ======================================= */
        constexpr static inline timer_type::duration_type get_duration( const std::string& )
        {
          return timer_type::zero_duration;
        }
        constexpr static inline timer_type::rep_type get_duration_rep( const std::string& )
        {
          return timer_type::zero_duration_rep;
        }
        constexpr static inline const char* get_duration_str( const std::string& )
        {
          return "0";
        }
        constexpr static inline timer_type::duration_type get_lap( const std::string& )
        {
          return timer_type::zero_duration;
        }
        constexpr static inline timer_type::rep_type get_lap_rep( const std::string& )
        {
          return timer_type::zero_duration_rep;
        }
        constexpr static inline const char* get_lap_str( const std::string& )
        {
          return "0";
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
          seed_mismatches( mismatches ), traverser( vargraph, seed_len )
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
          pick_paths( PathIndex< TGraph, TText, TIndexSpec, TSequenceDirection >& paths,
              int n, bool patched=true,
              std::function< void( std::string const&, int ) > callback=nullptr )
          {
            if ( n == 0 ) return;
            auto timer = stats_type( "pick-paths" );

            paths.reserve( n * this->vargraph->path_count );
            seqan::Iterator< VarGraph, Haplotyper<> >::Type hap_itr( this->vargraph );
            auto context = paths.get_context();
            for ( std::size_t rank = 1; rank <= this->vargraph->max_path_rank(); ++rank ) {
              const auto& path_name = this->vargraph->path_name( rank );
              auto s = this->vargraph->node_at_path_position( path_name, 0 );
              go_begin( hap_itr, s );
              for ( int i = 0; i < n; ++i ) {
                if ( callback ) callback( path_name, i + 1 );
                if ( patched ) get_uniq_patched_haplotype( paths, hap_itr, context );
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
          seeds_on_paths( PathIndex< TGraph, TText, TIndexSpec, TSequenceDirection >& paths,
              std::function< void(typename TTraverser::output_type const &) > callback )
          {
            typedef TopDownFine< seqan::ParentLinks<> > TIterSpec;
            typedef typename PathIndex< TGraph, TText, TIndexSpec, TSequenceDirection >::index_type TPIndex;
            typedef typename seqan::Iterator< TPIndex, TIterSpec >::Type TPIterator;
            typedef typename seqan::Iterator< readsindex_type, TIterSpec >::Type TRIterator;

            if ( length( indexText( paths.index ) ) == 0 ) return;

            auto timer = stats_type( "paths-seed-find" );

            TPIterator piter( paths.index );
            TRIterator riter( this->reads_index );
            kmer_exact_matches( piter, riter, &paths, &(this->reads), this->seed_len, callback );
          }  /* -----  end of method template Mapper::seeds_on_paths  ----- */

        template< typename TPath, typename TSpec >
            inline void
          add_all_loci( PathSet< TPath, TSpec >& paths, unsigned int k,
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
                    !covered_by( current_path, paths ) ) {
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

        template< typename TPath, typename TSpec >
            inline unsigned long long int
          nof_uncovered_kmers( PathSet< TPath, TSpec >& paths, unsigned int k )
          {
            if ( this->starting_loci.size() == 0 ) return 0;
            auto timer = stats_type( "count-uncovered-kmer" );

            seqan::Iterator< VarGraph, Backtracker >::Type bt_itr( this->vargraph );
            Path< VarGraph > trav_path( this->vargraph );
            Path< VarGraph > current_path( this->vargraph );
            unsigned long long int uncovered = 0;

            long long int prev_id = 0;
            for ( const vg::Position& l : this->starting_loci ) {
              if ( prev_id == l.node_id() ) continue;
              prev_id = l.node_id();
              auto label_len = this->vargraph->node_length( l.node_id() );

              go_begin( bt_itr, l.node_id() );
              while ( !at_end( bt_itr ) ) {
                std::make_unsigned< VarGraph::offset_type >::type offset = label_len;
                extend_to_k( trav_path, bt_itr, offset - 1 + k );
                if ( trav_path.get_sequence_len() >= k ) current_path = trav_path;
                while ( current_path.get_sequence_len() != 0 &&
                    !covered_by( current_path, paths ) ) {
                  auto trimmed_len = current_path.get_sequence_len()
                    - this->vargraph->node_length( current_path.get_nodes().back() );
                  if ( trimmed_len <= k - 1 ) {
                    offset = 0;
                    break;
                  }
                  offset = trimmed_len - k + 1;
                  trim_back( current_path );
                }
                uncovered += label_len - offset;
                auto ub = trav_path.get_sequence_len() + 1 - k;
                if ( offset < ub && ub < label_len ) uncovered -= ub - offset;

                --bt_itr;
                trim_back( trav_path, *bt_itr );
                clear( current_path );
              }

              clear( trav_path );
            }

            return uncovered;
          }

          inline bool
        open_starts( const std::string& prefix, unsigned int seed_len,
            unsigned int step_size )
        {
          std::string filepath = prefix + "_loci_"
            "e" + std::to_string( step_size ) + "l" + std::to_string( seed_len );
          std::ifstream ifs( filepath, std::ifstream::in | std::ifstream::binary );
          if ( !ifs ) return false;

          std::function< void( vg::Position& ) > push_back =
            [this]( vg::Position& pos ) { this->starting_loci.push_back( pos ); };

          try {
            stream::for_each( ifs, push_back );
          }
          catch ( const std::runtime_error& ) {
            return false;
          }

          return true;
        }

          inline bool
        save_starts( const std::string& prefix, unsigned int seed_len,
            unsigned int step_size )
        {
          std::string filepath = prefix + "_loci_"
            "e" + std::to_string( step_size ) + "l" + std::to_string( seed_len );
          std::ofstream ofs( filepath, std::ofstream::out | std::ofstream::binary );
          if ( !ofs ) return false;

          std::function< vg::Position( uint64_t ) > lambda =
            [this]( uint64_t i ) { return this->starting_loci.at( i ); };

          try {
            stream::write( ofs, this->starting_loci.size(), lambda );
          }
          catch ( const std::runtime_error& ) {
            return false;
          }

          return true;
        }

          inline std::size_t
        get_nof_uniq_nodes( )
        {
          std::unordered_set< VarGraph::nodeid_type > set;
          for ( const auto& l : this->starting_loci ) set.insert( l.node_id() );
          return set.size();
        }

          inline void
        traverse( std::function< void( typename TTraverser::output_type const& ) > callback )
        {
          auto timer = stats_type( "traverse" );
          stats_type::set_total_nof_loci( this->starting_loci.size() );

          this->traverser.set_reads( &(this->reads) );
          this->traverser.set_reads_index( &(this->reads_index) );
          for ( std::size_t idx = 0; idx < this->starting_loci.size(); ++idx )
          {
            const auto& locus = this->starting_loci[ idx ];
            this->traverser.add_locus( locus );
            if ( this->starting_loci[ idx + 1 ].node_id() == locus.node_id() ) continue;

            stats_type::set_lastproc_locus( locus );
            this->traverser.run( callback );
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
        TTraverser traverser;
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
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_MAPPER_HPP__ --- */
