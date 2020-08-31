/**
 *    @file  traverser.hpp
 *   @brief  TraverserBase template class.
 *
 *  TraverserBase template class definition and its template specialisations.
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

#ifndef PSI_TRAVERSER_BASE_HPP__
#define PSI_TRAVERSER_BASE_HPP__

// :TODO:Fri Aug 25 13:15:\@cartoonist: remove unused headers.
// :TODO:Fri Aug 25 13:15:\@cartoonist: add required headers explicitly.
#include <stdexcept>
#include <cstdint>
#include <atomic>
#include <mutex>
#include <cmath>
#include <vector>
#include <array>
#include <functional>

#include "graph.hpp"
#include "graph_iter.hpp"
#include "sequence.hpp"
#include "index.hpp"
#include "index_iter.hpp"
#include "seed.hpp"
#include "stats.hpp"

// TODO: refactor: types (const, * and &).

namespace psi {
  /* Forwards */
  template< class TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TStatsSpec >
    class TraverserBase;

  template< typename TGraph, typename TIter, std::size_t MaxMismatches >
    struct MatchingTraits {
      static const std::size_t max_mismatches = MaxMismatches;
      typedef struct State {
        typedef typename seqan::Container< TIter >::Type TIndex;

        TIter iter;
        unsigned char mismatches;
        vg::Position spos;
        vg::Position cpos;
        size_t depth;
        bool end;

        State( TIndex* index, unsigned char mm,
            typename TGraph::id_type sid, typename TGraph::offset_type soffset,
            typename TGraph::id_type cid, typename TGraph::offset_type coffset, size_t d )
          : iter( *index ), mismatches( mm ), depth( d ), end( false )
        {
          spos.set_node_id( sid );
          spos.set_offset( soffset );
          cpos.set_node_id( cid );
          cpos.set_offset( coffset );
        }

        State( TIndex* index, unsigned char mm,
            typename TGraph::id_type sid, typename TGraph::offset_type soffset, size_t d )
          : State( index, mm, sid, soffset, sid, soffset, d )
        { }

        State( TIndex* index, unsigned char mm,
            vg::Position sp, vg::Position cp, size_t d )
          : iter( *index ), mismatches( mm ),
          spos( sp ), cpos( cp ), depth( d ), end( false )
        { }

        State( TIndex* index, unsigned char mm, vg::Position sp, size_t d )
          : State( index, mm, sp, sp, d )
        { }

        State( State const& ) = default;
        State( State&& ) = delete;
        State& operator=( State const& ) = default;
        State& operator=( State&& ) = delete;
        ~State( ) = default;
      } TState;
    };

  template< typename TGraph, typename TIter >
    using ExactMatching = MatchingTraits< TGraph, TIter, 0 >;
  template< typename TGraph, typename TIter >
    using ApproxMatching = MatchingTraits< TGraph, TIter, 3 >;

  /**
   *  @brief  TraverserStats template class.
   *
   *  Collect statistics from a `TraverserBase` class instance(s) in running time.
   */
  template< typename TSpec = WithStats >
    class TraverserStats
    {
      private:
        /* ====================  CONSTANTS     ======================================= */
        constexpr static const unsigned int RETRY_THRESHOLD = 4;
        constexpr static const unsigned long int PARTIAL_PATHLEN_SUM_UBOUND = ULONG_MAX - 65536;
        constexpr static const unsigned long int PLEN_NONE = -1;
      public:
        /* ====================  ACCESSORS     ======================================= */
        static inline std::atomic_ullong& get_total_seeds_off_paths( )
        {
          static std::atomic_ullong total_seeds_off_paths( 0 );
          return total_seeds_off_paths;
        }
        static inline std::atomic_ullong& get_total_nof_godowns( )
        {
          static std::atomic_ullong total_nof_godowns( 0 );
          return total_nof_godowns;
        }
        static inline std::atomic_ulong& get_total_nof_paths( )
        {
          static std::atomic_ulong total_nof_paths( 0 );
          return total_nof_paths;
        }
        /* ====================  METHODS       ======================================= */
        static inline void inc_total_seeds_off_paths( unsigned long long int by=1 )
        {
          get_total_seeds_off_paths().fetch_add( by );
        }
        static inline void reset_total_seeds_off_paths( )
        {
          get_total_seeds_off_paths().store( 0 );
        }
        static inline void inc_total_nof_godowns( unsigned long long int by=1 )
        {
          get_total_nof_godowns().fetch_add( by );
        }
        static inline void reset_total_nof_godowns( )
        {
          get_total_nof_godowns().store( 0 );
        }
        static inline void inc_total_nof_paths( unsigned long int by=1 )
        {
          get_total_nof_paths().fetch_add( by );
        }
        static inline void reset_total_nof_paths( )
        {
          get_total_nof_paths().store( 0 );
        }
        static inline void add_pathlen( unsigned long int len=1 )
        {
          unsigned int retry = RETRY_THRESHOLD;
          while ( true ) {
            auto peek_sum = get_partial_pathlen_sum().load();
            if ( peek_sum >= PARTIAL_PATHLEN_SUM_UBOUND ) {
              UniqWriterLock reducer( get_rws_lock() );
              if ( reducer ) {
                update_avg_pathlen();
              }
              continue;
            }
            {
              ReaderLock adder( get_rws_lock() );
              if ( get_partial_pathlen_sum().compare_exchange_weak(
                       peek_sum,
                       peek_sum+len,
                       std::memory_order_release,
                       std::memory_order_relaxed ) ) {
                inc_partial_nof_paths();
                break;
              }
            }

            if ( --retry == 0 ) {
              retry = RETRY_THRESHOLD;
              std::this_thread::yield();
            }
          }
        }
        static inline double compute_avg_pathlen( )
        {
          WriterLock lock( get_rws_lock() );
          update_avg_pathlen();
          return get_avg_pathlen().load();
        }
      private:
        /* ====================  ACCESSORS     ======================================= */
        static inline RWSpinLock<>& get_rws_lock( )
        {
          static RWSpinLock<> rws_lock;
          return rws_lock;
        }
        static inline std::atomic< double >& get_avg_pathlen( )
        {
          static std::atomic< double > avg_pathlen( PLEN_NONE );
          return avg_pathlen;
        }
        static inline std::atomic_ulong& get_partial_pathlen_sum( )
        {
          static std::atomic_ulong partial_pathlen_sum( 0 );
          return partial_pathlen_sum;
        }
        static inline std::atomic_ulong& get_partial_nof_paths( )
        {
          static std::atomic_ulong partial_nof_paths( 0 );
          return partial_nof_paths;
        }
        /* ====================  METHODS       ======================================= */
        static inline void update_avg_pathlen( )
        {
          assert( !get_rws_lock().acquire_writer_weak() );
          assert( get_partial_pathlen_sum().load() >= PARTIAL_PATHLEN_SUM_UBOUND );

          auto partial_sum = get_partial_pathlen_sum().load();
          auto partial_total = get_partial_nof_paths().load();
          reset_partial_nof_paths();
          reset_partial_pathlen_sum();

          double pre_avg_pathlen = get_avg_pathlen().load();
          double new_avg_pathlen = partial_sum / static_cast< double >( partial_total );
          if ( pre_avg_pathlen == PLEN_NONE ) get_avg_pathlen().store( new_avg_pathlen );
          else get_avg_pathlen().store( ( new_avg_pathlen + pre_avg_pathlen ) / 2.0 );
        }
        static inline void reset_partial_pathlen_sum( )
        {
          get_partial_pathlen_sum().store( 0 );
        }
        static inline void inc_partial_nof_paths( unsigned int by=1 )
        {
          get_partial_nof_paths().fetch_add( by );
        }
        static inline void reset_partial_nof_paths( )
        {
          get_partial_nof_paths().store( 0 );
        }
    };  /* --- end of template class TraverserStats --- */

  /**
   *  @brief  TraverserStats template class no-stats specialization.
   *
   *  Do nothing.
   */
  template< >
    class TraverserStats< NoStats >
    {
      public:
        /* ====================  ACCESSORS     ======================================= */
        constexpr static inline unsigned long long int get_total_seeds_off_paths( ) { return 0; }
        constexpr static inline unsigned long long int get_total_nof_godowns( ) { return 0; }
        constexpr static inline unsigned long int get_total_nof_paths( ) { return 0; }
        /* ====================  METHODS       ======================================= */
        constexpr static inline void inc_total_seeds_off_paths( unsigned long long int by=1 ) { }
        constexpr static inline void reset_total_seeds_off_paths( ) { }
        constexpr static inline void inc_total_nof_godowns( unsigned long long int by=1 ) { }
        constexpr static inline void reset_total_nof_godowns( ) { }
        constexpr static inline void inc_total_nof_paths( unsigned long int by=1 ) { }
        constexpr static inline void reset_total_nof_paths( ) { }
        constexpr static inline void add_pathlen( unsigned long int len=1 ) { }
        constexpr static inline double compute_avg_pathlen( ) { return 0; }
    };  /* --- end of template class TraverserStats --- */

  /**
   *  @brief  Stats template class specialization for `TraverserBase`.
   */
  template< class TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TSpec >
    class Stats< TraverserBase< TGraph, TIndex, TStrategy, TMatchingTraits, TSpec > >
    {
      public:
        typedef TraverserStats< TSpec > Type;
    };

  template< class TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TStatsSpec >
    class TraverserBase
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef Seed<> output_type;
        typedef TGraph graph_type;
        typedef TIndex index_type;
        typedef typename graph_type::id_type id_type;
        typedef typename graph_type::offset_type offset_type;
        typedef typename seqan::Spec< TIndex >::Type indexspec_type;
        typedef typename seqan::Fibre< TIndex, seqan::FibreText >::Type stringset_type;
        typedef typename seqan::Value< stringset_type >::Type text_type;
        typedef Records< stringset_type > records_type;
        typedef TopDownFine< > iterspec_type;
        typedef TIndexIter< TIndex, iterspec_type > iterator_type;
        typedef TMatchingTraits< graph_type, iterator_type > traits_type;
        typedef typename seqan::SAValue< TIndex >::Type TSAValue;
        typedef typename Stats< TraverserBase >::Type stats_type;
        /* ====================  DATA MEMBERS  ======================================= */
        static const auto max_mismatches = traits_type::max_mismatches;
        /* ====================  LIFECYCLE      ====================================== */
        TraverserBase( const graph_type* g, const records_type* r, TIndex* index,
            unsigned int len )
          : graph_ptr( g ), reads( r ), reads_index( index ), seed_len( len )
        { }

        TraverserBase( const graph_type* g, unsigned int len )
          : TraverserBase( g, nullptr, nullptr, len )
        { }
        /* ====================  ACCESSORS      ====================================== */
        /**
         *  @brief  getter function for graph_ptr.
         */
          inline const graph_type*
        get_graph_ptr (  ) const
        {
          return this->graph_ptr;
        }

        /**
         *  @brief  getter function for reads.
         */
          inline const records_type*
        get_reads( ) const
        {
          return this->reads;
        }

        /**
         *  @brief  getter function for reads_index.
         */
          inline const TIndex*
        get_reads_index (  ) const
        {
          return this->reads_index;
        }

        /**
         *  @brief  getter function for seed_len.
         */
          inline unsigned int
        get_seed_len (  ) const
        {
          return this->seed_len;
        }
        /* ====================  MUTATORS       ====================================== */
        /**
         *  @brief  setter function for graph_ptr.
         */
          inline void
        set_graph_ptr ( const graph_type* value )
        {
          this->graph_ptr = value;
        }

        /**
         *  @brief  setter function for reads.
         */
          inline void
        set_reads( const records_type* value )
        {
          this->reads = value;
        }

        /**
         *  @brief  setter function for reads_index.
         */
          inline void
        set_reads_index ( TIndex* value )
        {
          this->reads_index = value;
        }

        /**
         *  @brief  setter function for seed_len.
         */
          inline void
        set_seed_len ( unsigned int value )
        {
          this->seed_len = value;
        }

          inline void
        add_locus( vg::Position p )
        {
          this->states.emplace_back(
              this->reads_index,
              max_mismatches + 1,
              std::move( p ),
              0 );
        }

          inline void
        add_locus( id_type id, offset_type offset )
        {
          this->states.emplace_back(
              this->reads_index,
              max_mismatches + 1,
              id,
              offset,
              0 );
        }

          inline void
        states_reserve( size_t size )
        {
          this->states.reserve( size );
        }
      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        const graph_type* graph_ptr;   /**< @brief Pointer to variation graph. */
        const records_type* reads;     /**< @brief Pointer to reads record. */
        TIndex* reads_index;           /**< @brief Pointer to reads index. */
        unsigned int seed_len;         /**< @brief Seed length. */
        std::vector< typename traits_type::TState > states;
    };  /* --- end of template class TraverserBase --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_TRAVERSER_BASE_HPP__ --- */
