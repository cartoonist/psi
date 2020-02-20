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
#include "logger.hpp"
#include "sequence.hpp"
#include "index.hpp"
#include "index_iter.hpp"
#include "seed.hpp"
#include "stat.hpp"

// TODO: refactor: types (const, * and &).

namespace psi {
  /* Forwards */
  template< class TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TStatSpec >
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
   *  @brief  TraverserStat template class.
   *
   *  Collect statistics from a `TraverserBase` class instance(s) in running time.
   */
  template< typename TSpec = void >
    class TraverserStat
    {
      public:
        /* ====================  ACCESSORS     ======================================= */
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
        static inline std::atomic< double >& get_avg_pathlens( )
        {
          static std::atomic< double > avg_pathlens( 0 );
          std::lock_guard< std::mutex > lock( get_avg_pathlens_mutex() );
          update_avg_pathlens();
          return avg_pathlens;
        }
        /* ====================  METHODS       ======================================= */
        static inline void inc_total_nof_godowns( unsigned int by=1 )
        {
          get_total_nof_godowns().fetch_add( by );
        }
        static inline void reset_total_nof_godowns( )
        {
          get_total_nof_godowns().store( 0 );
        }
        static inline void inc_total_nof_paths( unsigned int by=1 )
        {
          get_total_nof_paths().fetch_add( by );
          inc_partial_nof_paths();
        }
        static inline void reset_total_nof_paths( )
        {
          get_total_nof_paths().store( 0 );
        }
        static inline void inc_pathlens_partial_sum( unsigned int by=1 )
        {
          std::lock_guard< std::mutex > lock( get_avg_pathlens_mutex() );
          if ( get_pathlens_partial_sum() > ULONG_MAX - by ) {
            update_avg_pathlens();
          }
          get_pathlens_partial_sum().fetch_add( by );
          inc_total_nof_paths();
        }
      private:
        /* ====================  ACCESSORS     ======================================= */
        static inline std::mutex& get_avg_pathlens_mutex( )
        {
          static std::mutex avg_pathlens_mutex;
          return avg_pathlens_mutex;
        }
        static inline std::atomic_ulong& get_pathlens_partial_sum( )
        {
          static std::atomic_ulong pathlens_partial_sum( 0 );
          return pathlens_partial_sum;
        }
        static inline std::atomic_ulong& get_partial_nof_paths( )
        {
          static std::atomic_ulong partial_nof_paths( 0 );
          return partial_nof_paths;
        }
        /* ====================  METHODS       ======================================= */
        static inline void update_avg_pathlens( )
        {
          assert( !get_avg_pathlens_mutex().try_lock() );  // should be already locked.

          if ( get_partial_nof_paths() == 0 ) return;

          std::atomic< double >& avg_pathlens = get_avg_pathlens();
          double pre_avg_pathlens = avg_pathlens;
          avg_pathlens = get_pathlens_partial_sum() /
            static_cast< double >( get_partial_nof_paths() );
          if ( pre_avg_pathlens != 0 ) {
            avg_pathlens = ( avg_pathlens + pre_avg_pathlens ) / 2.0;
          }
          reset_partial_nof_paths();
          reset_pathlens_partial_sum();
        }
        static inline void reset_pathlens_partial_sum( )
        {
          get_pathlens_partial_sum().store( 0 );
        }
        static inline void inc_partial_nof_paths( unsigned int by=1 )
        {
          get_partial_nof_paths().fetch_add( by );
        }
        static inline void reset_partial_nof_paths( )
        {
          get_partial_nof_paths().store( 0 );
        }
    };  /* --- end of template class TraverserStat --- */

  /**
   *  @brief  TraverserStat template class no-stat specialization.
   *
   *  Do nothing.
   */
  template< >
    class TraverserStat< NoStat >
    {
      public:
        /* ====================  ACCESSORS     ======================================= */
        static inline unsigned long long int get_total_nof_godowns( ) { return 0; }
        static inline unsigned long int get_total_nof_paths( ) { return 0; }
        static inline double get_avg_pathlens( ) { return 0; }
        /* ====================  METHODS       ======================================= */
        static inline void inc_total_nof_godowns( unsigned int by=1 ) { }
        static inline void reset_total_nof_godowns( ) { }
        static inline void inc_total_nof_paths( unsigned int by=1 ) { }
        static inline void reset_total_nof_paths( ) { }
        static inline void inc_pathlens_partial_sum( unsigned int by=1 ) { }
    };  /* --- end of template class TraverserStat --- */

  /**
   *  @brief  Stat template class specialization for `TraverserBase`.
   */
  template< class TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TSpec >
    class Stat< TraverserBase< TGraph, TIndex, TStrategy, TMatchingTraits, TSpec > >
    {
      public:
        typedef TraverserStat< TSpec > Type;
    };

  template< class TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TStatSpec >
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
        typedef typename Stat< TraverserBase >::Type stats_type;
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
