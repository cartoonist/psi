/**
 *    @file  traverser.h
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

#ifndef TRAVERSER_BASE_H__
#define TRAVERSER_BASE_H__

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

#include <seqan/seeds.h>

#include "vargraph.h"
#include "logger.h"
#include "sequence.h"
#include "index.h"
#include "index_iter.h"
#include "seed.h"
#include "stat.h"

// TODO: refactor: types (const, * and &).

namespace grem
{
  /* Forwards */
  template< typename TIndexSpec,
    typename TStrategy,
    template<typename> class TMatchingTraits,
    typename TStatSpec >
    class TraverserBase;

  template< typename TIter, std::size_t MaxMismatches >
    struct MatchingTraits {
      static const std::size_t max_mismatches = MaxMismatches;
      typedef struct {
        TIter iter;
        unsigned char mismatches;
        vg::Position pos;
      } TState;
    };

  template< typename TIter >
    using ExactMatching = MatchingTraits< TIter, 0 >;
  template< typename TIter >
    using ApproxMatching = MatchingTraits< TIter, 3 >;

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
    };  /* ----------  end of template class TraverserStat  ---------- */

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
    };  /* ----------  end of template class TraverserStat  ---------- */

  /**
   *  @brief  Stat template class specialization for `TraverserBase`.
   */
  template< typename TIndexSpec,
    typename TStrategy,
    template<typename> class TMatchingTraits,
    typename TSpec >
    class Stat< TraverserBase< TIndexSpec, TStrategy, TMatchingTraits, TSpec > >
    {
      public:
        typedef TraverserStat< TSpec > Type;
    };

  template< typename TIndexSpec,
    typename TStrategy,
    template<typename> class TMatchingTraits,
    typename TStatSpec >
    class TraverserBase
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        // :TODO:Tue Aug 29 14:49:\@cartoonist: Use grem::Seed class instead of Output.
        /**< @brief The output type. */
        typedef seqan::Seed < seqan::Simple > output_type;
        /**< @brief Index iterator template specialization parameter. */
        typedef TIndexSpec indexspec_type;
        typedef TopDownFine< seqan::ParentLinks<> > iterspec_type;
        typedef Dna5QStringSetIndex< TIndexSpec > index_type;
        typedef TIndexIter< index_type, iterspec_type > iterator_type;
        typedef TMatchingTraits< iterator_type > traits_type;
        typedef typename seqan::SAValue< index_type >::Type TSAValue;
        typedef typename Stat< TraverserBase >::Type stats_type;
        /* ====================  DATA MEMBERS  ======================================= */
        static const auto max_mismatches = traits_type::max_mismatches;
        /* ====================  LIFECYCLE      ====================================== */
        TraverserBase( const VarGraph* graph, index_type* index, unsigned int len, vg::Position s )
          : vargraph( graph ), reads_index( index ), seed_len( len ), start_locus( s )
        { }

        TraverserBase( const VarGraph* graph, index_type* index, unsigned int len )
          : vargraph( graph ), reads_index( index ), seed_len( len )
        {
          this->set_start_locus( 0, 0 );
        }
        /* ====================  ACCESSORS      ====================================== */
        /**
         *  @brief  getter function for vargraph.
         */
          inline const VarGraph*
        get_vargraph (  ) const
        {
          return this->vargraph;
        }  /* -----  end of method get_vargraph  ----- */

        /**
         *  @brief  getter function for reads_index.
         */
          inline const index_type*
        get_reads_index (  ) const
        {
          return this->reads_index;
        }  /* -----  end of method get_reads_index  ----- */

        /**
         *  @brief  getter function for seed_len.
         */
          inline unsigned int
        get_seed_len (  ) const
        {
          return this->seed_len;
        }  /* -----  end of method get_seed_len  ----- */

        /**
         *  @brief  getter function for start_locus.
         */
          inline vg::Position
        get_start_locus (  ) const
        {
          return this->start_locus;
        }  /* -----  end of method get_start_locus  ----- */
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
         *  @brief  setter function for reads_index.
         */
          inline void
        set_reads_index ( index_type* value )
        {
          this->reads_index = value;
        }  /* -----  end of method set_reads_index  ----- */

        /**
         *  @brief  setter function for seed_len.
         */
          inline void
        set_seed_len ( unsigned int value )
        {
          this->seed_len = value;
        }  /* -----  end of method set_seed_len  ----- */

        /**
         *  @brief  setter function for start_locus.
         */
          inline void
        set_start_locus ( vg::Position value )
        {
          this->start_locus = value;
        }  /* -----  end of method set_start_locus  ----- */

        /**
         *  @brief  setter function for start_locus.
         */
          inline void
        set_start_locus ( VarGraph::nodeid_type node_id, VarGraph::offset_type offset )
        {
          this->start_locus.set_node_id( node_id );
          this->start_locus.set_offset( offset );
        }  /* -----  end of method set_start_locus  ----- */
      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        const VarGraph* vargraph;      /**< @brief Pointer to variation graph. */
        index_type* reads_index;       /**< @brief Pointer to reads index. */
        unsigned int seed_len;         /**< @brief Seed length. */
        vg::Position start_locus;      /**< @brief Starting point. */
        std::vector< typename traits_type::TState > frontier_states;
    };  /* ----------  end of template class TraverserBase  ---------- */
}

#endif  // end of TRAVERSER_BASE_H__
