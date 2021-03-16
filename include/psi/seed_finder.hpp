/**
 *    @file  seed_finder.hpp
 *   @brief  SeedFinder template class.
 *
 *  SeedFinder template class definition and its template specialisations.
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

#ifndef PSI_SEED_FINDER_HPP__
#define PSI_SEED_FINDER_HPP__

#include <fstream>
#include <type_traits>
#include <atomic>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <functional>
#include <algorithm>

#include <sdsl/bit_vectors.hpp>
#include <pairg/reachability.hpp>
#include <KokkosKernels_IOUtils.hpp>

#include "graph.hpp"
#include "traverser.hpp"
#include "sequence.hpp"
#include "index.hpp"
#include "index_iter.hpp"
#include "pathindex.hpp"
#include "crs_matrix.hpp"
#include "utils.hpp"
#include "stats.hpp"


namespace psi {
  struct SeedFinderStatsBase {
    enum progress_type : int {
      finder_off,
      instantiated,
      load_pindex,
      load_starts,
      load_dindex,
      select_paths,
      create_pindex,
      find_uncovered,
      create_dindex,
      write_pindex,
      write_starts,
      write_dindex,
      ready,
    };

    constexpr static const char* progress_table[] = {
      "Progress is not tracked",
      "Just instantiated",
      "Loading path index",
      "Loading starting loci",
      "Loading distance index",
      "Selecting paths from input graph",
      "Indexing paths sequences",
      "Finding uncovered loci",
      "Creating distance index",
      "Writing path index",
      "Writing starting loci",
      "Writing distance index",
      "Ready for seed finding",
    };

    enum thread_progress_type : int {
      thread_off,
      sleeping,
      seed_chunk,
      index_chunk,
      find_on_paths,
      find_off_paths,
      find_mems,
      query_dindex,
    };

    constexpr static const char* thread_progress_table[] = {
      "Progress is not tracked",
      "Zzz",
      "Seeding a read chunk",
      "Indexing a read chunk",
      "Finding seeds on paths",
      "Finding seeds off paths",
      "Finding MEMs on paths",
      "Querying distance index",
    };
  };  /* --- end of template class SeedFinderStats --- */

  /**
   *  @brief  SeedFinderStats template class.
   *
   *  Collect running-time statistics from a `SeedFinder` class instance.
   */
  template< class TSeedFinder, typename TSpec = WithStats >
    class SeedFinderStats
      : public SeedFinderStatsBase
    {
      public:
        /* === TYPE MEMBERS === */
        typedef TSeedFinder seedfinder_type;
        typedef SeedFinderStatsBase base_type;
        using base_type::progress_type;
        using base_type::thread_progress_type;

        struct ThreadStats {
          private:
            /* === CONSTANTS === */
            constexpr static const double GOCC_AVG_NONE = -1;
            constexpr static const unsigned long long int GOCC_UBOUND = ULLONG_MAX / 2;

            /* === DATA MEMBERS === */
            thread_progress_type progress;
            unsigned int chunks_done;
            std::size_t locus_idx;
            unsigned long long int gocc_sum;  /**< @brief Sum of genome occurrence counts */
            unsigned long long int gocc_tot;  /**< @brief No. of seeds contributed in the sum */
            double gocc_avg;  /**< @brief Average seed genome occurrence count */
            unsigned long long int gocc_skips;  /**< @brief No. of skipped seeds because of high gocc */

          public:
            /* === LIFECYCLE === */
            ThreadStats( )
              : progress( thread_progress_type::sleeping ), chunks_done( 0 ),
                locus_idx( 0 ), gocc_sum( 0 ), gocc_tot( 0 ), gocc_avg( GOCC_AVG_NONE ),
                gocc_skips( 0 )
            { }

            /* === ACCESSORS === */
              inline thread_progress_type
            get_progress( ) const
            {
              return this->progress;
            }

              inline const char*
            get_progress_str( ) const
            {
              return base_type::thread_progress_table[ this->progress ];
            }

              inline unsigned int
            get_chunks_done( ) const
            {
              return this->chunks_done;
            }

              inline std::size_t
            get_locus_idx( ) const
            {
              return this->locus_idx;
            }

              inline unsigned long long int
            get_gocc_skips( ) const
            {
              return this->gocc_skips;
            }

            /* === MUTATORS === */
              inline void
            set_progress( thread_progress_type value )
            {
              this->progress = value;
            }

              inline void
            set_chunks_done( unsigned int value )
            {
              this->chunks_done = value;
            }

              inline void
            inc_chunks_done( )
            {
              ++this->chunks_done;
            }

              inline void
            set_locus_idx( std::size_t value )
            {
              this->locus_idx = value;
            }

              inline void
            set_gocc_skips( unsigned long long int value )
            {
              this->gocc_skips = value;
            }

              inline void
            inc_gocc_skips( )
            {
              ++this->gocc_skips;
            }

            /* === METHODS === */
              inline void
            add_seed_gocc( unsigned long long int count )
            {
              if ( this->gocc_sum >= GOCC_UBOUND ) update_avg_seed_gocc();
              this->gocc_sum += count;
              this->gocc_tot++;
            }

              inline void
            update_avg_seed_gocc()
            {
              this->gocc_avg = this->avg_seed_gocc();
              this->gocc_sum = 0;
              this->gocc_tot = 0;
            }

              inline double
            avg_seed_gocc( ) const
            {
              if ( this->gocc_tot != 0 ) {
                double new_avg = this->gocc_sum / static_cast< double >( this->gocc_tot );
                if ( this->gocc_avg == GOCC_AVG_NONE ) return new_avg;
                else return ( this->gocc_avg + new_avg ) / 2.0;
              }
              return ( this->gocc_avg != GOCC_AVG_NONE ? this->gocc_avg : 0 );
            }
        };
        typedef std::unordered_map< std::string, ThreadStats > container_type;
        typedef typename container_type::size_type size_type;
        typedef Timer<> timer_type;
        /* === CONST MEMBERS === */
        constexpr static const int CLS_ID_LEN = 8;
      private:
        /* === DATA MEMBERS === */
        seedfinder_type const* finder_ptr;
        progress_type progress;
        container_type tstats;
        std::string id;

        /* === STATIC MEMBERS === */
          static inline SeedFinderStats*&
        get_instance_ptr( )
        {
          static SeedFinderStats* instance_ptr = nullptr;
          return instance_ptr;
        }

          static inline void
        set_instance_ptr( SeedFinderStats* value )
        {
          SeedFinderStats::get_instance_ptr() = value;
        }

      public:
        /* === STATIC MEMBERS === */
          static inline SeedFinderStats const*
        get_instance_const_ptr( )
        {
          return SeedFinderStats::get_instance_ptr();
        }

          static inline void
        signal_handler( int signo )
        {
          std::cout << "\n====  "
                    << "Received \"" << strsignal( signo ) << "\" (" << signo << ")"
                    << "  ====" << std::endl;
          if ( SeedFinderStats::get_instance_ptr() == nullptr ) {
            std::cout << "No tracking seed finder!" << std::endl;
            return;
          }
          std::cout << "PSI seed finder last status: "
                    << SeedFinderStats::get_instance_ptr()->get_progress_str() << std::endl;
          auto const& threads_stats = SeedFinderStats::get_instance_ptr()->get_threads_stats();
          std::cout << ( threads_stats.empty() ? "No" : std::to_string( threads_stats.size() ) )
                    << " running thread(s)" << ( threads_stats.empty() ? "." : ":" )
                    << std::endl;
          int tid = 0;
          for ( auto const& stats : threads_stats ) {
            std::cout << stats.first << " -- Thread: " << ++tid << std::endl;
            std::cout << stats.first << " -- Last status: " << stats.second.get_progress_str()
                      << std::endl;
            std::cout << stats.first << " -- Chunks done: "
                      << stats.second.get_chunks_done() << std::endl;
            std::cout << stats.first << " -- Average seed genome occurrence count: "
                      << stats.second.avg_seed_gocc() << std::endl;
            std::cout << stats.first << " -- Skipped seeds because of high genome occurrence count: "
                      << stats.second.get_gocc_skips() << std::endl;
            if ( stats.second.get_progress() == thread_progress_type::find_off_paths ) {
              auto loc_idx = stats.second.get_locus_idx();
              auto loc_num = SeedFinderStats::get_instance_ptr()->get_ptr()->get_starting_loci().size();
              auto wlen = std::to_string( loc_num ).length();
              std::cout << stats.first << " -- Traversed loci: " << std::setw(wlen)
                        << loc_idx << " / " << std::setw(wlen) << loc_num
                        << " [%" << std::setw(3) << loc_idx * 100 / loc_num << "]"
                        << std::endl;
              auto pos = SeedFinderStats::get_instance_ptr()->get_ptr()->get_starting_loci()[ loc_idx ];
              std::cout << stats.first << " -- Last traversed locus: "
                        << "(" << pos.node_id() << ", " << pos.offset() << ")" << std::endl;
            }
            SeedFinderStats::get_instance_ptr()->for_each_timer(
                stats.first,
                [&stats]( auto name, auto period ) {
                  std::cout << stats.first << " -- Timer '" << name << "': "
                            << period.get_lap().str() << std::endl;
                  return true;
                } );
            std::cout << std::endl;
          }
          bool first = true;
          SeedFinderStats::get_instance_ptr()->for_each_timer(
              [&first]( auto name, auto period ) {
                if ( first ) {
                  std::cout << "All timers" << std::endl;
                  std::cout << "----------" << std::endl;
                  first = false;
                }
                std::cout << "Timer '" << name << "': "
                          << period.get_lap().str() << std::endl;
                return true;
              } );
          if ( !first ) {
            std::cout << "----------" << std::endl;
          }
        }

        /* === LIFECYCLE === */
        SeedFinderStats( seedfinder_type const* ptr )
          : finder_ptr( ptr ), progress( progress_type::instantiated ),
            id( random::random_string( SeedFinderStats::CLS_ID_LEN ) )
        {
          SeedFinderStats::set_instance_ptr( this );  // keep the last instance as tracked one
        }

        ~SeedFinderStats( ) noexcept
        {
          if ( SeedFinderStats::get_instance_ptr() == this ) {
            SeedFinderStats::set_instance_ptr( nullptr );
          }
        }

        /* === ACCESSORS === */
        /**
         *  @brief  Get the pointer to the SeedFinder instance.
         */
          inline seedfinder_type const*
        get_ptr( ) const
        {
          return this->finder_ptr;
        }

        /**
         *  @brief  Get the progress code.
         */
          inline progress_type
        get_progress( ) const
        {
          return this->progress;
        }

        /**
         *  @brief  Get the progress string.
         */
          inline const char*
        get_progress_str( ) const
        {
          return base_type::progress_table[ this->progress ];
        }

        /**
         *  @brief  Get the statistics for the thread with the given `id`.
         */
          inline ThreadStats&
        get_thread_stats( std::string const& id )
        {
          return this->tstats[ id ];
        }

        /**
         *  @brief  Get the statistics for this thread.
         */
          inline ThreadStats&
        get_this_thread_stats( )
        {
          return this->get_thread_stats( this->get_this_thread_id() );
        }

        /**
         *  @brief  Get the statistics for all threads.
         */
          inline container_type const&
        get_threads_stats( ) const
        {
          return this->tstats;
        }

        /* === MUTATORS === */
        /**
         *  @brief  Set the progress code.
         */
          inline void
        set_progress( progress_type value )
        {
          this->progress = value;
        }

          inline void
        set_as_tracked( )
        {
          SeedFinderStats::set_instance_ptr( this );
        }

        /* === METHODS === */
          inline timer_type
        timeit( std::string const& name, std::string const& thread_id ) const
        {
          return timer_type( this->id + name + thread_id );
        }

        /* timeit thread-safe */
          inline timer_type
        timeit_ts( std::string const& name ) const
        {
          return this->timeit( name, this->get_this_thread_id() );
        }

          inline timer_type
        timeit( std::string const& name ) const
        {
          return timer_type( this->id + name );
        }

          inline timer_type::period_type
        get_timer( std::string const& name ) const
        {
          return timer_type::get_timers()[ this->id + name ];
        }

          inline timer_type::period_type
        get_timer( std::string const& name, std::string const& thread_id ) const
        {
          return timer_type::get_timers()[ this->id + name + thread_id ];
        }

        template< typename TCallback >
          inline bool
        for_each_timer( TCallback callback ) const
        {
          for ( auto const& timer : timer_type::get_timers() ) {
            if ( starts_with( timer.first, this->id ) &&
                 !callback( timer.first.substr( SeedFinderStats::CLS_ID_LEN ),
                            timer.second ) ) return false;
          }
          return true;
        }

        template< typename TCallback >
          inline bool
        for_each_timer( std::string const& thread_id, TCallback callback ) const
        {
          auto thread_id_len = thread_id.size();
          for ( auto const& timer : timer_type::get_timers() ) {
            if ( starts_with( timer.first, this->id ) &&
                 ends_with( timer.first, thread_id ) &&
                 !callback( timer.first.substr( SeedFinderStats::CLS_ID_LEN,
                                                timer.first.size()
                                                - thread_id_len
                                                - SeedFinderStats::CLS_ID_LEN ),
                            timer.second ) ) return false;
          }
          return true;
        }
      private:
        /* === METHODS === */
          inline std::string const&
        get_this_thread_id() const
        {
          thread_local static const std::string thread_id = get_thread_id();
          return thread_id;
        }
    };  /* --- end of template class SeedFinderStats --- */

  /**
   *  @brief  SeedFinderStats template specialization for no-stats.
   *
   *  Do nothing.
   */
  template< class TSeedFinder >
    class SeedFinderStats< TSeedFinder, NoStats >
      : public SeedFinderStatsBase
    {
      public:
        /* === TYPE MEMBERS === */
        typedef TSeedFinder seedfinder_type;
        typedef SeedFinderStatsBase base_type;
        using base_type::progress_type;
        using base_type::thread_progress_type;

        struct ThreadStats {
          public:
            /* === ACCESSORS === */
              constexpr inline thread_progress_type
            get_progress( ) const
            {
              return thread_progress_type::thread_off;
            }

              constexpr inline const char*
            get_progress_str( ) const
            {
              return base_type::thread_progress_table[ thread_progress_type::thread_off ];
            }

              constexpr inline unsigned int
            get_chunks_done( ) const
            {
              return 0;
            }

              constexpr inline std::size_t
            get_locus_idx( ) const
            {
              return 0;
            }

              constexpr inline unsigned long long int
            get_gocc_skips( ) const
            {
              return 0;
            }

            /* === MUTATORS === */
              constexpr inline void
            set_progress( thread_progress_type )
            { /* noop */ }

              constexpr inline void
            set_chunks_done( unsigned int )
            { /* noop */ }

              constexpr inline void
            inc_chunks_done( )
            { /* noop */ }

              constexpr inline void
            set_locus_idx( std::size_t )
            { /* noop */ }

              constexpr inline void
            set_gocc_skips( unsigned long long int )
            { /* noop */ }

              constexpr inline void
            inc_gocc_skips( )
            { /* noop */ }

            /* === METHODS === */
              constexpr inline void
            add_seed_gocc( unsigned long long int )
            { /* noop */}

              constexpr inline void
            update_avg_seed_gocc()
            { /* noop  */ }

              constexpr inline double
            avg_seed_gocc( )
            { return 0; }
        };
        typedef std::unordered_map< std::string, ThreadStats > container_type;
        typedef typename container_type::size_type size_type;
        typedef Timer< void > timer_type;
        /* === STATIC MEMBERS === */
          static inline void
        signal_handler( int signo )
        {
          std::cout << "\n====  "
                    << "Received \"" << strsignal( signo ) << "\" (" << signo << ")"
                    << "  ====" << std::endl;
          std::cout << "PSI seed finder last status: "
                    << base_type::progress_table[ progress_type::finder_off ] << std::endl;
        }

      private:
        /* === DATA MEMBERS === */
        ThreadStats null_stats;

        /* === STATIC MEMBERS === */
          static inline SeedFinderStats*&
        get_instance_ptr( )
        {
          static SeedFinderStats* instance_ptr = nullptr;
          return instance_ptr;
        }

          static inline void
        set_instance_ptr( SeedFinderStats* value )
        {
          SeedFinderStats::get_instance_ptr() = value;
        }

      public:
        /* === STATIC MEMBERS === */
        static inline SeedFinderStats const*
        get_instance_const_ptr( )
        {
          return SeedFinderStats::get_instance_ptr();
        }

        /* === LIFECYCLE === */
        SeedFinderStats( seedfinder_type const* )
        {
          SeedFinderStats::set_instance_ptr( this );
        }

        ~SeedFinderStats( ) noexcept
        {
          if ( SeedFinderStats::get_instance_ptr() == this ) {
            SeedFinderStats::set_instance_ptr( nullptr );
          }
        }

        /* === ACCESSORS === */
          constexpr inline seedfinder_type const*
        get_ptr( ) const
        {
          return nullptr;
        }
          constexpr inline progress_type
        get_progress( ) const
        {
          return progress_type::finder_off;
        }

          constexpr inline const char*
        get_progress_str( ) const
        {
          return base_type::progress_table[ progress_type::finder_off ];
        }

          inline ThreadStats&
        get_thread_stats( std::string const& )
        {
          return this->null_stats;
        }

          inline ThreadStats&
        get_this_thread_stats( )
        {
          return this->null_stats;
        }

          inline container_type
        get_threads_stats( ) const
        {
          return container_type();
        }

        /* === MUTATORS === */
          constexpr inline void
        set_progress( progress_type ) { /* noop */ }

          constexpr inline void
        set_as_tracked( )
        {
          SeedFinderStats::set_instance_ptr( this );
        }

        /* === METHODS === */
          constexpr inline timer_type
        timeit( std::string const&, std::string const& ) const
        {
          return timer_type( );
        }

          constexpr inline timer_type
        timeit_ts( std::string const& ) const
        {
          return timer_type( );
        }

          constexpr inline timer_type
        timeit( std::string const& ) const
        {
          return timer_type( );
        }

          constexpr inline timer_type::period_type
        get_timer( std::string const& ) const
        {
          return timer_type::period_type();
        }

          constexpr inline timer_type::period_type
        get_timer( std::string const&, std::string const& ) const
        {
          return timer_type::period_type();
        }

        template< typename TCallback >
          constexpr inline bool
        for_each_timer( TCallback ) const
        {
          return true;
        }

        template< typename TCallback >
          constexpr inline bool
        for_each_timer( std::string const&, TCallback ) const
        {
          return true;
        }
    };  /* --- end of template class SeedFinderStats --- */

  template< typename TGraphSpec = gum::Succinct,
            typename TReadsStringSet = Dna5QStringSet<>,
            typename TReadsIndexSpec = seqan::IndexWotd<>,
            typename TPathsStringSetSpec = DiskBased,
            typename TStrategy = BFS,
            template<typename, typename> class TMatchingTraits = ExactMatching >
  struct SeedFinderTraits {
    typedef gum::SeqGraph< TGraphSpec > graph_type;
    typedef seqan::Index< TReadsStringSet, TReadsIndexSpec > seedindex_type;

    template< typename TStatsSpec = NoStats >
      using traverser_type = typename Traverser< graph_type, seedindex_type,
                                                 TStrategy, TMatchingTraits,
                                                 TStatsSpec >::Type;

    typedef TPathsStringSetSpec pathstrsetspec_type;
  };

  template< typename TStatsSpec = NoStats,
            typename TTraits = SeedFinderTraits<> >
    class SeedFinder;

  /**
   *  @brief  `Stats` template class specialization for `SeedFinder`.
   */
  template< typename TStatsSpec, typename TTraits >
    class Stats< SeedFinder< TStatsSpec, TTraits > >
    {
      typedef SeedFinder< TStatsSpec, TTraits > seedfinder_type;
      public:
        typedef SeedFinderStats< seedfinder_type, TStatsSpec > type;
    };  /* --- end of template class Stats --- */

  template< typename TStatsSpec, typename TTraits >
    class SeedFinder
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TTraits traits_type;
        typedef typename traits_type::graph_type graph_type;
        typedef typename traits_type::template traverser_type< TStatsSpec > traverser_type;
        typedef typename traits_type::pathstrsetspec_type pathstrsetspec_type;
        typedef typename graph_type::id_type id_type;
        typedef typename graph_type::offset_type offset_type;
        typedef typename graph_type::rank_type rank_type;
        typedef StatsType< SeedFinder > stats_type;
        typedef typename stats_type::progress_type progress_type;
        typedef typename stats_type::thread_progress_type thread_progress_type;
        typedef Records< typename traverser_type::stringset_type > readsrecord_type;
        typedef typename traverser_type::index_type readsindex_type;
        typedef YaString< pathstrsetspec_type > text_type;
        typedef PathIndex< graph_type, text_type, psi::FMIndex<>, Reversed > pathindex_type;
        typedef pairg::matrixOps crs_traits_type;
        typedef CRSMatrix< crs_matrix::Compressed, bool, uint32_t, uint64_t > crsmat_type;
        typedef crs_matrix::Buffered mutable_crsmat_spec_type;
        typedef make_spec_t< mutable_crsmat_spec_type, crsmat_type > mutable_crsmat_type;

        class KokkosHandler {
        public:
          /* === LIFECYCLE === */
          KokkosHandler( bool fin=true ) : finaliser( fin )
          {
            /* Ensure no concurrent dtor is running */
            ReaderLock constructor( KokkosHandler::get_final_lock() );
            if ( !Kokkos::is_initialized() ) {
              /* Ensure only one of the concurrent ctor gets the lock */
              UniqWriterLock initialiser( KokkosHandler::get_init_lock() );
              if ( initialiser ) {
                /* Initialise Kokkos */
                Kokkos::initialize();
              }
            }
            /* Sync ctors with the initialiser ctor to ensure that the object is ready */
            WriterLock sync( KokkosHandler::get_init_lock() );
          }

          ~KokkosHandler() noexcept
          {
            /* Ensure no concurrent ctor or dtor is running */
            WriterLock lock( KokkosHandler::get_final_lock() );
            /* Finalise Kokkos if we are responsible/finaliser */
            /* NOTE: No initialisation can be done after finalising. */
            if ( Kokkos::is_initialized() && this->finaliser ) KokkosHandler::finalise();
          }
          /* === OPERATORS === */
          inline
          operator bool() const
          {
            return this->finaliser;
          }

          inline KokkosHandler&
          operator=( bool value )
          {
            this->finaliser = value;
            return *this;
          }
          /* === STATIC MEMBERS === */
          static inline void
          finalise( )
          {
            Kokkos::finalize();
          }

          static inline RWSpinLock<>&
          get_init_lock( )
          {
            static RWSpinLock<> init_lock;
            return init_lock;
          }

          static inline RWSpinLock<>&
          get_final_lock( )
          {
            static RWSpinLock<> final_lock;
            return final_lock;
          }
          /* === DATA MEMBERS === */
          bool finaliser;
        };
        /* ====================  LIFECYCLE      ====================================== */
        SeedFinder( const graph_type& g,
            unsigned int len,
            unsigned int gocc_thr = 0,
            unsigned char mismatches = 0 )
          : graph_ptr( &g ), pindex( g, true ), handler( true ), seed_len( len ),
          seed_mismatches( mismatches ),
          gocc_threshold( ( gocc_thr != 0 ? gocc_thr : UINT_MAX ) ),
          stats_ptr( std::make_unique< stats_type >( this ) )
        { }
        /* ====================  ACCESSORS      ====================================== */
        /**
         *  @brief  getter function for graph_ptr.
         */
          inline const graph_type*
        get_graph_ptr( ) const
        {
          return this->graph_ptr;
        }

        /**
         *  @brief  getter function for starting_loci.
         */
          inline const std::vector< vg::Position >&
        get_starting_loci( ) const
        {
          return this->starting_loci;
        }

        /**
         *  @brief  getter function for seed_len.
         */
          inline unsigned int
        get_seed_len( ) const
        {
          return this->seed_len;
        }

        /**
         *  @brief  getter function for seed_mismatches.
         */
          inline unsigned char
        get_seed_mismatches( ) const
        {
          return this->seed_mismatches;
        }

        /**
         *  @brief  getter function for gocc_threshold.
         */
          inline unsigned int
        get_gocc_threshold( ) const
        {
          return this->gocc_threshold;
        }

        /**
         *  @brief  getter function for `handler.finaliser`.
         */
          inline bool
        is_finaliser( ) const
        {
          return this->handler;
        }

        /**
         *  @brief  getter function for pindex.
         */
          inline pathindex_type const&
        get_pindex( ) const
        {
          return this->pindex;
        }

        /**
         *  @brief  getter function for distance index matrix.
         */
          inline crsmat_type const&
        get_distance_matrix( ) const
        {
          return this->distance_mat;
        }

        /**
         * @brief  getter function for stats_ptr.
         */
          inline stats_type const&
        get_stats() const
        {
          return *this->stats_ptr;
        }

        /* ====================  MUTATORS       ====================================== */
        /**
         *  @brief  setter function for graph_ptr.
         */
          inline void
        set_graph_ptr( const graph_type* value )
        {
          this->graph_ptr = value;
        }

        /**
         *  @brief  setter function for starting_loci.
         *
         *  Copy assignment.
         */
          inline void
        set_starting_loci( const std::vector< vg::Position >& loci )
        {
          this->starting_loci = loci;
        }

        /**
         *  @brief  setter function for starting_loci.
         *
         *  Move assignment.
         */
          inline void
        set_starting_loci( std::vector< vg::Position >&& loci )
        {
          this->starting_loci = std::move( loci );
        }

        /**
         *  @brief  setter function for seed_len.
         */
          inline void
        set_seed_len( unsigned int value )
        {
          this->seed_len = value;
        }

        /**
         *  @brief  setter function for seed_mismatches.
         */
          inline void
        set_seed_mismatches( unsigned char value )
        {
          this->seed_mismatches = value;
        }

        /**
         *  @brief setter function for gocc_threshold
         */
          inline void
        set_gocc_threshold( unsigned int value )
        {
          this->gocc_threshold = value;
        }

          inline void
        set_as_finaliser( )
        {
          this->handler = true;
        }

          inline void
        unset_as_finaliser( )
        {
          this->handler = false;
        }

          inline readsrecord_type
        create_readrecord( ) const
        {
          return readsrecord_type();
        }

          inline readsindex_type
        index_reads( readsrecord_type const& reads ) const
        {
          this->stats_ptr->get_this_thread_stats().set_progress(
              thread_progress_type::index_chunk );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "index-reads" );

          return readsindex_type( reads.str );
        }

        template< typename T >
          inline void
        get_seeds( readsrecord_type& seeds, readsrecord_type const& reads,
                   T distance ) const
        {
          this->stats_ptr->get_this_thread_stats().set_progress(
              thread_progress_type::seed_chunk );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "seeding" );

          seeding( seeds, reads, this->seed_len, distance );
        }

          inline void
        add_start( const vg::Position& locus )
        {
          this->starting_loci.push_back( locus );
        }

          inline void
        add_start( id_type node_id, offset_type offset )
        {
          vg::Position locus;
          locus.set_node_id( node_id );
          locus.set_offset( offset );
          this->add_start( locus );
        }
        /* ====================  METHODS        ====================================== */
        /**
         *  @brief  Pick n paths from the graph.
         *
         *  @param[out]  paths Set of paths to be generated.
         *  @param[in]  n Number of paths.
         *
         *  This method generates a set of (probably) unique whole-genome paths from the
         *  graph.
         *
         *  XXX: We assume that each connect component in the graph has one and only one
         *  path indicating a sample haplotype in that region.
         */
            inline void
          pick_paths( unsigned int n, bool patched=true, unsigned int context=0,
              std::function< void( std::string const&, int ) > callback=nullptr,
              std::function< void( std::string const& ) > info=nullptr,
              std::function< void( std::string const& ) > warn=nullptr )
          {
            if ( n == 0 ) return;

            this->stats_ptr->set_progress( progress_type::select_paths );
            [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "pick-paths" );

            this->pindex.reserve( n * this->graph_ptr->get_path_count() );
            auto hp_itr = begin( *this->graph_ptr, Haplotyper<>() );
            auto hp_end = end( *this->graph_ptr, Haplotyper<>() );
            context = this->set_context( context, patched, info, warn );
            this->graph_ptr->for_each_path(
                [&]( auto path_rank, auto path_id ) {
                  auto path_name = this->graph_ptr->path_name( path_id );
                  id_type s = *this->graph_ptr->path( path_id ).begin();
                  hp_itr.reset( s );
                  for ( unsigned int i = 0; i < n; ++i ) {
                    if ( callback ) callback( path_name, i + 1 );
                    get_uniq_haplotype( this->pindex, hp_itr, hp_end, context, patched );
                  }
                  return true;
                } );
          }

        inline void
        index_paths( )
        {
          this->stats_ptr->set_progress( progress_type::create_pindex );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "index-paths" );

          this->pindex.create_index();
        }

      /**
       *  @brief  Create distance index matrix.
       *
       *  NOTE: This method assumes that the input graph is sorted such that node rank
       *  ranges in components are disjoint.
       *
       *  NOTE: This function assumes that the graph is augmented by one path per region
       *        and nothing more.
       */
        inline void
        create_distance_index( unsigned int dmin=0, unsigned int dmax=0,
                               std::function< void( std::string const& ) > info=nullptr,
                               std::function< void( std::string const& ) > warn=nullptr )
        {
          if ( dmin == 0 || dmax < dmin ) return;  // not constructible
          if ( dmax == 0 ) dmax = dmin;

          this->stats_ptr->set_progress( progress_type::create_dindex );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "index-distances" );

          auto provider = [this, dmin, dmax, info]( auto callback ) {
            /* Extract component boundary nodes */
            auto comp_ranks = util::components_ranks( *this->graph_ptr );
            comp_ranks.push_back( 0 );  // add the upper bound of the last component

            if ( info ) {
              info( "Constructing distance index for " +
                    std::to_string( comp_ranks.size()-1 ) + " regions..." );
            }

            for ( std::size_t idx = 0; idx < comp_ranks.size()-1; ++idx ) {
              auto adj_mat = util::adjacency_matrix( *this->graph_ptr, crs_traits_type(),
                                                     comp_ranks[idx], comp_ranks[idx+1] );
              auto dist_mat = pairg::buildValidPairsMatrix( adj_mat, dmin, dmax );
              auto sid = this->graph_ptr->rank_to_id( comp_ranks[idx] );
              auto srow = gum::util::id_to_charorder( *this->graph_ptr, sid );
              callback( dist_mat, srow, srow );
              if ( info ) {
                info( "Created distance index for region " + std::to_string( idx+1 ) + "." );
              }
            }
          };
          auto nrows = util::total_nof_loci( *this->graph_ptr );
          auto nnz_est = ( nrows - this->graph_ptr->get_node_count() +
                           this->graph_ptr->get_edge_count() ) * ( dmax - dmin );
          mutable_crsmat_type udindex( nrows, nrows, provider, nnz_est );
          this->distance_mat.assign(
              util::compress_distance_index< mutable_crsmat_type >( udindex,
                                                                    *this->graph_ptr ) );
          this->d = std::make_pair( dmin, dmax );
        }

        inline bool
        save_distance_index( std::string prefix ) const
        {
          if ( this->distance_mat.numCols() == 0 ) return true;  // empty distance index

          auto fname = prefix + "_dist_mat_" + "m" + std::to_string( this->d.first ) +
              "M" + std::to_string( this->d.second );
          std::ofstream ofs( fname, std::ofstream::out | std::ofstream::binary );
          if ( !ofs ) return false;

          this->stats_ptr->set_progress( progress_type::write_dindex );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "save-dindex" );

          this->distance_mat.serialize( ofs );
          return true;
        }

        inline bool
        open_distance_index( std::string prefix, unsigned int dmin=0, unsigned int dmax=0 )
        {
          if ( dmax == 0 ) dmax = dmin;
          this->d = std::make_pair( dmin, dmax );

          auto fname = prefix + "_dist_mat_" + "m" + std::to_string( this->d.first ) +
              "M" + std::to_string( this->d.second );
          std::ifstream ifs( fname, std::ifstream::in | std::ifstream::binary );
          if ( !ifs ) return false;

          this->stats_ptr->set_progress( progress_type::load_dindex );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "load-dindex" );

          this->distance_mat.load( ifs );
          return true;
        }

        inline bool
        verify_distance( id_type v, offset_type o, id_type u, offset_type p ) const
        {
          this->stats_ptr->set_progress( progress_type::ready );
          this->stats_ptr->get_this_thread_stats().set_progress(
              thread_progress_type::query_dindex );

          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "query-dindex" );

          if ( v == u ) {  // intra-node distance
            if ( o > p ) return false;
            return this->d.first <= ( p - o ) && ( p - o ) <= this->d.second;
          }
          // inter-node distance
          auto v_charid = gum::util::id_to_charorder( *this->graph_ptr, v ) + o;
          auto u_charid = gum::util::id_to_charorder( *this->graph_ptr, u ) + p;
          return this->distance_mat( v_charid, u_charid );
        }

        /**
         *  @brief  Create path index.
         *
         *  @param  n The number of paths.
         *  @param  patched Whether patch the selected paths or not.
         *  @param  context The context size for patching.
         *  @param  step_size  The step size for starting loci sampling.
         *  @param  dmin  The distance index minimum read insert size.
         *  @param  dmax  The distance index maximum read insert size.
         *  @param  progress A callback function reporting the progress of path selection.
         */
        inline void
        create_path_index( unsigned int n, bool patched=true,
            unsigned int context=0, unsigned int step_size=1,
            unsigned int dmin=0, unsigned int dmax=0,
            std::function< void( std::string const& ) > info=nullptr,
            std::function< void( std::string const& ) > warn=nullptr )
        {
          /* Select the requested number of genome-wide paths. */
          std::function< void( std::string const&, int ) > progress = nullptr;
          if ( info ) {
            progress =
                [&info]( std::string const& name, int i ) {
                  info( "Selecting path " + std::to_string( i ) +
                        " of region " + name + "..." );
                };
          }
          this->pick_paths( n, patched, context, progress, info, warn );
          if ( info ) info( "Indexing the selected paths..." );
          this->index_paths();
          if ( info ) info( "Detecting uncovered loci..." );
          this->add_uncovered_loci( step_size );
          if ( info ) info( "Constructing distance index for pair distance queries..." );
          this->create_distance_index( dmin, dmax, info, warn );
        }

        inline bool
        serialize_path_index_only( std::string const& fpath )
        {
          if ( fpath.empty() ) return false;

          this->stats_ptr->set_progress( progress_type::write_pindex );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "save-pindex" );

          return this->pindex.serialize( fpath );
        }

        /**
         *  @brief  Serialise path index.
         *
         */
        inline bool
        serialize_path_index( std::string const& fpath, unsigned int step_size=1 )
        {
          return this->serialize_path_index_only( fpath ) &&
              this->save_starts( fpath, this->seed_len, step_size ) &&
              this->save_distance_index( fpath );
        }

        /**
         *  @brief  Load path index.
         *
         */
        inline bool
        load_path_index_only( std::string const& fpath, unsigned int context=0 )
        {
          if ( fpath.empty() ) return false;

          this->stats_ptr->set_progress( progress_type::load_pindex );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "load-pindex" );

          this->pindex.set_context( context );
          return this->pindex.load( fpath );
        }

        inline bool
        load_path_index( std::string const& fpath, unsigned int context=0,
                         unsigned int step_size=1, unsigned int dmin=0, unsigned int dmax=0 )
        {
          if ( !this->load_path_index_only( fpath, context ) ) return false;
          if ( !this->open_starts( fpath, this->seed_len, step_size ) ) {
            this->add_uncovered_loci( step_size );
            this->save_starts( fpath, this->seed_len, step_size );
          }
          if ( !this->open_distance_index( fpath, dmin, dmax ) ) {
            this->create_distance_index( dmin, dmax );
            this->save_distance_index( fpath );
          }
          return true;
        }

        /**
         *  @brief  Find seeds on a set of whole-genome paths for the input reads chunk.
         *
         *  @param[in]  reads The set of input reads.
         *  @param[in]  reads_index The reads index.
         *  @param[in]  callback The call back function applied on the found seeds.
         *
         *  This function uses a set of paths from the graph to find seeds of
         *  the input set of reads on these paths by traversing the virtual
         *  suffix tree of both indexes of reads chunk and whole-genome paths.
         */
            inline void
          seeds_on_paths( readsrecord_type const& reads, readsindex_type& reads_index,
                          std::function< void(typename traverser_type::output_type const &) > callback,
                          bool reversed_on_odds=false ) const
          {
            typedef TopDownFine< seqan::ParentLinks<> > TIterSpec;
            typedef typename seqan::Iterator< typename pathindex_type::index_type, TIterSpec >::Type TPIterator;
            typedef typename seqan::Iterator< readsindex_type, TIterSpec >::Type TRIterator;

            auto context = this->pindex.get_context();
            if (  context != 0 /* means patched */ && context < this->seed_len ) {
              throw std::runtime_error( "seed length should not be larger than context size" );
            }

            if ( length( this->pindex.index ) == 0 ) return;

            this->stats_ptr->set_progress( progress_type::ready );
            auto&& thread_stats = this->stats_ptr->get_this_thread_stats();
            thread_stats.set_progress( thread_progress_type::find_on_paths );

            [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "seeds-on-paths" );

            TPIterator piter( this->pindex.index );
            TRIterator riter( reads_index );
            auto collect_stats =
                [&thread_stats]( std::size_t count, bool skipped ) {
                  thread_stats.add_seed_gocc( count );
                  if ( skipped ) thread_stats.inc_gocc_skips();
                };

            kmer_exact_matches( piter, riter, &this->pindex, &reads, this->seed_len,
                                callback, this->gocc_threshold, reversed_on_odds, collect_stats );
          }

          template< typename TString >
          inline void
          seeds_on_paths( TString const& sequence,
                          std::function< void(typename traverser_type::output_type const &) > callback,
                          bool reversed=false ) const
          {
            typedef TopDownFine<> TIterSpec;
            typedef typename seqan::Iterator< typename pathindex_type::index_type, TIterSpec >::Type TPIterator;

            if ( length( indexText( this->pindex.index ) ) == 0 ) return;

            this->stats_ptr->set_progress( progress_type::ready );
            auto&& thread_stats = this->stats_ptr->get_this_thread_stats();
            thread_stats.set_progress( thread_progress_type::find_mems );

            [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "query-paths" );

            TPIterator piter( this->pindex.index );
            auto context = this->pindex.get_context();
            find_mems( sequence, piter, &this->pindex, this->seed_len, context, callback,
                       this->gocc_threshold, reversed );
          }

            inline void
          add_uncovered_loci( unsigned int step=1 )
          {
            auto&& pathset = this->pindex.get_paths_set();
            if ( pathset.size() == 0 )
            {
              this->add_all_loci( step );
              return;
            }

            this->stats_ptr->set_progress( progress_type::find_uncovered );
            [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "find-uncovered" );

            auto bt_itr = begin( *this->graph_ptr, Backtracker() );
            auto bt_end = end( *this->graph_ptr, Backtracker() );
            Path< graph_type > trav_path( this->graph_ptr );
            Path< graph_type > current_path( this->graph_ptr );
            sdsl::bit_vector bv_starts( util::max_node_len( *this->graph_ptr ), 0 );

            this->graph_ptr->for_each_node(
                [&]( rank_type rank, id_type id ) {
                  auto label_len = this->graph_ptr->node_length( id );
                  offset_type offset = label_len;

                  bt_itr.reset( id );
                  while ( bt_itr != bt_end && offset != 0 ) {
                    util::extend_to_k( trav_path, bt_itr, bt_end, offset - 1 + this->seed_len );
                    if ( trav_path.get_sequence_len() >= this->seed_len ) current_path = trav_path;
                    while ( current_path.get_sequence_len() != 0 &&
                        !covered_by( current_path, pathset ) ) {
                      auto trimmed_len = current_path.get_sequence_len()
                        - this->graph_ptr->node_length( current_path.get_nodes().back() );
                      if ( trimmed_len <= this->seed_len - 1 ) {
                        offset = 0;
                        break;
                      }
                      offset = trimmed_len - this->seed_len + 1;
                      trim_back( current_path );
                    }
                    for ( auto f = offset;
                        f < label_len && f + this->seed_len < trav_path.get_sequence_len() + 1;
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
                  return true;
                } );
          }

        inline void add_all_loci( unsigned int step=1 )
        {
          // TODO: Add documentation.
          // TODO: mention in the documentation that the `step` is approximately preserved in
          //       the whole graph.
          this->stats_ptr->set_progress( progress_type::find_uncovered );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "find-uncovered" );

          auto bfs_itr = begin( *this->graph_ptr,  BFS() );
          auto bfs_end = end( *this->graph_ptr,  BFS() );

          unsigned long int prenode_remain = 0;
          unsigned long int remain_estimate = 0;
          id_type prenode_level = 0;
          while ( bfs_itr != bfs_end ) {
            if ( prenode_level != bfs_itr.level() ) {
              prenode_remain = remain_estimate;
              remain_estimate = 0;
              prenode_level = bfs_itr.level();
            }

            auto seq_len = this->graph_ptr->node_length( *bfs_itr );
            unsigned long int cursor = ( step - prenode_remain ) % step;
            while ( cursor < seq_len ) {
              this->add_start( *bfs_itr, cursor );
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

            ++bfs_itr;
          }
        }

        template< typename TPath, typename TSpec >
            inline unsigned long long int
          nof_uncovered_kmers( PathSet< TPath, TSpec >& paths, unsigned int k )
          {
            if ( this->starting_loci.size() == 0 ) return 0;
            [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "count-uncovered-kmer" );

            auto bt_itr = begin( *this->graph_ptr, Backtracker() );
            auto bt_end = end( *this->graph_ptr, Backtracker() );
            Path< graph_type > trav_path( this->graph_ptr );
            Path< graph_type > current_path( this->graph_ptr );
            unsigned long long int uncovered = 0;

            long long int prev_id = 0;
            for ( const vg::Position& l : this->starting_loci ) {
              if ( prev_id == l.node_id() ) continue;
              prev_id = l.node_id();
              auto label_len = this->graph_ptr->node_length( l.node_id() );

              bt_itr.reset( l.node_id() );
              while ( bt_itr != bt_end ) {
                offset_type offset = label_len;
                util::extend_to_k( trav_path, bt_itr, bt_end, offset - 1 + k );
                if ( trav_path.get_sequence_len() >= k ) current_path = trav_path;
                while ( current_path.get_sequence_len() != 0 &&
                    !covered_by( current_path, paths ) ) {
                  auto trimmed_len = current_path.get_sequence_len()
                    - this->graph_ptr->node_length( current_path.get_nodes().back() );
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

          this->stats_ptr->set_progress( progress_type::load_starts );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "load-starts" );

          std::function< void( vg::Position& ) > push_back =
              [this]( vg::Position& pos ) {
                pos.set_node_id( this->graph_ptr->id_by_coordinate( pos.node_id() ) );
                this->starting_loci.push_back( pos );
              };

          try {
            vg::io::for_each( ifs, push_back );
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

          this->stats_ptr->set_progress( progress_type::write_starts );
          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "save-starts" );

          std::function< vg::Position( uint64_t ) > lambda =
              [this]( uint64_t i ) {
                auto pos = this->starting_loci.at( i );
                pos.set_node_id( this->graph_ptr->coordinate_id( pos.node_id() ) );
                return pos;
              };

          try {
            vg::io::write( ofs, this->starting_loci.size(), lambda );
          }
          catch ( const std::runtime_error& ) {
            return false;
          }

          return true;
        }

          inline std::size_t
        get_nof_uniq_nodes( )
        {
          std::unordered_set< id_type > set;
          for ( const auto& l : this->starting_loci ) set.insert( l.node_id() );
          return set.size();
        }

          inline traverser_type
        create_traverser( bool reversed_on_odds=false ) const
        {
          return traverser_type( this->graph_ptr, this->seed_len, reversed_on_odds );
        }

          inline void
        setup_traverser( traverser_type& traverser, readsrecord_type const& reads,
                         readsindex_type& reads_index ) const
        {
          traverser.set_reads( &reads );
          traverser.set_reads_index( &reads_index );
        }

          inline void
        seeds_off_paths( traverser_type& traverser,
                         std::function< void( typename traverser_type::output_type const& ) > callback ) const
        {
          this->stats_ptr->set_progress( progress_type::ready );
          this->stats_ptr->get_this_thread_stats().set_progress(
              thread_progress_type::find_off_paths );

          [[maybe_unused]] auto timer = this->stats_ptr->timeit_ts( "seeds-off-path" );

          for ( std::size_t idx = 0; idx < this->starting_loci.size(); ++idx )
          {
            const auto& locus = this->starting_loci[ idx ];
            traverser.add_locus( locus );
            if ( this->starting_loci[ idx + 1 ].node_id() == locus.node_id() ) continue;

            traverser.run( callback );
            this->stats_ptr->get_this_thread_stats().set_locus_idx( idx );
          }
        }

          inline void
        seeds_all( readsrecord_type const& reads, readsindex_type& reads_index, traverser_type& traverser,
                   std::function< void(typename traverser_type::output_type const &) > callback,
                   bool reversed_on_odds=false ) const
        {
          this->seeds_on_paths( reads, reads_index, callback, reversed_on_odds );
          this->setup_traverser( traverser, reads, reads_index );
          this->seeds_off_paths( traverser, callback );
          this->stats_ptr->get_this_thread_stats().inc_chunks_done( );
        }

          inline void
        seeds_all( readsrecord_type const& reads, readsindex_type& reads_index, traverser_type& traverser,
                   std::function< void(typename traverser_type::output_type const &) > callback1,
                   std::function< void(typename traverser_type::output_type const &) > callback2,
                   bool reversed_on_odds=false ) const
        {
          this->seeds_on_paths( reads, reads_index, callback1, reversed_on_odds );
          this->setup_traverser( traverser, reads, reads_index );
          this->seeds_off_paths( traverser, callback2 );
          this->stats_ptr->get_this_thread_stats().inc_chunks_done( );
        }

      private:
        /* ====================  DATA MEMBERS  ======================================= */
        const graph_type* graph_ptr;
        std::vector< vg::Position > starting_loci;
        pathindex_type pindex;  /**< @brief Genome-wide path index in lazy mode. */
        KokkosHandler handler;
        crsmat_type distance_mat;
        unsigned int seed_len;
        unsigned char seed_mismatches;  /**< @brief Allowed mismatches in a seed hit. */
        unsigned int gocc_threshold;  /**< @brief Seed genome occurrence count threshold. */
        std::pair< unsigned int, unsigned int > d; /**< @brief distance constraints. */
        std::unique_ptr< stats_type > stats_ptr;
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  Set the context size for patching.
         *
         *  @param  context The context size for patching.
         *  @param  patched Whether patch the selected paths or not.
         *
         *  @return the actual context length set for the path index.
         *
         *  NOTE: If `patched` is set, the context size cannot be zero. If the
         *  given value for context size is zero, it will be reset to seed
         *  length value. Otherwise (`patched` is `false`), the context doesn't
         *  matter and will be set to zero.
         */
        inline unsigned int
        set_context( unsigned int context, bool patched,
            std::function< void( std::string const& ) > info=nullptr,
            std::function< void( std::string const& ) > warn=nullptr )
        {
          /* See the NOTE section above. */
          if ( !patched ) context = 0;
          if ( patched && context == 0 ) {
            if ( warn ) warn( "The context size cannot be zero for patching. "
                              "Assuming the seed length as the context size..." );
            context = this->seed_len;
          }
          /* Set the context size for path index. */
          this->pindex.set_context( context );
          return context;
        }
    };
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_SEED_FINDER_HPP__ --- */
