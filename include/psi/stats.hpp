/**
 *    @file  stats.hpp
 *   @brief  Module for performance measurement and running-time statistics.
 *
 *  This module provides necessary tools for measuring performance and capturing
 *  running-time statistics.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Aug 31, 2017  00:40
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_STATS_HPP__
#define PSI_STATS_HPP__

#include <ctime>
#include <chrono>
#include <unordered_map>
#include <cassert>

#include "base.hpp"


namespace psi {
  struct NoStats;                        /**< @brief No-stats mode tag. */
  struct WithStats;                      /**< @brief With-stats mode tag. */

  /**
   *  @brief  Observer template class to collect statistics.
   *
   *  This class is an observer class to collect some statistics from host class; i.e.
   *  `TObject` template parameter.
   */
  template < typename TObject >
    class Stats;

  template< class TObject >
    using StatsType = typename Stats< TObject >::type;

  typedef clock_t CpuClock;
  typedef std::chrono::steady_clock SteadyClock;
  typedef struct timespec BetterClock;

  template< typename TSpec = CpuClock >
    class TimerTraits;

  template< >
    class TimerTraits< std::chrono::steady_clock > {
      public:
        /* ====================  TYPE MEMBERS  ======================================= */
        typedef std::chrono::steady_clock clock_type;
        typedef std::chrono::microseconds duration_type;
        typedef duration_type::rep rep_type;
        /* ====================  DATA MEMBERS  ======================================= */
        constexpr static const char* unit_repr = "us";
        constexpr static const duration_type zero_duration =
          std::chrono::duration_cast< duration_type >( clock_type::duration::zero() );
        constexpr static const rep_type zero_duration_rep =
          TimerTraits::zero_duration.count();
        /* ====================  METHODS       ======================================= */
          static inline duration_type
        duration( clock_type::time_point end, clock_type::time_point start,
                  duration_type pre_elapsed=zero_duration )
        {
          return std::chrono::duration_cast< duration_type >( end - start ) + pre_elapsed;
        }

          static inline rep_type
        duration_rep( clock_type::time_point end, clock_type::time_point start,
                      duration_type pre_elapsed=zero_duration )
        {
          return TimerTraits::duration( end, start, pre_elapsed ).count();
        }

          static inline std::string
        duration_str( clock_type::time_point end, clock_type::time_point start,
                      duration_type pre_elapsed=zero_duration )
        {
          return std::to_string( TimerTraits::duration_rep( end, start, pre_elapsed ) ) +
              " " + TimerTraits::unit_repr;
        }
    };  /* ---  end of class TimerTraits  --- */

  template< >
    class TimerTraits< struct timespec > {
      public:
        /* ====================  TYPE MEMBERS  ======================================= */
        struct Clock {
          struct TimePoint {
            typedef double duration_type;
            constexpr static const long int NSEC_PER_USEC = 1000;
            constexpr static const long int NSEC_PER_SEC = 1000000000;
            constexpr static const long int USEC_PER_SEC = NSEC_PER_SEC / NSEC_PER_USEC;

            TimePoint() = default;

            TimePoint( duration_type d )
            {
              this->p.tv_sec = d / USEC_PER_SEC;
              this->p.tv_nsec = ( d - this->p.tv_sec * USEC_PER_SEC ) * NSEC_PER_USEC;
            }

            inline
            operator duration_type() const
            {
              return this->p.tv_sec * USEC_PER_SEC + this->p.tv_nsec / NSEC_PER_USEC;
            }

              inline bool
            operator>( TimePoint const& y ) const
            {
              if ( this->p.tv_sec == y.p.tv_sec ) return this->p.tv_nsec > y.p.tv_nsec;
              return this->p.tv_sec > y.p.tv_sec;
            }

              inline bool
            operator<( TimePoint const& y ) const
            {
              return y > *this;
            }

            /**
             *  @brief Subtract the ‘struct timespec’ values x and y.
             *
             *  NOTE: x should be larger than y!
             */
              inline TimePoint
            operator-( TimePoint y ) const
            {
              TimePoint result;
              /* Perform the carry for the later subtraction by updating y. */
              if ( this->p.tv_nsec < y.p.tv_nsec ) {
                int sec = ( y.p.tv_nsec - this->p.tv_nsec ) / NSEC_PER_SEC + 1;
                y.p.tv_nsec -= sec * NSEC_PER_SEC;
                y.p.tv_sec += sec;
              }
              if ( this->p.tv_nsec - y.p.tv_nsec > NSEC_PER_SEC ) {
                int sec = ( this->p.tv_nsec - y.p.tv_nsec ) / NSEC_PER_SEC;
                y.p.tv_nsec += sec * NSEC_PER_SEC;
                y.p.tv_sec -= sec;
              }
              /* Compute the time remaining to wait. tv_nsec is certainly positive. */
              result.p.tv_sec = this->p.tv_sec - y.p.tv_sec;
              result.p.tv_nsec = this->p.tv_nsec - y.p.tv_nsec;

              return result;
            }
            struct timespec p;
          };
          typedef TimePoint time_point;
            static inline time_point
          now( clockid_t cid=0 )
          {
            time_point tp;
            if ( cid == 0 ) clock_getcpuclockid( pthread_self(), &cid );
            clock_gettime( cid, &tp.p );
            return tp;
          }
        };
        typedef Clock clock_type;
        typedef clock_type::time_point::duration_type duration_type;
        typedef duration_type rep_type;
        /* ====================  DATA MEMBERS  ======================================= */
        constexpr static const char* unit_repr = "us";
        constexpr static const duration_type zero_duration = 0;
        constexpr static const rep_type zero_duration_rep = 0;
        /* ====================  METHODS       ======================================= */
          static inline duration_type
        duration( clock_type::time_point end, clock_type::time_point start,
                  duration_type pre_elapsed=zero_duration )
        {
          return static_cast< duration_type >( end - start ) + pre_elapsed;
        }

          static inline rep_type
        duration_rep( clock_type::time_point end, clock_type::time_point start,
                      duration_type pre_elapsed=zero_duration )
        {
          return TimerTraits::duration( end, start, pre_elapsed );
        }

          static inline std::string
        duration_str( clock_type::time_point end, clock_type::time_point start,
                      duration_type pre_elapsed=zero_duration )
        {
          return std::to_string( TimerTraits::duration_rep( end, start, pre_elapsed ) ) +
              " " + TimerTraits::unit_repr;
        }
    };  /* ---  end of class TimerTraits  --- */

  template< >
    class TimerTraits< clock_t > {
      public:
        /* ====================  TYPE MEMBERS  ======================================= */
        typedef struct {
          typedef clock_t time_point;
            static inline time_point
          now()
          {
            return clock();
          }
        } clock_type;
        typedef float duration_type;
        typedef duration_type rep_type;
        /* ====================  DATA MEMBERS  ======================================= */
        constexpr static const char* unit_repr = "s";
        constexpr static const duration_type zero_duration = 0;
        constexpr static const rep_type zero_duration_rep = 0;
        /* ====================  METHODS       ======================================= */
          static inline duration_type
        duration( clock_type::time_point end, clock_type::time_point start,
                  duration_type pre_elapsed=zero_duration )
        {
          return static_cast< float >( end - start ) / CLOCKS_PER_SEC + pre_elapsed;
        }

          static inline rep_type
        duration_rep( clock_type::time_point end, clock_type::time_point start,
                      duration_type pre_elapsed=zero_duration )
        {
          return TimerTraits::duration( end, start, pre_elapsed );
        }

          static inline std::string
        duration_str( clock_type::time_point end, clock_type::time_point start,
                      duration_type pre_elapsed=zero_duration )
        {
          return std::to_string( TimerTraits::duration_rep( end, start, pre_elapsed ) ) +
              " " + TimerTraits::unit_repr;
        }
    };  /* ---  end of class TimerTraits  --- */

  template< >
    class TimerTraits< void > {
      public:
        /* === TYPE MEMBERS === */
        typedef struct {
          typedef void* time_point;
          static inline time_point now() { return nullptr; }
        } clock_type;
        typedef float duration_type;
        typedef duration_type rep_type;
        /* === DATA MEMBERS === */
        constexpr static const char* unit_repr = "s";
        constexpr static const duration_type zero_duration = 0;
        constexpr static const rep_type zero_duration_rep = 0;
        /* === METHODS === */
        constexpr static inline duration_type
        duration( clock_type::time_point, clock_type::time_point )
        {
          return TimerTraits::zero_duration;
        }

        constexpr static inline duration_type
        duration( clock_type::time_point, clock_type::time_point, duration_type )
        {
          return TimerTraits::zero_duration;
        }

        constexpr static inline rep_type
        duration_rep( clock_type::time_point, clock_type::time_point )
        {
          return TimerTraits::zero_duration_rep;
        }

        constexpr static inline rep_type
        duration_rep( clock_type::time_point, clock_type::time_point, duration_type )
        {
          return TimerTraits::zero_duration_rep;
        }

        constexpr static inline const char*
        duration_str( clock_type::time_point, clock_type::time_point )
        {
          return "0";
        }

        constexpr static inline const char*
        duration_str( clock_type::time_point, clock_type::time_point, duration_type )
        {
          return "0";
        }
    };  /* ---  end of class TimerTraits  --- */

  /**
   *  @brief  Timers for measuring execution time.
   *
   *  Measure the time period between its instantiation and destruction. The timers are
   *  kept in static table hashed by the timer name.
   */
  template< typename TClock = CpuClock >
    class Timer
    {
      public:
        /* ====================  TYPE MEMBERS  ======================================= */
        typedef TimerTraits< TClock > traits_type;
        typedef typename traits_type::clock_type clock_type;
        typedef typename traits_type::duration_type duration_type;
        typedef typename traits_type::rep_type rep_type;
        struct TimePeriod {
          typedef Timer timer_type;

          duration_type pre_elapsed;
          typename clock_type::time_point start;
          typename clock_type::time_point end;
          TimePeriod( ) : pre_elapsed( traits_type::zero_duration ),
                          start( traits_type::zero_duration ),
                          end( traits_type::zero_duration )
          { }

          /**
           *  @brief  Get the measured duration.
           */
            inline duration_type
          duration( ) const
          {
            return traits_type::duration( this->end, this->start, this->pre_elapsed );
          }

          /**
           *  @brief  Get the measured duration (arithmetic representation).
           */
            inline rep_type
          rep( ) const
          {
            return traits_type::duration_rep( this->end, this->start, this->pre_elapsed );
          }

          /**
           *  @brief  Get the measured duration (string representation).
           */
            inline std::string
          str( ) const
          {
            return traits_type::duration_str( this->end, this->start, this->pre_elapsed );
          }

          /**
           *  @brief  Get time lap for an ongoing timer.
           */
            inline TimePeriod
          get_lap( ) const
          {
            TimePeriod period( *this );
            if ( period.end <= period.start ) period.end = clock_type::now();
            return period;
          }
        };
        typedef TimePeriod period_type;
        /* ====================  STATIC DATA   ======================================= */
        constexpr static const char* unit_repr = traits_type::unit_repr;
        constexpr static const duration_type zero_duration = traits_type::zero_duration;
        constexpr static const rep_type zero_duration_rep = traits_type::zero_duration_rep;
        /* ====================  LIFECYCLE     ======================================= */
        /**
         *  @brief  Timer constructor.
         *
         *  @param  name The name of the timer to start.
         *
         *  If timer does not exist it will be created.
         */
        Timer( const std::string& name )
        {
          this->timer_name = name;
          auto found = get_timers().find( this->timer_name );
          if ( found != get_timers().end() ) {
            assert( found->second.end > found->second.start );
            found->second.pre_elapsed = found->second.duration();
            found->second.start = clock_type::now();
          } else get_timers()[ this->timer_name ].start = clock_type::now();
        }  /* -----  end of method Timer  (constructor)  ----- */

        /**
         *  @brief  Timer destructor.
         *
         *  Stop the timer as the Timer object dies.
         */
        ~Timer()
        {
          get_timers()[ this->timer_name ].end = clock_type::now();
        }  /* -----  end of method ~Timer  (destructor)  ----- */
        /* ====================  METHODS       ======================================= */
        /**
         *  @brief  static getter function for static timers.
         */
          static inline std::unordered_map< std::string, TimePeriod >&
        get_timers( )
        {
          static std::unordered_map< std::string, TimePeriod > timers;
          return timers;
        }  /* -----  end of method get_timers  ----- */

        /**
         *  @brief  Get the timer duration by name.
         *
         *  @param  name The name of the timer to start.
         *  @return the duration represented by requested timer.
         *
         *  Get the duration represented by the timer.
         */
          static inline duration_type
        get_duration( const std::string& name )
        {
          return get_timers()[ name ].duration();
        }  /* -----  end of method get_duration  ----- */

        /**
         *  @brief  Get the timer duration (arithmetic representation) by name.
         *
         *  @param  name The name of the timer to start.
         *  @return arithmetic representation of the requested timer duration.
         *
         *  Get the arithmetic representation of the requested timer duration.
         */
          static inline rep_type
        get_duration_rep( const std::string& name )
        {
          return get_timers()[ name ].rep();
        }  /* -----  end of method get_duration  ----- */

        /**
         *  @brief  Get the timer duration (string representation) by name.
         *
         *  @param  name The name of the timer to start.
         *  @return string representation of the requested timer duration.
         *
         *  Get the string representation of the requested timer duration.
         */
          static inline std::string
        get_duration_str( const std::string& name )
        {
          return get_timers()[ name ].str();
        }  /* -----  end of method get_duration  ----- */

        /**
         *  @brief  Get time lap for an ongoing timer.
         *
         *  @param  name The name of the timer.
         *  @return the duration since 'start' to 'now' if the timer is not finished;
         *          otherwise it returns the duration of the timer.
         */
          static inline duration_type
        get_lap_duration( const std::string& name )
        {
          return get_timers()[ name ].get_lap().duration();
        }  /* -----  end of method get_lap  ----- */

        /**
         *  @brief  Get time lap for an ongoing timer (arithmetic representation).
         *
         *  @param  name The name of the timer.
         *  @return the duration since 'start' to 'now' if the timer is not finished;
         *          otherwise it returns the duration of the timer.
         */
          static inline rep_type
        get_lap_rep( const std::string& name )
        {
          return get_timers()[ name ].get_lap().rep();
        }  /* -----  end of method get_lap  ----- */

        /**
         *  @brief  Get time lap for an ongoing timer (string representation).
         *
         *  @param  name The name of the timer.
         *  @return the duration since 'start' to 'now' if the timer is not finished;
         *          otherwise it returns the duration of the timer.
         */
          static inline std::string
        get_lap_str( const std::string& name )
        {
          return get_timers()[ name ].get_lap().str();
        }  /* -----  end of method get_lap  ----- */
      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        std::string timer_name;    /**< @brief The timer name of the current instance. */
    };  /* ---  end of class Timer  --- */

  template< >
    class Timer< void >
    {
    public:
      /* === TYPE MEMBERS === */
      typedef TimerTraits< void > traits_type;
      typedef typename traits_type::clock_type clock_type;
      typedef typename traits_type::duration_type duration_type;
      typedef typename traits_type::rep_type rep_type;
      struct TimePeriod {
        typedef Timer timer_type;

          constexpr inline duration_type
        duration( ) const
        {
          return traits_type::zero_duration;
        }

          constexpr inline rep_type
        rep( ) const
        {
          return traits_type::zero_duration_rep;
        }

          constexpr inline const char*
        str( ) const
        {
          return "0";
        }

          constexpr inline TimePeriod
        get_lap( ) const
        {
          return TimePeriod();
        }
      };
      typedef TimePeriod period_type;
      /* === LIFECYCLE === */
      Timer( std::string const& ) { }
      Timer( ) { }
      /* === METHODS === */
      static inline std::unordered_map< std::string, TimePeriod >
      get_timers( )
      {
        return std::unordered_map< std::string, TimePeriod >{};
      }
      constexpr static inline duration_type get_duration( const std::string& )
      {
        return traits_type::zero_duration;
      }
      constexpr static inline rep_type get_duration_rep( const std::string& )
      {
        return traits_type::zero_duration_rep;
      }
      constexpr static inline const char* get_duration_str( const std::string& )
      {
        return "0";
      }
      constexpr static inline duration_type get_lap_duration( const std::string& )
      {
        return traits_type::zero_duration;
      }
      constexpr static inline rep_type get_lap_rep( const std::string& )
      {
        return traits_type::zero_duration_rep;
      }
      constexpr static inline const char* get_lap_str( const std::string& )
      {
        return "0";
      }
    };  /*  --- end of class Timer --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_STATS_HPP__ --- */
