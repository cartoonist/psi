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

#include <chrono>
#include <unordered_map>

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
    class Stat;

  typedef clock_t CpuClock;
  typedef std::chrono::steady_clock SteadyClock;

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
        duration( clock_type::time_point end, clock_type::time_point start )
        {
          return std::chrono::duration_cast< duration_type >( end - start );
        }

          static inline rep_type
        duration_rep( clock_type::time_point end, clock_type::time_point start )
        {
          return TimerTraits::duration( end, start ).count();
        }

          static inline std::string
        duration_str( clock_type::time_point end, clock_type::time_point start )
        {
          return
            std::to_string( TimerTraits::duration_rep( end, start ) ) + " "
            + TimerTraits::unit_repr;
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
        duration( clock_type::time_point end, clock_type::time_point start )
        {
          return static_cast< float >( end - start ) / CLOCKS_PER_SEC;
        }

          static inline rep_type
        duration_rep( clock_type::time_point end, clock_type::time_point start )
        {
          return TimerTraits::duration( end, start );
        }

          static inline std::string
        duration_str( clock_type::time_point end, clock_type::time_point start )
        {
          return std::to_string( TimerTraits::duration_rep( end, start ) ) + " "
            + TimerTraits::unit_repr;
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

        constexpr static inline rep_type
        duration_rep( clock_type::time_point, clock_type::time_point )
        {
          return TimerTraits::zero_duration_rep;
        }

        constexpr static inline const char*
        duration_str( clock_type::time_point, clock_type::time_point )
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
  template< typename TClock=CpuClock >
    class Timer
    {
      public:
        /* ====================  TYPE MEMBERS  ======================================= */
        typedef TimerTraits< TClock > traits_type;
        typedef typename traits_type::clock_type clock_type;
        typedef typename traits_type::duration_type duration_type;
        typedef typename traits_type::rep_type rep_type;
        struct TimePeriod {
          typename clock_type::time_point start;
          typename clock_type::time_point end;
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
          get_timers()[ this->timer_name ].start = clock_type::now();
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
        get_duration( TimePeriod const& period )
        {
          return traits_type::duration( period.end, period.start );
        }

          static inline duration_type
        get_duration( const std::string& name )
        {
          return get_duration( get_timers()[ name ] );
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
        get_duration_rep( TimePeriod const& period )
        {
          return traits_type::duration_rep( period.end, period.start );
        }

          static inline rep_type
        get_duration_rep( const std::string& name )
        {
          return get_duration_rep( get_timers()[ name ] );
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
        get_duration_str( TimePeriod const& period )
        {
          return traits_type::duration_str( period.end, period.start );
        }

          static inline std::string
        get_duration_str( const std::string& name )
        {
          return get_duration_str( get_timers()[ name ] );
        }  /* -----  end of method get_duration  ----- */

        /**
         *  @brief  Get time lap for an ongoing timer.
         *
         *  @param  name The name of the timer.
         *  @return the duration since 'start' to 'now' if the timer is not finished;
         *          otherwise it returns the duration of the timer by calling method
         *          `get_duration`.
         *
         *  It first checks if the timer specified by its `name` is finished or not. If
         *  so, it return the duration by calling `get_duration` method. Otherwise, it
         *  returns the duration between start time to 'now'.
         */
          static inline duration_type
        get_lap( TimePeriod const& period )
        {
          if ( period.end > period.start ) {
            return get_duration( period );
          }
          return traits_type::duration( clock_type::now(), period.start );
        }

          static inline duration_type
        get_lap( const std::string& name )
        {
          return get_lap( get_timers()[ name ] );
        }  /* -----  end of method get_lap  ----- */

        /**
         *  @brief  Get time lap for an ongoing timer.
         *
         *  @param  name The name of the timer.
         *  @return the duration since 'start' to 'now' if the timer is not finished;
         *          otherwise it returns the duration of the timer by calling method
         *          `get_duration`.
         *
         *  It first checks if the timer specified by its `name` is finished or not. If
         *  so, it return the duration by calling `get_duration` method. Otherwise, it
         *  returns the duration between start time to 'now'.
         */
          static inline rep_type
        get_lap_rep( TimePeriod const& period )
        {
          if ( period.end > period.start ) {
            return get_duration_rep( period );
          }
          return traits_type::duration_rep( clock_type::now(), period.start );
        }

          static inline rep_type
        get_lap_rep( const std::string& name )
        {
          return get_lap_rep( get_timers()[ name ] );
        }  /* -----  end of method get_lap  ----- */

        /**
         *  @brief  Get time lap for an ongoing timer.
         *
         *  @param  name The name of the timer.
         *  @return the duration since 'start' to 'now' if the timer is not finished;
         *          otherwise it returns the duration of the timer by calling method
         *          `get_duration`.
         *
         *  It first checks if the timer specified by its `name` is finished or not. If
         *  so, it return the duration by calling `get_duration` method. Otherwise, it
         *  returns the duration between start time to 'now'.
         */
          static inline std::string
        get_lap_str( TimePeriod const& period )
        {
          if ( period.end > period.start ) {
            return get_duration_str( period );
          }
          return traits_type::duration_str( clock_type::now(), period.start );
        }

          static inline std::string
        get_lap_str( const std::string& name )
        {
          return get_lap_str( get_timers()[ name ] );
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
      typedef void* period_type;
      /* === LIFECYCLE === */
      Timer( std::string const& ) { }
      Timer( ) { }
      /* === METHODS === */
      constexpr static inline traits_type::duration_type get_duration( period_type )
      {
        return traits_type::zero_duration;
      }
      constexpr static inline traits_type::duration_type get_duration( const std::string& )
      {
        return traits_type::zero_duration;
      }
      constexpr static inline traits_type::rep_type get_duration_rep( period_type )
      {
        return traits_type::zero_duration_rep;
      }
      constexpr static inline traits_type::rep_type get_duration_rep( const std::string& )
      {
        return traits_type::zero_duration_rep;
      }
      constexpr static inline const char* get_duration_str( period_type )
      {
        return "0";
      }
      constexpr static inline const char* get_duration_str( const std::string& )
      {
        return "0";
      }
      constexpr static inline traits_type::duration_type get_lap( period_type )
      {
        return traits_type::zero_duration;
      }
      constexpr static inline traits_type::duration_type get_lap( const std::string& )
      {
        return traits_type::zero_duration;
      }
      constexpr static inline traits_type::rep_type get_lap_rep( period_type )
      {
        return traits_type::zero_duration_rep;
      }
      constexpr static inline traits_type::rep_type get_lap_rep( const std::string& )
      {
        return traits_type::zero_duration_rep;
      }
      constexpr static inline const char* get_lap_str( period_type )
      {
        return "0";
      }
      constexpr static inline const char* get_lap_str( const std::string& )
      {
        return "0";
      }
    };  /*  --- end of class Timer --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_STATS_HPP__ --- */
