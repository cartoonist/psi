/**
 *    @file  stat.hpp
 *   @brief  Stat template class.
 *
 *  Stat template class definition.
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

#ifndef PSI_STAT_HPP__
#define PSI_STAT_HPP__

#include <chrono>
#include <unordered_map>

#include <seqan/basic.h>

#include "base.hpp"


namespace grem {
  struct NoStatStrategy;                        /**< @brief No-stat mode strategy. */
  typedef seqan::Tag< NoStatStrategy > NoStat;  /**< @brief No-stat mode tag. */
  /**
   *  @brief  Observer template class to collect statistics.
   *
   *  This class is an observer class to collect some statistics from host class; i.e.
   *  `TObject` template parameter.
   */
  template < typename TObject >
    class Stat;

  template< typename TSpec >
    class TimerTraits;

  typedef clock_t CpuClock;
  typedef std::chrono::steady_clock SteadyClock;

  template< >
    class TimerTraits< std::chrono::steady_clock > {
      public:
        /* ====================  MEMBER TYPES  ======================================= */
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
    };

  template< >
    class TimerTraits< clock_t > {
      public:
        /* ====================  MEMBER TYPES  ======================================= */
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
    };

  /**
   *  @brief  Timers for measuring execution time.
   *
   *  Measure the time period between its instantiation and destruction. The timers are
   *  kept in static table hashed by the timer name.
   */
  template< typename TClock=clock_t >
    class Timer
    {
      public:
        /* ====================  MEMBER TYPES  ======================================= */
        typedef TimerTraits< TClock > trait_type;
        typedef typename trait_type::clock_type clock_type;
        typedef typename trait_type::duration_type duration_type;
        typedef typename trait_type::rep_type rep_type;
        struct TimePeriod {
          typename clock_type::time_point start;
          typename clock_type::time_point end;
        };
        /* ====================  STATIC DATA   ======================================= */
        constexpr static const char* unit_repr = trait_type::unit_repr;
        constexpr static const duration_type zero_duration = trait_type::zero_duration;
        constexpr static const rep_type zero_duration_rep = trait_type::zero_duration_rep;
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
        get_duration( const std::string& name )
        {
          return trait_type::duration( get_timers()[ name ].end,
              get_timers()[ name ].start );
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
          return trait_type::duration_rep( get_timers()[ name ].end,
              get_timers()[ name ].start );
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
          return trait_type::duration_str( get_timers()[ name ].end,
              get_timers()[ name ].start );
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
        get_lap( const std::string& name )
        {
          if ( get_timers()[ name ].end > get_timers()[ name ].start ) {
            return get_duration( name );
          }
          return trait_type::duration( clock_type::now(), get_timers()[ name ].start );
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
        get_lap_rep( const std::string& name )
        {
          if ( get_timers()[ name ].end > get_timers()[ name ].start ) {
            return get_duration_rep( name );
          }
          return trait_type::duration_rep( clock_type::now(), get_timers()[ name ].start );
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
        get_lap_str( const std::string& name )
        {
          if ( get_timers()[ name ].end > get_timers()[ name ].start ) {
            return get_duration_str( name );
          }
          return trait_type::duration_str( clock_type::now(), get_timers()[ name ].start );
        }  /* -----  end of method get_lap  ----- */
      protected:
        /* ====================  DATA MEMBERS  ======================================= */
        std::string timer_name;    /**< @brief The timer name of the current instance. */
    };  /* -----  end of class Timer  ----- */

}  /* -----  end of namespace grem  ----- */

#endif  /* --- #ifndef PSI_STAT_HPP__ --- */
