/**
 *    @file  stat.h
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

#ifndef STAT_H__
#define STAT_H__

#include <chrono>
#include <unordered_map>

#include <seqan/basic.h>


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

  /**
   *  @brief  Timers for measuring execution time.
   *
   *  Measure the time period between its instantiation and destruction. The timers are
   *  kept in static table hashed by the timer name.
   */
  class Timer
  {
    public:
      /* ====================  TYPEDEFS      ======================================= */
      typedef std::chrono::steady_clock clock_type;
      /* ====================  CLASSES       ======================================= */
      struct TimePeriod {
        clock_type::time_point start;
        clock_type::time_point end;
      };
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
       *  @brief  Get the timer by name.
       *
       *  @param  name The name of the timer to start.
       *  @return the duration represented by requested timer in microsecond.
       *
       *  Get the duration represented by the timer in microsecond.
       */
        static inline std::chrono::microseconds
      get_duration( const std::string& name )
      {
        auto duration = get_timers()[ name ].end - get_timers()[ name ].start;
        return std::chrono::duration_cast< std::chrono::microseconds >( duration );
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
       *  returns the duration between start time to 'now' in microseconds.
       */
        static inline std::chrono::microseconds
      get_lap( const std::string& name )
      {
        if ( get_timers()[ name ].end > get_timers()[ name ].start ) {
          return get_duration( name );
        }
        auto duration = clock_type::now() - get_timers()[ name ].start;
        return std::chrono::duration_cast< std::chrono::microseconds >( duration );
      }  /* -----  end of method get_lap  ----- */
    protected:
      /* ====================  DATA MEMBERS  ======================================= */
      std::string timer_name;    /**< @brief The timer name of the current instance. */
  };  /* -----  end of class Timer  ----- */

}  /* -----  end of namespace grem  ----- */

#endif  // end of STAT_H__
