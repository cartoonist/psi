/**
 *    @file  logger.hpp
 *   @brief  Logger header file.
 *
 *  This header file wraps the external logging library and sets up the general
 *  configuration of the library required for the entire package.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Feb 26, 2017  12:59
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_LOGGER_HPP__
#define PSI_LOGGER_HPP__

#include <spdlog/spdlog.h>
#include <spdlog/async.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "options.hpp"


namespace psi {
  /**
   *  @brief  Configure the main logger.
   *
   *  @param  nolog Disable logging completely.
   *  @param  quiet Disable logging to console.
   *  @param  nocolor Disable coloured output.
   *  @param  verbose Enable verbosity.
   *  @param  nologfile Disable logging to file.
   *  @param  log_path Path to log file.
   *
   *  It configures the "main" logger. The logger can be accessed by `get_logger`
   *  method.
   *
   *  NOTE: The logger are multi-threaded.
   */
    inline void
  config_logger( bool nolog, bool quiet, bool nocolor, bool verbose, bool nologfile,
     const std::string& log_path )
  {
    typedef spdlog::sinks::stdout_sink_mt console_sink;
    typedef spdlog::sinks::stdout_color_sink_mt console_color_sink;
    typedef spdlog::sinks::basic_file_sink_mt file_sink;
    typedef spdlog::async_logger logger_type;

    std::vector< spdlog::sink_ptr > sinks;

    if ( !nolog && !quiet ) {
      if ( nocolor ) sinks.push_back( std::make_shared< console_sink >() );
      else sinks.push_back( std::make_shared< console_color_sink >() );

      if ( verbose ) sinks.back()->set_level( spdlog::level::info );
      else sinks.back()->set_level( spdlog::level::warn );
    }

    if ( !nolog && !nologfile ) {
      sinks.push_back( std::make_shared< file_sink >( log_path ) );
      // Always verbose for file sink.
      sinks.back()->set_level( spdlog::level::info );
    }

    spdlog::init_thread_pool(8192, 2);
    auto main_logger =
        std::make_shared< logger_type >( "main",
                                         begin(sinks),
                                         end(sinks),
                                         spdlog::thread_pool(),
                                         spdlog::async_overflow_policy::overrun_oldest );
    spdlog::register_logger( main_logger );
  }

  /**
   *  @brief  Configure the main logger.
   *
   *  @param  options Program options.
   *
   *  It configures the "main" logger. The logger can be accessed by `get_logger`
   *  method.
   *
   *  NOTE: The logger are multi-threaded.
   */
    inline void
  config_logger( const Options& options )
  {
    config_logger( options.nolog, options.quiet, options.nocolor, options.verbose,
        options.nologfile, options.log_path );
  }

  /**
   *  @brief  Get the logger by name.
   *
   *  @param  name The name of the requested logger.
   *  @return shared pointer to the logger.
   *
   *  It is a wrapper function for `spdlog::get` method.
   */
    inline std::shared_ptr<spdlog::logger>
  get_logger( const std::string& name )
  {
    return spdlog::get( name );
  }  /* -----  end of function get_logger  ----- */

  /**
   *  @brief  Release and close a logger by name.
   *
   *  @param  name The name of the logger.
   *
   *  It is a wrapper function for `spdlog::drop` method.
   */
    inline void
  drop_logger( const std::string& name )
  {
    spdlog::drop( name );
  }  /* -----  end of function drop_loggers  ----- */

  /**
   *  @brief  Release and close all loggers.
   *
   *  It is a wrapper function for `spdlog::drop_all` method.
   */
    inline void
  drop_all_loggers( )
  {
    spdlog::drop_all();
  }  /* -----  end of function drop_all_loggers  ----- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_LOGGER_HPP__ --- */
