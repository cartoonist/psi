/**
 *    @file  logger.h
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

#ifndef LOGGER_H__
#define LOGGER_H__

#include "base.h"
#include "options.h"

#include "spdlog/spdlog.h"


namespace grem {
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
  config_logger ( const Options& options )
  {
    // Set to asynchronous mode.
    spdlog::set_async_mode(8192);
    std::vector< spdlog::sink_ptr > sinks;

    if ( !options.nolog && !options.quiet ) {
      if ( options.nocolor ) {
        auto stdout_sink = spdlog::sinks::stdout_sink_mt::instance();
        sinks.push_back( stdout_sink );
      }
      else {
        auto color_sink = std::make_shared< spdlog::sinks::ansicolor_stdout_sink_mt >();
        sinks.push_back( color_sink );
      }
    }

    if ( !options.nolog && !options.nologfile ) {
      auto simple
        = std::make_shared< spdlog::sinks::simple_file_sink_mt >( options.log_path );
      sinks.push_back( simple );
    }

    auto main_logger
      = std::make_shared< spdlog::logger >( "main", begin(sinks), end(sinks) );
    if ( options.verbose ) {
      main_logger->set_level( spdlog::level::info );
    }
    else {
      main_logger->set_level( spdlog::level::warn );
    }

    spdlog::register_logger( main_logger );
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
}  /* -----  end of namespace grem  ----- */
#endif  // LOGGER_H__
