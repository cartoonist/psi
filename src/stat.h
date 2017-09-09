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
}  /* -----  end of namespace grem  ----- */

#endif  // end of STAT_H__
