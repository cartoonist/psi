/**
 *    @file  options.h
 *   @brief  Options class definition.
 *
 *  It contains data structures storing program options.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Fri Nov 11, 2016  09:40
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef TYPES_H__
#define TYPES_H__

#include <seqan/sequence.h>

namespace grem
{
  // :TODO:Tue Mar 07 20:34:\@cartoonist: Option class.
  enum class IndexType {
    Sa = 1,               /**< @brief Suffix array index. */
    Esa,                  /**< @brief Enhanced suffix array index. */
    Wotd,                 /**< @brief Lazy suffix tree (write only, top down) index. */
    Dfi,                  /**< @brief Deferred frequency index. */
    QGram,                /**< @brief An index based on an array of sorted q-grams. */
    FM                    /**< @brief FM index. */
  };

  typedef struct
  {
    unsigned int seed_len;
    unsigned int chunk_size;
    unsigned int start_every;
    IndexType index;
    seqan::CharString rf_path;
    seqan::CharString fq_path;
    seqan::CharString log_path;
    bool nologfile;
    bool nolog;
    bool quiet;
    bool nocolor;
  } GremOptions;
}

#endif  // TYPES_H__
