/**
 *    @file  types.h
 *   @brief  Type definitions header file.
 *
 *  All global type definitions (`typedef`s) go here.
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
#include <seqan/index.h>

namespace grem
{
  typedef seqan::StringSet< seqan::CharString >                     CharStringSet;
  typedef seqan::Dna5QString                                        DnaSeq;
  typedef seqan::StringSet< DnaSeq >                                DnaSeqSet;
  // TODO: Move to a class represents a string set and its index.
  template< typename TIndexSpec >
    using DnaSeqSetIndex = seqan::Index< DnaSeqSet, TIndexSpec >;
  template< typename TIndex, typename TSpec >
    using TIndexIterator = typename seqan::Iterator< TIndex, TSpec >::Type;

  typedef struct
  {
    CharStringSet ids;
    DnaSeqSet     seqs;
    CharStringSet quals;
  } ReadsChunk;
  // END OF TODO

  // TODO: Move to a class for command-line options.
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
  // END OF TODO
}

namespace seqan
{
  /* Saving memory by overriding SAValue type:
   *
   * NOTE -- CURRENT LIMITATION:
   *   (unlimited) number of reads of length at most (2^16=65536).
   *
   */
  template<>
  struct SAValue< grem::DnaSeqSet >
  {
    typedef Pair<long unsigned int, uint16_t, Tag<Pack_> > Type;
  };
}

#endif  // TYPES_H__
