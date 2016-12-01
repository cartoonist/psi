/*
 * =====================================================================================
 *
 * Filename: types.h
 *
 * Created: Fri Nov 11, 2016  09:40
 * Last modified: Thu Dec 01, 2016  16:08
 *
 * Description: Types header file.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef TYPES_H__
#define TYPES_H__

#include <cstdint>

#include <seqan/sequence.h>
#include <seqan/index.h>

namespace grem
{
  typedef int64_t id_t;

  typedef seqan::StringSet< seqan::CharString >                     CharStringSet;
  typedef seqan::Dna5QString                                        DnaSeq;
  typedef seqan::StringSet< DnaSeq >                                DnaSeqSet;
  typedef seqan::Index< DnaSeqSet, seqan::IndexWotd<> >             DnaSeqSetIndex;
  typedef seqan::Iterator< DnaSeqSetIndex, seqan::TopDown<> >::Type DnaSSIndexIter;

  typedef struct
  {
    CharStringSet ids;
    DnaSeqSet     seqs;
    CharStringSet quals;
  } ReadsChunk;

  typedef struct
  {
    unsigned int seed_len;
    unsigned int chunk_size;
    seqan::CharString rf_path;
    seqan::CharString fq_path;
  } GremOptions;
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
