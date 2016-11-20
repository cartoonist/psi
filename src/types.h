/*
 * =====================================================================================
 *
 * Filename: types.h
 *
 * Created: Fri Nov 11, 2016  09:40
 * Last modified: Tue Nov 15, 2016  15:35
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
}

#endif  // TYPES_H__
