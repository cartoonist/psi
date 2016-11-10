/*
 * =====================================================================================
 *
 * Filename: types.h
 *
 * Created: Fri Nov 11, 2016  09:40
 * Last modified: Mon Nov 14, 2016  00:38
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

namespace grem
{
  typedef int64_t id_t;

  typedef struct
  {
    seqan::StringSet< seqan::CharString > ids;
    seqan::StringSet< seqan::Dna5String > seqs;
    seqan::StringSet< seqan::CharString > quals;
  } ReadsChunk;
}

#endif  // TYPES_H__
