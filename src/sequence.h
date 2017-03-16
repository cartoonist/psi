/**
 *    @file  sequence.h
 *   @brief  Sequence abstract data types.
 *
 *  Sequence header file contains type definitions, abstract data types, and helper meta-
 *  functions (mostly implemented at the top of SeqAn sequence library; i.e.
 *  `seqan/sequence.h`) to work with sequences and sequence sets.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 07, 2017  12:50
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  SEQUENCE_H__
#define  SEQUENCE_H__

#include <seqan/sequence.h>

namespace grem {
  /* Typedefs  ------------------------------------------------------------------- */
  typedef seqan::StringSet< seqan::CharString >    CharStringSet;
  typedef seqan::StringSet< seqan::Dna5QString >   Dna5QStringSet;
  typedef seqan::Position < Dna5QStringSet >::Type Dna5QStringSetPosition;
  /* END OF Typedefs  ------------------------------------------------------------ */

  /* Data structures  ------------------------------------------------------------ */

  /**
   *  @brief  String set with an ID associated with each string.
   *
   *  It is a wrapper class on StringSet associating an ID to each string.
   */
  template < typename TStringSet >
    class NamedStringSet
    {
      public:
        /* ====================  DATA MEMBERS  ======================================= */
        TStringSet str;
        CharStringSet id;
    };  /* ----------  end of template class NamedStringSet  ---------- */

  typedef NamedStringSet < Dna5QStringSet > Dna5QRecords;

  /* END OF Data structures  ----------------------------------------------------- */

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef SEQUENCE_H__  ----- */
