/**
 *    @file  index.h
 *   @brief  Index header file.
 *
 *  This header file contains type definitions, abstract data types, and helper meta-
 *  functions (mostly implemented on top of SeqAn index library; i.e. `seqan/index.h`)
 *  for working with string indexes.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Tue Mar 07, 2017  20:20
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  INDEX_H__
#define  INDEX_H__

#include <seqan/index.h>

#include "sequence.h"

namespace grem {
  /* Typedefs  ------------------------------------------------------------------- */
  template < typename TIndexSpec >
    using Dna5QStringSetIndex = seqan::Index< Dna5QStringSet, TIndexSpec >;

  typedef seqan::FastFMIndexConfig<void, uint64_t, 2, 1> TFMIndexConfig;

  template < typename TText >
    using TBiFMIndex = seqan::Index < TText,
          seqan::BidirectionalIndex < seqan::FMIndex < void, TFMIndexConfig > > >;

  template < typename TText >
    using TFMIndex = seqan::Index < TText, seqan::FMIndex < void, TFMIndexConfig > >;
  /* END OF Typedefs  ------------------------------------------------------------ */

  template < typename TText >
    void
  create_index ( seqan::Index < TText, seqan::IndexEsa<> > & index )
  {
    indexRequire ( index, seqan::EsaSA() );
    indexRequire ( index, seqan::EsaLcp() );
    indexRequire ( index, seqan::EsaChildtab() );
    indexRequire ( index, seqan::EsaBwt() );
  }

}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Saving memory by overriding SAValue type.
   *
   *  @note This change limit the length of the reads to 2^16 (=65536).
   */
  template < >
  struct SAValue< StringSet < Dna5QString > >
  {
    typedef Pair<long unsigned int, long unsigned int, Tag<Pack_> > Type;
  };
}  /* -----  end of namespace seqan  ----- */

#endif  /* ----- #ifndef INDEX_H__  ----- */