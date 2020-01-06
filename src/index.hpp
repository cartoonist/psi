/**
 *    @file  index.hpp
 *   @brief  Index header file.
 *
 *  This header file contains type definitions, abstract data types, and helper
 *  interface functions (mostly implemented on top of SeqAn index library; i.e.
 *  `seqan/index.h`) for working with string indexes.
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

#ifndef  PSI_INDEX_HPP__
#define  PSI_INDEX_HPP__

#include <seqan/index.h>

#include "sequence.hpp"
#include "fmindex.hpp"


namespace grem {
  /* Typedefs  ------------------------------------------------------------------- */

  // Configured FMIndex template specialization tags.
  // :TODO:Mon Oct 23 02:11:\@cartoonist: Change all `typedef`s to `using`.
  using TFMIndexConfig = seqan::FastFMIndexConfig< void, uint64_t, 2, 1 >;
  using CFMIndex = seqan::FMIndex< void, TFMIndexConfig >;
  using CBiFMIndex = seqan::BidirectionalIndex< grem::CFMIndex >;

  template< >
    class is_fmindex< CFMIndex > : public std::true_type {
    };

  /* END OF Typedefs  ------------------------------------------------------------ */

  template< typename TText >
      inline void
    create_index( seqan::Index< TText, seqan::IndexEsa<> >& index )
    {
      indexRequire( index, seqan::FibreSA() );
      indexRequire( index, seqan::FibreLcp() );
      indexRequire( index, seqan::FibreChildtab() );
    }

  template< typename TIndex >
      inline void
    _create_fm_index( TIndex& index )
    {
      indexRequire( index, seqan::FibreSALF() );
    }

  template< typename TText >
      inline void
    create_index( seqan::Index< TText, grem::CFMIndex >& index )
    {
      _create_fm_index( index );
    }

  template< typename TText >
      inline void
    create_index( seqan::Index< TText, grem::CBiFMIndex >& index )
    {
      _create_fm_index( index );
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    create_index( seqan::Index< TText, FMIndex< TWT, TDens, TInvDens > >& index )
    {
      _create_fm_index( index );
    }

  template< typename TText, typename TIndexSpec >
      inline bool
    open( seqan::Index< TText, TIndexSpec >& index, const std::string& file_name )
    {
      return seqan::open( index, file_name.c_str() );
    }

  template< typename TText, typename TIndexSpec >
      inline bool
    save( seqan::Index< TText, TIndexSpec >& index, const std::string& file_name )
    {
      return seqan::save( index, file_name.c_str() );
    }
}  /* -----  end of namespace grem  ----- */

namespace seqan {
  /**
   *  @brief  Saving memory by overriding SAValue type.
   *
   *  @note This change limit the length of the reads to 2^32 (=4,294,967,295) in 32-bit
   *        systems and to 2^64 (=18,446,744,073,709,551,615) in 64-bit ones.
   */
  template< typename TSpec >
    struct SAValue< grem::Dna5QStringSet< TSpec > >
    {
      typedef Pair< long unsigned int, long unsigned int, Tag<Pack_> > Type;
    };
}  /* -----  end of namespace seqan  ----- */

#endif  /* --- #ifndef PSI_INDEX_HPP__ --- */
