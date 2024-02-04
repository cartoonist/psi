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


namespace psi {
  /* Typedefs  ------------------------------------------------------------------- */

  // Configured FMIndex template specialization tags.
  // :TODO:Mon Oct 23 02:11:\@cartoonist: Change all `typedef`s to `using`.
  using TFMIndexConfig = seqan2::FastFMIndexConfig< void, uint64_t, 2, 1 >;
  using CFMIndex = seqan2::FMIndex< void, TFMIndexConfig >;
  using CBiFMIndex = seqan2::BidirectionalIndex< psi::CFMIndex >;

  template< >
    class is_fmindex< CFMIndex > : public std::true_type {
    };

  /* END OF Typedefs  ------------------------------------------------------------ */

  template< typename TText >
      inline void
    create_index( seqan2::Index< TText, seqan2::IndexEsa<> >& index )
    {
      indexRequire( index, seqan2::FibreSA() );
      indexRequire( index, seqan2::FibreLcp() );
      indexRequire( index, seqan2::FibreChildtab() );
    }

  template< typename TIndex >
      inline void
    _create_fm_index( TIndex& index )
    {
      indexRequire( index, seqan2::FibreSALF() );
    }

  template< typename TText >
      inline void
    create_index( seqan2::Index< TText, psi::CFMIndex >& index )
    {
      _create_fm_index( index );
    }

  template< typename TText >
      inline void
    create_index( seqan2::Index< TText, psi::CBiFMIndex >& index )
    {
      _create_fm_index( index );
    }

  template< typename TText, class TWT, uint32_t TDens, uint32_t TInvDens >
      inline void
    create_index( seqan2::Index< TText, FMIndex< TWT, TDens, TInvDens > >& index )
    {
      _create_fm_index( index );
    }

  template< typename TText, typename TIndexSpec >
      inline bool
    open( seqan2::Index< TText, TIndexSpec >& index, const std::string& file_name )
    {
      return seqan2::open( index, file_name.c_str() );
    }

  template< typename TText, typename TIndexSpec >
      inline bool
    save( seqan2::Index< TText, TIndexSpec >& index, const std::string& file_name )
    {
      return seqan2::save( index, file_name.c_str() );
    }
}  /* --- end of namespace psi --- */

namespace seqan2 {
  /**
   *  @brief  Saving memory by overriding SAValue type.
   *
   *  @note This change limit the length of the reads to 2^32 (=4,294,967,295) in 32-bit
   *        systems and to 2^64 (=18,446,744,073,709,551,615) in 64-bit ones.
   */
  template< typename TSpec >
    struct SAValue< psi::Dna5QStringSet< TSpec > >
    {
      typedef Pair< long unsigned int, long unsigned int, Tag<Pack_> > Type;
    };
}  /* -----  end of namespace seqan2  ----- */

#endif  /* --- #ifndef PSI_INDEX_HPP__ --- */
