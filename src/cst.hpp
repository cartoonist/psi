/**
 *    @file  cst.h
 *   @brief  Compressed Suffix Tree wrapper on `sdsl::cst_*`
 *
 *  This is a wrapper module for `sdsl::cst_*`.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Wed Aug 01, 2018  11:24
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2018, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  CST_H__
#define  CST_H__

#include <seqan/index.h>
#include <sdsl/suffix_trees.hpp>

#include "sequence.h"
#include "utils.h"


namespace grem {
  template< class TCSA = sdsl::csa_sada<>,
    class TLCP = lcp_support_sada<>,
    class TBPS = bp_support_sada<>,
    class TRank10 = rank_support_v5<10,2>,
    class TSelect10 = select_support_mcl<10,2> >
      struct CSTSada;
}  /* -----  end of namespace grem  ----- */

namespace seqan {
  template< class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10 >
      void
    indexRequire(
        Index< grem::DiskString, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >& index,
        FibreSA );

  template< class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10 >
      void
    indexRequire(
        Index< grem::MemString, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >& index,
        FibreSA );

  template< class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10, typename TFibre >
      void
    indexRequire(
        Index< StringSet< grem::DiskString >, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >& index,
        FibreSA );

  template< class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10, typename TFibre >
      void
    indexRequire(
        Index< StringSet< grem::MemString >, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >& index,
        FibreSA );

  template< typename TText, class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10 >
      typename Index< TText, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >::text_type&
    getFibre( Index< TText, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >& index, FibreText );

  template< typename TText, class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10 >
      typename Index< TText, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > >::text_type const&
    getFibre( Index< TText, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > > const& index, FibreText );

  template< typename TText, class TCSA, class TLCP, class TBPS, class TRank10, class TSelect10 >
    class Index< TText, grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > > {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TText text_type;                  /**< @brief input text type */
        typedef typename text_type::pos_type pos_type;
        typedef grem::CSTSada< TCSA, TLCP, TBPS, TRank10, TSelect10 > spec_type;
        typedef sdsl::cst_sada< TCSA, TLCP, TBPS, TRank10, TSelect10 > value_type;
        typedef typename value_type::string_type string_type;  /**< @brief output string type */
        typedef typename value_type::size_type savalue_type;
        typedef typename value_type::index_category index_category;
        typedef typename value_type::char_type char_type;
        typedef typename value_type::comp_char_type comp_char_type;
        /* ====================  LIFECYCLE     ======================================= */
        Index ( )
          : text_p( nullptr ), owner( true ) { }
        Index( text_type* t )
          : text_p( t ), owner( false ) { }
        Index( text_type& t )
          : Index( &t ) { }

        Index( Index const& ) = delete;
        Index( Index&& ) = default;
        Index& operator=( Index const& ) = delete;
        Index& operator=( Index&& ) = default;
        ~Index( ) noexcept
        {
          this->clear();
        }
        /* ====================  ACCESSORS     ======================================= */
          inline bool
        owns_text( ) const
        {
          return this->owner;
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        value_type cst;
        text_type* text_p;
        bool owner;
        /* ====================  METHODS              ================================ */
        // :TODO:Wed Aug 01 15:09:\@cartoonist: FIXME: a Holder class should be
        //     responsible for the text fibre not client class.
        //     === PROBLEM ===
        //     The index has a pointer to the text fibre which it might not own it. It
        //     is valid until the text object lives. When a class aggregates an index
        //     (particularly, CSTSada) and its text fibre when an object of that class
        //     is moved to another one. The text fibre of new object's index is no
        //     longer valid since it points to the text fibre copy of the moved object
        //     not the new one.
        //     This is the problem for all classes that aggregates the index and its
        //     text fibre. The text fibre should be owned by a Holder class not the
        //     client class (e.g. PathSet, or PathIndex here).
        //     === SOLUTION ===
        //     This function is a work-around for this issue. The true solution would be
        //     wrapping the text fibre in a Holder class which is responsible for the
        //     memory management of the text fibre not the class aggregating the index
        //     and its text. In other words, the text fibre object will be handed to
        //     objects on the fly and finaly the last holder would free the memory.
        //     This Holder class can be developed from scratch or use smart pointers in
        //     C++ Standard Library.
          inline void
        set_text_fibre( text_type* ext_p, bool update=true )
        {
          if ( update ) this->clear_fibres();
          this->text_p = ext_p;
          this->owner = false;
        }
        /* ====================  INTERFACE FUNCTIONS  ================================ */
          friend void
        indexRequire< TCSA, TLCP, TBPS, TRank10, TSelect10 >( Index& index, FibreSA );
          friend text_type&
        getFibre< TText, TCSA, TLCP, TBPS, TRank10, TSelect10 >( Index&, FibreText );
          friend text_type const&
        getFibre< TText, TCSA, TLCP, TBPS, TRank10, TSelect10 >( Index const&, FibreText );
        /* ====================  FRIENDSHIPS   ======================================= */
        friend class Finder< Index >;
        friend class Iter< Index, seqan::TopDown<> >;
        friend class Iter< Index, seqan::TopDown< seqan::ParentLinks<> > >;
    };
}  /* -----  end of namespace seqan  ----- */

#endif  /* ----- #ifndef CST_H__  ----- */
