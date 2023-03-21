/**
 *    @file  basic_types.hpp
 *   @brief  Basic type definitions.
 *
 *  This header file defines some basic types.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Apr 03, 2023  18:48
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2019, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_BASIC_TYPES_HPP__
#define PSI_BASIC_TYPES_HPP__

namespace psi {
  template< typename TId = typename gum::GraphBaseTrait< gum::Dynamic >::id_type,
            typename TOffset = typename gum::GraphBaseTrait< gum::Dynamic >::offset_type >
  class Position {
    public:
      /* === TYPE MEMBERS === */
      using id_type = TId;
      using offset_type = TOffset;
      /* === LIFECYCLE === */
      Position( id_type id=0, offset_type offset=0 )
        : m_id( id ), m_offset( offset )
      { }

      Position( const Position& ) = default;
      Position( Position&& ) = default;
      ~Position( ) = default;
      /* === OPERATORS === */
      Position& operator=( const Position& ) = default;
      Position& operator=( Position&& ) = default;
      /* === ACCESSORS === */
      inline id_type
      node_id( ) const
      {
        return this->m_id;
      }

      inline offset_type
      offset( ) const
      {
        return this->m_offset;
      }
      /* === MUTATORS === */
      inline void
      set_node_id( id_type e_id )
      {
        this->m_id = e_id;
      }

      inline void
      set_offset( offset_type e_offset )
      {
        this->m_offset = e_offset;
      }
    private:
      /* === DATA MEMBERS === */
      TId m_id;
      TOffset m_offset;
  };
}  /* --- end of namespace psi --- */

#endif // PSI_BASIC_TYPES_HPP__
