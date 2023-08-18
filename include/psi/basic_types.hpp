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
  template< typename TId, typename TOffset >
  class PositionBase {
    public:
      /* === TYPE MEMBERS === */
      using id_type = TId;
      using offset_type = TOffset;
      /* === LIFECYCLE === */
      PositionBase( id_type id=0, offset_type offset=0 )
        : m_id( id ), m_offset( offset )
      { }

      PositionBase( const PositionBase& ) = default;
      PositionBase( PositionBase&& ) = default;
      ~PositionBase( ) = default;
      /* === OPERATORS === */
      PositionBase& operator=( const PositionBase& ) = default;
      PositionBase& operator=( PositionBase&& ) = default;
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

  /* https://stackoverflow.com/a/67642371/357257 */
  template< typename I >
  class RangeIterator
  {
    private:
      /* === DATA MEMBERS === */
      I i;
    public:
      typedef I difference_type;
      typedef I value_type;
      typedef I pointer;
      typedef I reference;
      typedef std::random_access_iterator_tag iterator_category;

      RangeIterator( I i ) : i( i ) { }

      bool operator==( RangeIterator<I>& other ) { return i == other.i; }
      I operator-( RangeIterator<I>& other ) { return i - other.i; }
      I operator++() { return i++; }
      I operator*() { return i; }
  };
}  /* --- end of namespace psi --- */

#endif // PSI_BASIC_TYPES_HPP__
