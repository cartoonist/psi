/**
 *    @file  hbitvector.hpp
 *   @brief  HBitVector class definition header file.
 *
 *  This header file contains `HBitVector` class definition and related
 *  interface functions and type definitions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Tue Jun 20, 2023  22:23
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2019, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_HBITVECTOR_HPP_
#define PSI_HBITVECTOR_HPP_

#include <climits>
#include <cstdint>

#include <gum/basic_types.hpp>
#include <sdsl/bits.hpp>
#include <Kokkos_Core.hpp>


namespace psi {
  /**
   *  @brief  Hierarchical (two-level) bit vector
   *
   *  @tparam  TL1_Size Size of the first level bit vector
   *
   *  The arrangement of two bit arrays for each level is shown below:
   *
   *  g = global index, b = L1 begin position, r = relative index, l = local index
   *  g: b b+1  ...  b+|L1|-1   b+|L1| ... n 0 1 ...  b-1
   *     | |               |    |          | | |       |
   *   [     L1 region      ] [        L2 region         ]
   *     | |               |    | | |                  |
   *  l: 0 1    ...    |L1|-1   0 1 2      ...     |L2|-1
   *
   *  r = ( ( g-b ) % n )
   *  l = r < |L1| ? r : r - |L1|
   */
  template< int TL1_Size=2048, /*bits*/
            typename TSize=uint32_t,
            typename TBitset=uint64_t,
            typename TDevice=Kokkos::DefaultExecutionSpace >
  class HBitVector {
    public:
      /* === MEMBER TYPES === */
      using device_type = TDevice;
      using execution_space = typename device_type::execution_space;
      using scratch_space = typename execution_space::scratch_memory_space;
      using member_type = typename Kokkos::TeamPolicy< execution_space >::member_type;
      using size_type = TSize;
      using bitset_type = TBitset;
      using view_type = bitset_type*;
      using value_type = bool;
      /* === MEMBER CONSTANTS === */
      static constexpr const size_type BITSET_WIDTH            = gum::widthof< bitset_type >::value;  // 64 bits (if uint64_t)
      static constexpr const unsigned short int BITSET_SHIFT   = sdsl::bits::hi( BITSET_WIDTH );      // 6       (if uint64_t)
      static constexpr const bitset_type BITSET_MASK           = BITSET_WIDTH - 1u;        // 0x000000000000003f (if uint64_t)
      static constexpr const bitset_type BITSET_N_MASK         = ~BITSET_MASK;             // 0xffffffffffffffc0 (if uint64_t)
      static constexpr const bitset_type BITSET_ALL_NIL        = 0;                        // 0x0000000000000000 (if uint64_t)
      static constexpr const bitset_type BITSET_ALL_SET        = ~BITSET_ALL_NIL;          // 0x1111111111111111 (if uint64_t)
      static constexpr const unsigned short int L1_SIZE        = TL1_Size;
      static constexpr const unsigned short int L1_SIZE_BYTES  = TL1_Size >> 3;
      static constexpr const unsigned short int L1_NUM_BITSETS = TL1_Size >> BITSET_SHIFT;
      static constexpr std::size_t value_alignment             = std::max( sizeof( bitset_type ), alignof( bitset_type ) );
      static constexpr std::size_t space_alignment             = std::max( value_alignment, static_cast< size_t >( scratch_space::ALIGN ) );
      /* === STATIC ASSERTS === */
      // Accepting 64-bit L1 size requires spending extra time checking for corner cases
      static_assert( ( TL1_Size >= ( BITSET_WIDTH << 1 ) ), "L1 size should be at least twice larger than bitset width" );
      static_assert( ( sdsl::bits::cnt( TL1_Size ) == 1 ), "L1 size should be a power of 2" );
      static_assert( ( sdsl::bits::cnt( BITSET_WIDTH ) == 1 ), "Bitset width should be a power of 2" );
      static_assert( ( TL1_Size <= std::numeric_limits< size_type >::max() ), "L1 size cannot fit in size type" );
      /* === DATA MEMBERS === */
      size_type m_size;    //!< Size of the bit vector
      size_type l1_begin;  //!< Index of the first bit reside in L1 (inclusive)
      view_type l1_data;   //!< First level bit vector view
      view_type l2_data;   //!< Second level bit vector view
      /* === LIFECYCLE === */
      HBitVector( size_type n, size_type centre, const member_type& tm )
        : m_size( n ), l1_begin( 0 ), l1_data( ), l2_data( )
      {
        assert( centre < n );

        auto bitset = ( centre >> BITSET_SHIFT );
        // range is [bitset-(L1_NUM_BITSETS/2)+1...bitset+(L1_NUM_BITSETS/2)] (inclusive)
        auto nof_bitsets = HBitVector::num_bitsets( n );
        if ( L1_NUM_BITSETS < nof_bitsets ) {
          unsigned short int l_pad = ( L1_NUM_BITSETS >> 1 ) - 1;
          auto r_fit = nof_bitsets - L1_NUM_BITSETS;
          auto begin_bitset = ( ( bitset > l_pad ) ? ( bitset - l_pad ) : 0 );
          // for values of `centre` being closer to the end, l1 covers the last `TL1_Size` bits
          begin_bitset = ( ( r_fit < begin_bitset ) ? r_fit : begin_bitset );
          this->l1_begin = begin_bitset << BITSET_SHIFT;
        }
        else this->l1_begin = 0;

        auto l1size = HBitVector::l1_scratch_size();
        this->l1_data = ( view_type )
          ( tm.team_scratch( 0 ).get_shmem_aligned( l1size, space_alignment ) );
        auto l2size = HBitVector::l2_scratch_size( n );
        if ( l2size != 0 ) {
          this->l2_data = ( view_type )
            ( tm.team_scratch( 1 ).get_shmem_aligned( l2size, space_alignment ) );
        }
      }
      /* === STATIC MEMBERS === */
      static constexpr inline size_type
      l1_num_bitsets( )
      {
        return L1_NUM_BITSETS;
      }

      static constexpr inline size_type
      l1_size( )
      {
        return TL1_Size;
      }

      static constexpr inline size_type
      l1_scratch_size( )
      {
        return L1_SIZE_BYTES;
      }

      static constexpr inline size_type
      num_bitsets( size_type n )
      {
        return ( n >> BITSET_SHIFT ) +     // returns `n >> BITSET_SHIFT` if `n` is divisible by `BITSET_WIDTH`
          ( (bool)( n % BITSET_WIDTH ) );  // otherwise, returns `( n >> BITSET_SHIFT ) + 1`
      }

      static inline size_type
      l2_num_bitsets( size_type n )
      {
        auto nbitsets = HBitVector::num_bitsets( n );
        return ( nbitsets  > HBitVector::l1_num_bitsets() ) ? nbitsets - HBitVector::l1_num_bitsets() : 0;
      }

      static inline size_type
      l2_scratch_size( size_type n )
      {
        return HBitVector::l2_num_bitsets( n ) * sizeof( bitset_type );
      }

      static inline size_type
      l2_size( size_type n )
      {
        return HBitVector::l2_scratch_size( n ) * CHAR_BIT;
      }

      static inline size_type
      capacity( size_type n )
      {
        return HBitVector::l2_scratch_size( n ) +
               HBitVector::l1_scratch_size();
      }

      template< typename TPolicy >
      static inline TPolicy
      set_scratch_size( TPolicy& policy, size_type n )
      {
        auto l1size = HBitVector::l1_scratch_size();
        auto l2size = HBitVector::l2_scratch_size( n );
        auto new_policy = policy.set_scratch_size( 0, Kokkos::PerTeam( l1size ) );
        if ( l2size == 0 ) return new_policy;
        return new_policy.set_scratch_size( 1, Kokkos::PerTeam( l2size ) );
      }
      /* === METHODS === */
      inline size_type
      size( ) const
      {
        return this->m_size;
      }

      inline size_type
      capacity( ) const
      {
        return HBitVector::capacity( this->m_size );
      }

      inline size_type
      num_bitsets( ) const
      {
        return HBitVector::num_bitsets( this->m_size );
      }

      inline size_type
      l2_size() const
      {
        return HBitVector::l2_size( this->m_size );
      }

      inline size_type
      l2_num_bitsets( ) const
      {
        return HBitVector::l2_num_bitsets( this->m_size );
      }

      inline size_type
      l2_scratch_size( ) const
      {
        return HBitVector::l2_scratch_size( this->m_size );
      }

      inline size_type
      relative_idx( size_type i )
      {
        return ( static_cast< std::ptrdiff_t >( i ) - static_cast< std::ptrdiff_t >( this->l1_begin ) ) % this->m_size;
      }

      inline size_type
      relative_bitset( size_type bitset_idx ) const
      {
        auto start_idx = bitset_idx << BITSET_SHIFT;
        auto ridx = this->relative_idx( start_idx );
        return ridx >> BITSET_SHIFT;
      }

      inline void
      clear_l1( const member_type& tm )
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, 0, l1_num_bitsets() ),
            [=]( const uint64_t j ) { this->l1_data[ j ] = 0; } );
      }

      /**
       *   @brief Clear the range of bitsets [start, end) which occurs on L2
       *
       *   @param tm  Team Policy member
       *   @param start Absolute bitset start index in the vector
       *   @param end Absolute bitset end index in the vector
       *
       *   NOTE: `start` and `end` are absolute bitset indexes (not relative).
       *   NOTE: The caller should make sure that both [start, end) does not
       *         span over L1.
       */
      inline void
      clear_l2_by_bidx( const member_type& tm, size_type start, size_type end )
      {
        assert( end != 0 );

        auto rstart = this->relative_bitset( start ) - l1_num_bitsets();
        auto rend = this->relative_bitset( end - 1 ) + 1 - l1_num_bitsets();
        Kokkos::parallel_for( Kokkos::TeamVectorRange( tm, rstart, rend ),
                              [=]( const uint64_t j ) {
                                this->l2_data[ j ] = 0;
                              } );
      }

      /**
       *   @brief Clear the range of bitsets [start, end) which occurs on L2
       *
       *   @param tm  Team Policy member
       *   @param start Absolute start index in the vector
       *   @param end Absolute end index in the vector
       *
       *   NOTE: `start` and `end` are absolute bit indexes (not relative).
       *   NOTE: The caller should make sure that both [start, end) does not
       *         span over L1.
       */
      inline void
      clear_l2_by_idx( const member_type& tm, size_type start, size_type end )
      {
        assert( end != 0 );

        auto rstart = this->relative_idx( start ) - l1_size();
        auto rend = this->relative_idx( end - 1 ) + 1 - l1_size();
        auto lbs = ( rstart >> BITSET_SHIFT );
        auto rbs = ( rend >> BITSET_SHIFT );
        Kokkos::parallel_for( Kokkos::TeamVectorRange( tm, lbs, rbs ),
                              [=]( const uint64_t j ) {
                                this->l2_data[ j ] = 0;
                              } );
      }

      /**
       *   @brief Clear all bitsets in L2
       */
      inline void
      clear_l2( const member_type& tm )
      {
        Kokkos::parallel_for( Kokkos::TeamVectorRange( tm, this->l2_num_bitsets() ),
                              [=]( const uint64_t j ) {
                                this->l2_data[ j ] = 0;
                              } );
      }

      inline void
      set( size_type i )
      {
        assert( i < this->m_size );

        auto ridx = this->relative_idx( i );
        auto bitset = ( ridx >> BITSET_SHIFT );
        auto mask = 1ULL << ( ridx & BITSET_MASK );

        if ( ridx < this->l1_size() ) {  // most probable
          this->l1_data[ bitset ] |= mask;
        }
        else {
          bitset -= this->l1_num_bitsets();
          this->l2_data[ bitset ] |= mask;
        }
      }

      inline void
      set( size_type s, size_type f )
      {
        assert( s <= f );

        if ( s == f ) return this->set( f );

        auto s_ridx = this->relative_idx( s );
        auto f_ridx = this->relative_idx( f );

        auto setbits =
          []( auto data_ptr, auto s_lidx, auto f_lidx ) {
            auto s_bitset  = ( s_lidx >> BITSET_SHIFT );
            auto f_bitset = ( f_lidx >> BITSET_SHIFT );
            auto s_offset = ( s_lidx & BITSET_MASK );
            auto f_offset = ( f_lidx & BITSET_MASK );

            if ( s_bitset != f_bitset ) {
              auto mask = ( BITSET_ALL_SET << s_offset );
              data_ptr[ s_bitset ] |= mask;
              for ( auto i = s_bitset + 1; i < f_bitset; ++i ) data_ptr[ i ] |= BITSET_ALL_SET;
              mask = ( BITSET_ALL_SET >> ( BITSET_MASK - f_offset ) );
              data_ptr[ f_bitset ] |= mask;
            }
            else {
              auto mask = ( ( 1ULL << ( f_offset - s_offset + 1 ) ) - 1 ) << s_offset;
              data_ptr[ s_bitset ] |= mask;
            }
          };

        if ( s_ridx < this->l1_size() && f_ridx < this->l1_size() ) {  // the range is in L1
          setbits( this->l1_data, s_ridx, f_ridx );
        }
        else if ( s_ridx < this->l1_size() && this->l1_size() <= f_ridx ) {
          auto f_lidx = f_ridx - this->l1_size();
          setbits( this->l1_data, s_ridx, this->l1_size() - 1 );
          setbits( this->l2_data, 0, f_lidx );
        }
        else if ( this->l1_size() <= s_ridx && s_ridx <= f_ridx ) {
          auto s_lidx = s_ridx - this->l1_size();
          auto f_lidx = f_ridx - this->l1_size();
          setbits( this->l2_data, s_lidx, f_lidx );
        }
        else if ( this->l1_size() <= s_ridx && f_ridx < this->l1_size() ) {
          auto s_lidx = s_ridx - this->l1_size();
          setbits( this->l1_data, 0, f_ridx );
          setbits( this->l2_data, s_lidx, this->l2_size() - 1 );
        }
        else {
          auto s_lidx = s_ridx - this->l1_size();
          auto f_lidx = f_ridx - this->l1_size();
          setbits( this->l1_data, 0, this->l1_size() - 1 );
          setbits( this->l2_data, s_lidx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, f_lidx );
        }
      }

      inline void
      set( const member_type& tm, size_type s, size_type f )
      {
        assert( s <= f );

        if ( s == f ) return this->set( f );

        auto s_ridx = this->relative_idx( s );
        auto f_ridx = this->relative_idx( f );

        auto setbits =
          [&tm]( auto data_ptr, auto s_lidx, auto f_lidx ) {
            auto s_bitset  = ( s_lidx >> BITSET_SHIFT );
            auto f_bitset = ( f_lidx >> BITSET_SHIFT );
            auto s_offset = ( s_lidx & BITSET_MASK );
            auto f_offset = ( f_lidx & BITSET_MASK );

            if ( s_bitset != f_bitset ) {
              auto mask = ( BITSET_ALL_SET << s_offset );
              data_ptr[ s_bitset ] |= mask;

              Kokkos::parallel_for( Kokkos::ThreadVectorRange( tm, s_bitset+1, f_bitset ),
                                    [=] ( const uint64_t k ) {
                                      data_ptr[ k ] |= BITSET_ALL_SET;
                                    } );

              mask = ( BITSET_ALL_SET >> ( BITSET_MASK - f_offset ) );
              data_ptr[ f_bitset ] |= mask;
            }
            else {
              auto mask = ( ( 1ULL << ( f_offset - s_offset + 1 ) ) - 1 ) << s_offset;
              data_ptr[ s_bitset ] |= mask;
            }
          };

        if ( s_ridx < this->l1_size() && f_ridx < this->l1_size() ) {  // the range is in L1
          setbits( this->l1_data, s_ridx, f_ridx );
        }
        else if ( s_ridx < this->l1_size() && this->l1_size() <= f_ridx ) {
          auto f_lidx = f_ridx - this->l1_size();
          setbits( this->l1_data, s_ridx, this->l1_size() - 1 );
          setbits( this->l2_data, 0, f_lidx );
        }
        else if ( this->l1_size() <= s_ridx && s_ridx <= f_ridx ) {
          auto s_lidx = s_ridx - this->l1_size();
          auto f_lidx = f_ridx - this->l1_size();
          setbits( this->l2_data, s_lidx, f_lidx );
        }
        else if ( this->l1_size() <= s_ridx && f_ridx < this->l1_size() ) {
          auto s_lidx = s_ridx - this->l1_size();
          setbits( this->l1_data, 0, f_ridx );
          setbits( this->l2_data, s_lidx, this->l2_size() - 1 );
        }
        else {
          auto s_lidx = s_ridx - this->l1_size();
          auto f_lidx = f_ridx - this->l1_size();
          setbits( this->l1_data, 0, this->l1_size() - 1 );
          setbits( this->l2_data, s_lidx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, f_lidx );
        }
      }
  };
}  /* --- end of namespace psi --- */


#endif // PSI_HBITVECTOR_HPP_
