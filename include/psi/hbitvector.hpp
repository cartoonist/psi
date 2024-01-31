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
#include <limits>

#include <gum/basic_types.hpp>
#include <sdsl/bits.hpp>
#include <Kokkos_Core.hpp>

#include "range_sparse_base.hpp"


namespace psi {
  template< typename TDevice >
  struct HBitVectorTraits {
    using size_type = uint32_t;
    using bitset_type = uint64_t;
  };

#if defined( KOKKOS_ENABLE_CUDA )
  template< >
  struct HBitVectorTraits< Kokkos::Cuda > {
    using size_type = uint32_t;
    using bitset_type = uint32_t;
  };
#endif

  /**
   *  @brief  Hierarchical (two-level) bit vector
   *
   *  @tparam  TL1_Size Size of the first level bit vector (faster)
   *
   *  The arrangement of two bit arrays for each level is shown below:
   *
   *  n = vector size, N = aligned vector size (multiple of bitset width)
   *  g = global index, b = L1 begin position, r = relative index, l = local index
   *  g: b b+1  ...  b+|L1|-1   b+|L1| ... n 0 1 ...  b-1
   *     | |               |    |          | | |       |
   *   [     L1 region      ] [        L2 region         ]
   *     | |               |    | | |                  |
   *  l: 0 1    ...    |L1|-1   0 1 2      ...     |L2|-1
   *
   *  r = ( ( N + g - b ) % N )
   *  l = r < |L1| ? r : r - |L1|
   */
  template< unsigned int TL1_Size = 2048, /*bits*/
            typename TDevice = Kokkos::DefaultExecutionSpace,
            typename TTrait = HBitVectorTraits< TDevice > >
  class HBitVector {
    public:
      /* === MEMBER TYPES === */
      using device_type = TDevice;
      using trait_type = TTrait;
      using execution_space = typename device_type::execution_space;
      using scratch_space = typename execution_space::scratch_memory_space;
      using policy_type = typename Kokkos::TeamPolicy< execution_space >;
      using member_type = typename policy_type::member_type;
      using size_type = typename trait_type::size_type;
      using bitset_type = typename trait_type::bitset_type;
      using view_type = bitset_type*;
      using value_type = bool;
      /* === MEMBER CONSTANTS === */
      static constexpr const size_type BITSET_WIDTH            = gum::widthof< bitset_type >::value;  // 64 (if uint64_t)
      static constexpr const unsigned short int BINDEX_SHIFT   = sdsl::bits::hi( BITSET_WIDTH );      // 6  (if uint64_t)
      static constexpr const size_type BOFFSET_MASK            = BITSET_WIDTH - 1u; // 0x0000001f (if size_type=uint32_t)
      static constexpr const size_type INDEX_ALIGN_MASK        = ~BOFFSET_MASK;     // 0xffffffe0 (if size_type=uint32_t)
      static constexpr const bitset_type BITSET_ALL_NIL        = 0;                   // 0x0000000000000000 (if uint64_t)
      static constexpr const bitset_type BITSET_ALL_SET        = ~BITSET_ALL_NIL;     // 0xffffffffffffffff (if uint64_t)
      static constexpr const bitset_type BITSET_ONE            = 1u;
      static constexpr const size_type L1_SIZE                 = TL1_Size;
      static constexpr const size_type L1_SIZE_BYTES           = TL1_Size / CHAR_BIT;
      static constexpr const size_type L1_NUM_BITSETS          = TL1_Size >> BINDEX_SHIFT;
      static constexpr const std::size_t value_alignment       = Kokkos::max( sizeof( bitset_type ), alignof( bitset_type ) );
      static constexpr const std::size_t space_alignment       = Kokkos::max( value_alignment, static_cast< size_t >( scratch_space::ALIGN ) );
      /* === STATIC ASSERTS === */
      // Accepting 64-bit L1 size requires spending extra time checking for corner cases
      static_assert( ( TL1_Size >= ( BITSET_WIDTH << 1 ) ), "L1 size should be at least twice larger than bitset width" );
      static_assert( ( sdsl::bits::cnt( TL1_Size ) == 1 ), "L1 size should be a power of 2" );
      static_assert( ( sdsl::bits::cnt( BITSET_WIDTH ) == 1 ), "Bitset width should be a power of 2" );
      static_assert( ( TL1_Size <= std::numeric_limits< size_type >::max() ), "L1 size cannot fit in size type" );
      /* === DATA MEMBERS === */
      //size_type m_size;         //!< Size of the bit vector
      size_type m_x_size;       //!< Allocated size in bits (always a multiple of `BITSET_WIDTH`)
      size_type m_num_bitsets;  //!< Total number of bitsets
      size_type l1_begin_bidx;  //!< Bitset index of the first bitset in L1 (inclusive)
      size_type l1_begin;       //!< Index of the first bit reside in L1 (inclusive)
      view_type l1_data;        //!< First level bit vector view
      view_type l2_data;        //!< Second level bit vector view
      /* === LIFECYCLE === */
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, size_type n, size_type centre )
        : /*m_size( n ),*/ l1_begin( 0 ), l1_data( ), l2_data( )
      {
        assert( centre < n );

        this->m_x_size = HBitVector::aligned_size( n );
        this->m_num_bitsets = HBitVector::bindex( this->m_x_size );
        auto ctr_bidx = HBitVector::bindex( centre );
        // The begin index is set such that the L1 range would be (inclusive):
        //   [ctr_bidx-(L1_NUM_BITSETS/2)+1...ctr_bidx+(L1_NUM_BITSETS/2)]
        if ( L1_NUM_BITSETS < this->m_num_bitsets ) {
          // lflank: L1 left flank size relative to centre
          auto lflank = ( L1_NUM_BITSETS >> 1 ) - 1;
          // rfit_bidx: right-most bitset index fitting the whole L1
          auto rfit_bidx = this->m_num_bitsets - L1_NUM_BITSETS;
          auto pb_bidx = ( ( ctr_bidx > lflank ) ? ( ctr_bidx - lflank ) : 0 );
          // for values of `centre` being closer to the end, l1 covers the last `TL1_Size` bits
          this->l1_begin_bidx = Kokkos::min( rfit_bidx, pb_bidx );
          this->l1_begin = HBitVector::start_index( this->l1_begin_bidx );
        }
        else {
          this->l1_begin_bidx = 0;
          this->l1_begin = 0;
        }

        auto l1size = HBitVector::l1_scratch_size();
        this->l1_data = ( view_type )
          ( tm.team_scratch( 0 ).get_shmem_aligned( l1size, space_alignment ) );
        auto l2size = this->l2_scratch_size();
        if ( l2size != 0 ) {
          this->l2_data = ( view_type )
            ( tm.team_scratch( 1 ).get_shmem_aligned( l2size, space_alignment ) );
        }
      }
      /* === STATIC MEMBERS === */
      /**
       *   @brief Get bitset index of the bit at `idx`
       *
       *   @param idx Bit index
       */
      static KOKKOS_INLINE_FUNCTION size_type
      bindex( size_type idx ) noexcept
      {
        return idx >> BINDEX_SHIFT;
      }

      /**
       *   @brief Get bitset offset of the bit at `idx`
       *
       *   @param idx Bit index
       */
      static KOKKOS_INLINE_FUNCTION size_type
      boffset( size_type idx ) noexcept
      {
        return idx & BOFFSET_MASK;
      }

      /**
       *   @brief Get index of the starting bit of the bitset with index `bidx`
       *
       *   @param bidx Bitset index
       */
      static KOKKOS_INLINE_FUNCTION size_type
      start_index( size_type bidx ) noexcept
      {
        return bidx << BINDEX_SHIFT;
      }

      /**
       *   @brief Get the left closest 'aligned index' relative to index `idx`
       *
       *   @param idx Bit index
       *
       *   Aligned index of an index, lets say `idx`, is the index of the
       *   starting bit of the bitset which includes the bit at `idx`.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      aligned_index( size_type idx ) noexcept
      {
        return idx & INDEX_ALIGN_MASK;
      }

      /**
       *   @brief Get the right closest 'aligned index' relative to index `idx`
       *
       *   @param idx Bit index
       *
       *   @return the index of the starting bit of the bitset immediately
       *   after the one including the bit at `idx` unless `idx` is itself a
       *   starting bit of a bitset.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      aligned_index_ceil( size_type idx ) noexcept
      {
        return ( idx + ( BITSET_WIDTH - 1 ) ) & ( INDEX_ALIGN_MASK );
      }

      /**
       *   @brief Return aligned size for a vector of `n` bits
       *
       *   The aligned size is the smallest multiple of `BITSET_WIDTH` which is
       *   larger than the actual size of the bitvector. The minimum aligned
       *   size is the size of L1.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      aligned_size( size_type n )
      {
        auto asize = ( n + ( BITSET_WIDTH - 1 ) ) & ( INDEX_ALIGN_MASK );  // aligned_index_ceil( n );
        return Kokkos::max( asize, TL1_Size );
      }

      /**
       *   @brief Return the number of allocated bitsets for L1
       *
       *   NOTE: The vector itself might occupy less bitsets than the number of
       *   allocated ones as the L1 allocated size is fixed in compile time.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l1_num_bitsets( ) noexcept
      {
        return L1_NUM_BITSETS;
      }

      /**
       *   @brief Return the allocated size of L1 in bits
       *
       *   NOTE: The vector itself might occupy less bits than the value return
       *   by this function as the L1 allocated size is fixed in compile time.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l1_size( ) noexcept
      {
        return TL1_Size;
      }

      /**
       *   @brief Return the allocated size of L1 in bytes
       *
       *   NOTE: The vector itself might occupy less bytes than the value
       *   return by this function as the L1 allocated size is fixed in compile
       *   time.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l1_scratch_size( ) noexcept
      {
        return L1_SIZE_BYTES;
      }

      /**
       *   @brief Return the required number of bitsets for a vector of `n` bits
       */
      static KOKKOS_INLINE_FUNCTION size_type
      num_bitsets( size_type n ) noexcept
      {
        auto x_size = HBitVector::aligned_size( n );
        return HBitVector::bindex( x_size );
      }

      /**
       *   @brief Return the required number of bitsets in L2 for a vector of
       *          `n` bits
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l2_num_bitsets( size_type n ) noexcept
      {
        auto nbitsets = HBitVector::num_bitsets( n );
        return ( nbitsets > HBitVector::l1_num_bitsets() )
                   ? nbitsets - HBitVector::l1_num_bitsets()
                   : 0;
      }

      /**
       *   @brief Return the required L2 size in bits for a vector of `n` bits
       *
       *   NOTE: The vector itself might occupy less bits in L2 than the value
       *   return by this function which is the size need to be *allocated*.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l2_size( size_type n ) noexcept
      {
        auto x_size = HBitVector::aligned_size( n );
        return ( x_size > HBitVector::l1_size() )
                   ? x_size - HBitVector::l1_size()
                   : 0;
      }

      /**
       *   @brief Return the required L2 size in bytes for a vector of `n` bits
       *
       *   NOTE: The vector itself might occupy less bytes in L2 than the value
       *   return by this function which is the size need to be *allocated*.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l2_scratch_size( size_type n ) noexcept
      {
        return HBitVector::l2_size( n ) / CHAR_BIT;
      }

      /**
       *   @brief Return the required number of bytes for a vector of `n` bits
       */
      static KOKKOS_INLINE_FUNCTION size_type
      capacity( size_type n ) noexcept
      {
        return HBitVector::aligned_size( n ) / CHAR_BIT;
      }

      template< typename TPolicy >
      static inline TPolicy
      set_scratch_size( TPolicy& policy, size_type n )
      {
        auto l1size = HBitVector::l1_scratch_size();
        auto l2size = HBitVector::l2_scratch_size( n );
        policy.set_scratch_size( 0, Kokkos::PerTeam( l1size ) );
        if ( l2size != 0 ) policy.set_scratch_size( 1, Kokkos::PerTeam( l2size ) );
        return policy;
      }

      static KOKKOS_INLINE_FUNCTION size_type
      cnt( bitset_type x ) noexcept
      {
        return Kokkos::Experimental::popcount_builtin( x );  // Requires Kokkos >=4.1.00
      }

      /**
       *   @brief Calculate the position of the i-th rightmost 1 bit in `x`
       *
       *   @param x Input bitset value
       *   @param i Argument i must be in the range [1..cnt(x)].
       */
      static KOKKOS_INLINE_FUNCTION size_type
      sel( bitset_type x, size_type i )
      {
        KOKKOS_IF_ON_DEVICE( ( return HBitVector::sel_device( x, i ); ) )
        KOKKOS_IF_ON_HOST( ( return sdsl::bits::sel( x, i ); ) )
      }

#if defined( KOKKOS_ENABLE_CUDA )
      /**
       *   @brief Calculate the position of the i-th rightmost 1 bit in `x`
       *          (on CUDA)
       *
       *   @param x Input bitset value
       *   @param i Argument i must be in the range [1..cnt(x)].
       */
      static KOKKOS_INLINE_FUNCTION size_type
      sel_device( bitset_type x, size_type i )
      {
        if constexpr ( BITSET_WIDTH == 64 ) {
          uint32_t lsw = x & 0xffffffff;
          auto cnt_lsw = HBitVector::cnt( lsw );
          if ( i <= cnt_lsw ) return __fns( lsw, 0, i );
          else {
            uint32_t msw = ( x >> 32 ) & 0xffffffff;
            return __fns( msw, 0, i - cnt_lsw );
          }
        }
        else {
          return __fns( x, 0, i );
        }
      }
#endif

      static KOKKOS_INLINE_FUNCTION bitset_type
      msb( bitset_type x ) noexcept
      {
        return x >> ( BITSET_WIDTH - 1 );
      }

      static KOKKOS_INLINE_FUNCTION bitset_type
      lsb( bitset_type x ) noexcept
      {
        return x % 2;
      }

      // Following methods are taken from sdsl::bits
      static KOKKOS_INLINE_FUNCTION size_type
      cnt10( bitset_type x, bitset_type c ) noexcept
      {
        return HBitVector::cnt( ( ( x << 1 ) | c ) & ( ~x ) );
      }

      static KOKKOS_INLINE_FUNCTION bitset_type
      map10( bitset_type x, bitset_type c ) noexcept
      {
        return ( ( ( x << 1 ) | c ) & ( ~x ) );
      }

      static KOKKOS_INLINE_FUNCTION size_type
      cnt01( bitset_type x, bitset_type c ) noexcept
      {
        return HBitVector::cnt( ( x ^ ( ( x << 1 ) | c ) ) & x );
      }

      static KOKKOS_INLINE_FUNCTION bitset_type
      map01( bitset_type x, bitset_type c ) noexcept
      {
        return ( ( x ^ ( ( x << 1 ) | c ) ) & x );
      }
      /* === OPERATORS === */
      /**
       *   @brief Get bitset by *global* bitset index.
       */
      KOKKOS_INLINE_FUNCTION bitset_type&
      operator()( size_type bidx ) const
      {
        auto r_bidx = this->relative_bitset( bidx );
        if ( r_bidx < HBitVector::l1_num_bitsets() ) {
          return this->l1_data[ r_bidx ];
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          return this->l2_data[ r_bidx ];
        }
      }

      KOKKOS_INLINE_FUNCTION value_type
      operator[]( size_type idx ) const
      {
        auto r_idx = this->relative_idx( idx );
        auto r_bidx = HBitVector::bindex( r_idx );
        auto offset = HBitVector::boffset( r_idx );
        if ( r_idx < HBitVector::l1_size() ) {
          return ( this->l1_data[ r_bidx ] >> offset ) /*& BITSET_ONE*/;
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          return ( this->l2_data[ r_bidx ] >> offset ) /*& BITSET_ONE*/;
        }
      }
      /* === ACCESSORS === */
      /*
      KOKKOS_INLINE_FUNCTION size_type
      size( ) const
      {
        return this->m_size;
      }
      */

      /**
       *   @brief Return aligned size of the bitvector
       *
       *   The aligned size is the smallest multiple of `BITSET_WIDTH` which is
       *   larger than the actual size of the bitvector.
       */
      KOKKOS_INLINE_FUNCTION size_type
      aligned_size( ) const
      {
        return this->m_x_size;
      }

      KOKKOS_INLINE_FUNCTION size_type
      num_bitsets( ) const
      {
        return this->m_num_bitsets;
      }

      KOKKOS_INLINE_FUNCTION size_type
      l1_begin_bindex( ) const
      {
        return this->l1_begin_bidx;
      }

      KOKKOS_INLINE_FUNCTION size_type
      l1_begin_idx( ) const
      {
        return this->l1_begin;
      }
      /* === METHODS === */
      /**
       *   @brief Returns the allocated size (L1+L2) in bytes
       */
      KOKKOS_INLINE_FUNCTION size_type
      capacity( ) const noexcept
      {
        return this->m_x_size / CHAR_BIT;
      }

      KOKKOS_INLINE_FUNCTION size_type
      l2_num_bitsets( ) const noexcept
      {
        auto nbitsets = this->num_bitsets();
        return ( nbitsets > HBitVector::l1_num_bitsets() )
                   ? nbitsets - HBitVector::l1_num_bitsets()
                   : 0;
      }

      /**
       *   @brief Return the allocated size of L2 in bits
       *
       *   NOTE: The vector itself might occupy less bits in L2 than the value
       *   return by this function which is the *allocated* size.
       */
      KOKKOS_INLINE_FUNCTION size_type
      l2_size( ) const noexcept
      {
        return ( this->m_x_size > HBitVector::l1_size() )
                   ? this->m_x_size - HBitVector::l1_size()
                   : 0;
      }

      /**
       *   @brief Return the allocated size of L2 in bytes
       *
       *   NOTE: The vector itself might occupy less bytes in L2 than the value
       *   return by this function which is the *allocated* size.
       */
      KOKKOS_INLINE_FUNCTION size_type
      l2_scratch_size( ) const noexcept
      {
        return this->l2_size() / CHAR_BIT;
      }

      /**
       *   @brief Return the relative index of `idx`
       *
       *   NOTE: The bitvector leaves the unused bit space [size, aln_size)
       *   untouched meaning that rearranging the bits for multi-level
       *   decompistion is done as if all allocated bits are used.
       */
      KOKKOS_INLINE_FUNCTION size_type
      relative_idx( size_type idx ) const noexcept
      {
        // NOTE: Using if-else is faster than computing modulo on both CPU and GPU
        //return ( this->m_x_size + idx - this->l1_begin ) % this->m_x_size;
        if ( this->l1_begin <= idx ) return idx - this->l1_begin;
        else return this->m_x_size + idx - this->l1_begin;
      }

      KOKKOS_INLINE_FUNCTION size_type
      relative_bitset( size_type bidx ) const noexcept
      {
        // NOTE: Using if-else is faster than computing modulo on both CPU and GPU
        //return ( this->m_num_bitsets + bidx - this->l1_begin_bidx )
        //       % this->m_num_bitsets;
        if ( this->l1_begin_bidx <= bidx ) return bidx - this->l1_begin_bidx;
        else return this->m_num_bitsets + bidx - this->l1_begin_bidx;
      }

      /**
       *   @brief Zero all bitsets in L1
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l1( const member_type& tm ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, HBitVector::l1_num_bitsets() ),
            [=]( const uint64_t j ) { this->l1_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero all bitsets in L2
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, this->l2_num_bitsets() ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets [ls_bidx, lf_bidx) on L2
       *
       *   @param tm  Team Policy member
       *   @param ls_bidx Local bitset index on L2 (inclusive)
       *   @param lf_bidx Local bitset index on L2 (exclusive)
       *
       *   NOTE: Both input indices are local bitset indices.
       *
       *   NOTE: Caller should make sure that the range does not span over L1.
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, size_type ls_bidx,
                size_type lf_bidx ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets [s_bidx, f_bidx) which occurs on L2
       *
       *   @param tm  Team Policy member
       *   @param s_bidx Global bitset index (inclusive)
       *   @param f_bidx Global bitset index (exclusive)
       *
       *   NOTE: Both input indices are global bitset indices (i.e. not local).
       *
       *   NOTE: Caller should make sure that the range does not span over L1.
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_bidx( const member_type& tm, size_type s_bidx,
                        size_type f_bidx ) noexcept
      {
        assert( f_bidx != 0 );

        auto rs_bidx = this->relative_bitset( s_bidx );
        auto rf_bidx = this->relative_bitset( f_bidx - 1 ) + 1;
        auto ls_bidx = rs_bidx - HBitVector::l1_num_bitsets();
        auto lf_bidx = rf_bidx - HBitVector::l1_num_bitsets();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets on L2 covering bit range [s_idx, f_idx)
       *
       *   @param tm  Team Policy member
       *   @param s_idx Global bit index in the vector (inclusive)
       *   @param f_idx Global bit index in the vector (exclusive)
       *
       *   NOTE: `s_idx` and `f_idx` are global *bit* indices (i.e. not local).
       *
       *   NOTE: Caller should make sure that the range does not span over L1.
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_idx( const member_type& tm, size_type s_idx,
                       size_type f_idx ) noexcept
      {
        assert( f_idx != 0 );

        auto rs_idx = this->relative_idx( s_idx );
        auto rf_idx = this->relative_idx( f_idx - 1 ) + 1;
        auto ls_idx = rs_idx - HBitVector::l1_size();
        auto lf_idx = rf_idx - HBitVector::l1_size();
        auto ls_bidx = HBitVector::bindex( ls_idx );
        auto lf_bidx = HBitVector::bindex( lf_idx );
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Set a bit in the vector (ThreadSequential)
       *
       *   NOTE: The function should be called by a single thread/lane once.
       *
       *   NOTE: All write access are NON-atomic.
       */
      KOKKOS_INLINE_FUNCTION void
      set( size_type idx, ThreadSequentialPartition ) noexcept
      {
        assert( idx < this->m_x_size );

        auto r_idx = this->relative_idx( idx );
        auto r_bidx = HBitVector::bindex( r_idx );
        auto mask = BITSET_ONE << HBitVector::boffset( r_idx );

        if ( r_idx < HBitVector::l1_size() ) {  // most probable
          this->l1_data[ r_bidx ] |= mask;
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          this->l2_data[ r_bidx ] |= mask;
        }
      }

      /**
       *   @brief Set a bit in the vector
       *          (TeamSequential/ThreadParallel/TeamFlatParallel)
       *
       *   NOTE: The function should be called by a single thread to be run on
       *         a vector lane (sequential).
       *
       *   NOTE: All write access are atomic.
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION std::enable_if_t<
          !std::is_same< TSpec, ThreadSequentialTag >::value >
      set( size_type idx, ExecPartition< TSpec > ) noexcept
      {
        assert( idx < this->m_x_size );

        auto r_idx = this->relative_idx( idx );
        auto r_bidx = HBitVector::bindex( r_idx );
        auto mask = BITSET_ONE << HBitVector::boffset( r_idx );

        if ( r_idx < HBitVector::l1_size() ) {  // most probable
          Kokkos::atomic_or( &this->l1_data[ r_bidx ], mask );
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          Kokkos::atomic_or( &this->l2_data[ r_bidx ], mask );
        }
      }

      /**
       *   @brief Set a range of bits in the vector (TeamSequential)
       *
       *   NOTE: The function should be called by a team to be run by a single
       *         thread (vector parallelism).
       *
       *   NOTE: This function assumes that no other thread setting any bits in
       *         range [s_idx, f_idx] except for bitsets at end points.
       */
      KOKKOS_INLINE_FUNCTION void
      set( const member_type& tm, size_type s_idx, size_type f_idx,
           TeamSequentialPartition tag ) noexcept
      {
        assert( s_idx <= f_idx );
        assert( f_idx < this->m_x_size );

        if ( s_idx == f_idx ) {
          Kokkos::single( Kokkos::PerThread( tm ), [=]() {
            this->set( s_idx, tag );
          } );
          return;
        }

        auto rs_idx = this->relative_idx( s_idx );
        auto rf_idx = this->relative_idx( f_idx );

        auto setbits =
          [&tm]( auto data_ptr, auto ls_idx, auto lf_idx ) {
            auto ls_bidx = HBitVector::bindex( ls_idx );
            auto lf_bidx = HBitVector::bindex( lf_idx );

            if ( ls_bidx != lf_bidx ) {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto mask = ( BITSET_ALL_SET << s_offset );
                Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              } );

              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange( tm, ls_bidx + 1, lf_bidx ),
                  [=]( const uint64_t k ) {
                    data_ptr[ k ] |= BITSET_ALL_SET;
                  } );

              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ALL_SET >> ( BITSET_WIDTH - 1u - f_offset ) );
                Kokkos::atomic_or( &data_ptr[ lf_bidx ], mask );
              } );
            }
            else {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ONE << ( f_offset - s_offset ) );
                mask = ( ( ( mask << 1 ) - 1 ) << s_offset );
                Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              } );
            }
          };

        if ( rs_idx < l1_size() && rf_idx < l1_size() ) {  // the range is in L1
          setbits( this->l1_data, rs_idx, rf_idx );
        }
        else if ( rs_idx < l1_size() && l1_size() <= rf_idx ) {
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, rs_idx, l1_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rs_idx <= rf_idx ) {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l2_data, ls_idx, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rf_idx < l1_size() ) {
          auto ls_idx = rs_idx - l1_size();
          setbits( this->l1_data, 0, rf_idx );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
        }
        else {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, 0, l1_size() - 1 );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
      }

      KOKKOS_INLINE_FUNCTION void
      set( size_type s_idx, size_type f_idx )
      {
        assert( s_idx <= f_idx );
        assert( f_idx < this->m_x_size );

        if ( s_idx == f_idx )
          return this->set( f_idx, TeamSequentialPartition{} );

        auto rs_idx = this->relative_idx( s_idx );
        auto rf_idx = this->relative_idx( f_idx );

        auto setbits =
          []( auto data_ptr, auto ls_idx, auto lf_idx ) {
            auto ls_bidx = HBitVector::bindex( ls_idx );
            auto lf_bidx = HBitVector::bindex( lf_idx );
            auto s_offset = HBitVector::boffset( ls_idx );
            auto f_offset = HBitVector::boffset( lf_idx );

            if ( ls_bidx != lf_bidx ) {
              auto mask = ( BITSET_ALL_SET << s_offset );
              Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              for ( auto i = ls_bidx + 1; i < lf_bidx; ++i ) data_ptr[ i ] |= BITSET_ALL_SET;
              mask = ( BITSET_ALL_SET >> ( BITSET_WIDTH - 1u - f_offset ) );
              Kokkos::atomic_or( &data_ptr[ lf_bidx ], mask );
            }
            else {
              auto mask = ( BITSET_ONE << ( f_offset - s_offset ) );
              mask = ( ( ( mask << 1 ) - 1 ) << s_offset );
              Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
            }
          };

        if ( rs_idx < l1_size() && rf_idx < l1_size() ) {  // the range is in L1
          setbits( this->l1_data, rs_idx, rf_idx );
        }
        else if ( rs_idx < l1_size() && l1_size() <= rf_idx ) {
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, rs_idx, l1_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rs_idx <= rf_idx ) {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l2_data, ls_idx, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rf_idx < l1_size() ) {
          auto ls_idx = rs_idx - l1_size();
          setbits( this->l1_data, 0, rf_idx );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
        }
        else {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, 0, l1_size() - 1 );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
      }
  };
}  /* --- end of namespace psi --- */


#endif // PSI_HBITVECTOR_HPP_
