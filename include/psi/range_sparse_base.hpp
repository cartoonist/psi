/**
 *    @file  range_sparse_base.hpp
 *   @brief  Base header for range sparse matrix operations module
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu Sep 07, 2023  21:37
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_RANGE_SPARSE_BASE_HPP_
#define PSI_RANGE_SPARSE_BASE_HPP_

#include <cstddef>
#include <stdexcept>
#include <algorithm>

#include <Kokkos_Core.hpp>
#include <Kokkos_NestedSort.hpp>


namespace psi {
  // Accumulator tag
  template< typename TSpec >
  struct Accumulator {
    using spec_type = TSpec;
  };

  // Accumulator specialisation tag
  struct BTreeTag {};

  template< unsigned int TL1Size >
  struct HBitVectorTag {
    static constexpr const unsigned int L1SizeParam = TL1Size;
  };

  using NoAccumulator = Accumulator< void >;
  using BTreeAccumulator = Accumulator< BTreeTag >;

  template< unsigned int TL1Size=2048 >
  using HBitVectorAccumulator = Accumulator< HBitVectorTag< TL1Size > >;

  // Partition tag
  template< typename TSpec >
  struct ExecPartition {
    using type = TSpec;
  };

  // Partition specialisation tag
  struct ThreadRangeTag {};
  struct ThreadSequentialTag {};
  struct TeamSequentialTag {};
  struct ThreadParallelTag {};
  struct TeamFlatParallelTag {};

  using ThreadRangePartition = ExecPartition< ThreadRangeTag >;
  using ThreadSequentialPartition = ExecPartition< ThreadSequentialTag >;
  using TeamSequentialPartition = ExecPartition< TeamSequentialTag >;
  using ThreadParallelPartition = ExecPartition< ThreadParallelTag >;
  using TeamFlatParallelPartition = ExecPartition< TeamFlatParallelTag >;

  // Supported execution space by accumulators
  template< typename TAccumulator >
  struct AccumulatorExecSpace {
    using type = Kokkos::DefaultExecutionSpace;      // By default run on device
  };

  template<>
  struct AccumulatorExecSpace< BTreeAccumulator > {
    using type = Kokkos::DefaultHostExecutionSpace;  // Can only be executed on host
  };

  // Default partitioning by accumulators
  template< typename TAccumulator >
  struct AccumulatorDefaultPartition;

  template< >
  struct AccumulatorDefaultPartition< BTreeAccumulator >
  {
    using type = ThreadRangePartition;
  };

  template< unsigned int TL1Size >
  struct AccumulatorDefaultPartition< HBitVectorAccumulator< TL1Size > >
  {
    using type = TeamSequentialPartition;
  };

  /**
   *   @brief Suggested grid dimensions as a form of config class.
   */
  template< typename TExecSpace >
  struct SuggestedExecGrid {
    static constexpr inline int
    vector_size( )
    {
      return 1;
    }

    static constexpr inline int
    vector_size( const int row_density )
    {
      return 1;
    }

    static constexpr inline int
    vector_size( const int nnz, const int nr )
    {
      return 1;
    }

    static constexpr inline int
    team_size( )
    {
      return 1;
    }

    static constexpr inline int
    team_size( const int vector_size )
    {
      return 1;
    }

    static constexpr inline int
    team_work_size( )
    {
      return 16;
    }

    static constexpr inline int
    team_work_size( const int team_size )
    {
      return 16;
    }
  };

  template< typename TExecSpace >
  struct AutoExecGrid {
    static constexpr inline auto
    vector_size( )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    vector_size( const int row_density )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    vector_size( const int nnz, const int nr )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( const int vector_size )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline int
    team_work_size( )
    {
      return 16;
    }

    static constexpr inline int
    team_work_size( const int team_size )
    {
      return 16;
    }
  };

  #if defined(KOKKOS_ENABLE_CUDA)
  template< >
  struct SuggestedExecGrid< Kokkos::Cuda > {
    /* === STATIC MEMBERS === */
    static constexpr const int MAX_VECTOR_SIZE = 32;
    /* === STATIC METHODS === */
    static constexpr inline int
    row_density( const std::size_t nnz, const std::size_t nr )
    {
      if ( nr > 0 ) return nnz / double( nr ) + 0.5;
      return 1;
    }

    static constexpr inline int
    vector_size( const int rdense )
    {
      int vsize = rdense;
      if ( vsize < 3 ) {
        vsize = 2;
      } else if ( vsize <= 6 ) {
        vsize = 4;
      } else if ( vsize <= 12 ) {
        vsize = 8;
      } else if ( vsize <= 24 ) {
        vsize = 16;
      } else if ( vsize <= 48 ) {
        vsize = 32;
      } else {
        vsize = 64;
      }
      vsize = Kokkos::min( vsize, MAX_VECTOR_SIZE );
      return vsize;
    }

    static constexpr inline int
    vector_size( const std::size_t nnz, const std::size_t nr )
    {
      return vector_size( row_density( nnz, nr ) );
    }

    static constexpr inline int
    team_size( const int vector_size )
    {
      // TODO: where this is used, tune the target value for
      // threads per block (but 256 is probably OK for CUDA and HIP)
      return 256 / vector_size;
    }

    static constexpr inline int
    team_work_size( const int team_size )
    {
      return team_size;
    }
  };

  template< >
  struct AutoExecGrid< Kokkos::Cuda > {
    static constexpr inline auto
    vector_size( )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    vector_size( const int row_density )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    vector_size( const int nnz, const int nr )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( const int vector_size )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline int
    team_work_size( const int team_size )
    {
      return team_size;
    }
  };
  #endif

  // Configuration tag
  template< typename TAccumulator,
      typename TPartition = typename AccumulatorDefaultPartition< TAccumulator >::type,
      template< typename > typename TGrid = AutoExecGrid >
  struct SparseConfig {
    using partition_type = TPartition;
    using accumulator_type = TAccumulator;
    using execution_space = typename AccumulatorExecSpace< accumulator_type >::type;
    using grid_type = TGrid< execution_space >;
    /* === MEMBERS === */
    partition_type   part;
    accumulator_type accm;
    execution_space  space;
    grid_type        grid;
  };

  using DefaultSparseConfiguration = SparseConfig< BTreeAccumulator >;

  template< typename TRCRSMatrix >
  struct SparseRangeHandle {
    /* === MEMBERS TYPES === */
    using ordinal_type = typename TRCRSMatrix::ordinal_type;
    using size_type = typename TRCRSMatrix::size_type;
    /* === LIFECYCLE === */
    SparseRangeHandle() = delete;

    SparseRangeHandle( TRCRSMatrix const& a, TRCRSMatrix const& b )
    {
      this->a_ncols = a.numCols();
      this->b_ncols = b.numCols();
    }
    /* === DATA MEMBERS === */
    ordinal_type a_ncols;
    ordinal_type b_ncols;
  };

  template<
      typename TRowMapView,
      typename TEntriesView,
      typename TExecSpace=typename TRowMapView::execution_space >
  struct SortEntriesFunctor {
    /* === TYPE MEMBERS === */
    using policy_type = Kokkos::RangePolicy< TExecSpace >;
    using size_type = typename TRowMapView::non_const_value_type;
    /* === STATIC ASSERTS === */
    static_assert(
        std::is_same< typename TRowMapView::execution_space,
                      TExecSpace >::value,
        "execution space parameter is incompatible with input views" );
    static_assert(
        std::is_same< typename TEntriesView::execution_space,
                      TExecSpace >::value,
        "execution space parameter is incompatible with input views" );
    /* === DATA MEMBERS === */
    TRowMapView row_map;
    TEntriesView entries;
    /* === LIFE CYCLE === */
    SortEntriesFunctor( TRowMapView r, TEntriesView e )
      : row_map( r ), entries( e )
    { }
    /* === METHODS === */
    inline auto
    policy( const size_type nrows ) const
    {
      return policy_type( 0, nrows );
    }
    /* === OPERATORS === */
    KOKKOS_INLINE_FUNCTION void
    operator()( const uint64_t i ) const
    {
      auto begin = this->entries.data() + this->row_map( i );
      auto end = this->entries.data() + this->row_map( i + 1 );
      std::sort( begin, end );
    }
  };

#if defined(KOKKOS_ENABLE_CUDA)
  template< typename TRowMapView, typename TEntriesView >
  struct SortEntriesFunctor< TRowMapView, TEntriesView, Kokkos::Cuda > {
    /* === TYPE MEMBERS === */
    using execution_space = Kokkos::Cuda;
    using policy_type = Kokkos::TeamPolicy< execution_space >;
    using member_type = typename policy_type::member_type;
    using size_type = typename TRowMapView::non_const_value_type;
    /* === STATIC ASSERTS === */
    static_assert(
        std::is_same< typename TRowMapView::execution_space,
                      execution_space >::value,
        "execution space parameter is incompatible with input views" );
    static_assert(
        std::is_same< typename TEntriesView::execution_space,
                      execution_space >::value,
        "execution space parameter is incompatible with input views" );
    /* === DATA MEMBERS === */
    TRowMapView row_map;
    TEntriesView entries;
    /* === LIFE CYCLE === */
    SortEntriesFunctor( TRowMapView r, TEntriesView e )
      : row_map( r ), entries( e )
    { }
    /* === METHODS === */
    inline auto
    policy( const size_type nrows ) const
    {
      return policy_type( nrows, Kokkos::AUTO );
    }
    /* === OPERATORS === */
    KOKKOS_INLINE_FUNCTION void
    operator()( member_type const& tm ) const
    {
      auto i = tm.league_rank();
      auto l = this->row_map( i );
      auto u = this->row_map( i + 1 );
      assert( u >= l );
      if ( u != l ) {
        auto subview = Kokkos::subview( this->entries, Kokkos::pair( l, u ) );
        Kokkos::Experimental::sort_team( tm, subview );
      }
    }
  };
#endif

  /* borrowed from kokkos-kernel */
  template< typename TExecSpace >
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             TExecSpace )
  {
    throw std::runtime_error( "`get_suggested_vector_size` is not implemented for requested execution space" );
  }

  #if defined(KOKKOS_ENABLE_CUDA)
  inline int
  get_suggested_vector_size( const std::size_t nr,
                             const std::size_t nnz,
                             Kokkos::Cuda )
  {
    int suggested_vector_size_ = 1;
    int max_vector_size        = 32;
    if ( nr > 0 ) suggested_vector_size_ = nnz / double( nr ) + 0.5;
    if ( suggested_vector_size_ < 3 ) {
      suggested_vector_size_ = 2;
    } else if ( suggested_vector_size_ <= 6 ) {
      suggested_vector_size_ = 4;
    } else if ( suggested_vector_size_ <= 12 ) {
      suggested_vector_size_ = 8;
    } else if ( suggested_vector_size_ <= 24 ) {
      suggested_vector_size_ = 16;
    } else if ( suggested_vector_size_ <= 48 ) {
      suggested_vector_size_ = 32;
    } else {
      suggested_vector_size_ = 64;
    }
    if ( suggested_vector_size_ > max_vector_size )
      suggested_vector_size_ = max_vector_size;
    return suggested_vector_size_;
  }
  #endif

  #if defined(KOKKOS_ENABLE_OPENMP)
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             Kokkos::OpenMP )
  {
    //int suggested_vector_size_ = 1;
    //int max_vector_size        = 1;
    return /*suggested_vector_size_*/ 1;
  }
  #endif

  #if defined(KOKKOS_ENABLE_SERIAL)
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             Kokkos::Serial )
  {
    //int suggested_vector_size_ = 1;
    //int max_vector_size        = 1;
    return /*suggested_vector_size_*/ 1;
  }
  #endif

  #if defined(KOKKOS_ENABLE_THREADS)
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             Kokkos::Threads )
  {
    //int suggested_vector_size_ = 1;
    //int max_vector_size        = 1;
    return /*suggested_vector_size_*/ 1;
  }
  #endif

  template< typename TExecSpace, int TVectorSize >
  struct SuggestedTeamSize;

  template< typename TExecSpace >
  inline int get_suggested_team_size( const int vector_size, TExecSpace )
  {
    throw std::runtime_error( "`get_suggested_team_size` is not implemented for requested execution space" );
  }

  #if defined(KOKKOS_ENABLE_CUDA)
  inline int get_suggested_team_size( const int vector_size, Kokkos::Cuda )
  {
    // TODO: where this is used, tune the target value for
    // threads per block (but 256 is probably OK for CUDA and HIP)
    return 256 / vector_size;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::Cuda, TVectorSize > {
    constexpr static int value = 256 / TVectorSize;
  };
  #endif

  #if defined(KOKKOS_ENABLE_OPENMP)
  inline int get_suggested_team_size( const int vector_size, Kokkos::OpenMP )
  {
    return 1;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::OpenMP, TVectorSize > {
    constexpr static int value = 1;
  };
  #endif

  #if defined(KOKKOS_ENABLE_SERIAL)
  inline int get_suggested_team_size( const int vector_size, Kokkos::Serial )
  {
    return 1;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::Serial, TVectorSize > {
    constexpr static int value = 1;
  };
  #endif

  #if defined(KOKKOS_ENABLE_THREADS)
  inline int get_suggested_team_size( const int vector_size, Kokkos::Threads )
  {
    return 1;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::Threads, TVectorSize > {
    constexpr static int value = 1;
  };
  #endif

  template< typename TExecSpace, int TTeamSize >
  struct SuggestedTeamWorkSize;

  template< typename TExecSpace >
  inline int get_team_work_size( const int team_size, TExecSpace )
  {
    throw std::runtime_error( "`get_team_work_size` is not implemented for requested execution space" );
  }

  #if defined(KOKKOS_ENABLE_CUDA)
  inline int get_team_work_size( const int team_size, Kokkos::Cuda )
  {
    return team_size;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::Cuda, TTeamSize >
  {
    constexpr static int value = TTeamSize;
  };
  #endif

  #if defined(KOKKOS_ENABLE_OPENMP)
  inline int get_team_work_size( const int team_size, Kokkos::OpenMP )
  {
    return 16;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::OpenMP, TTeamSize >
  {
    constexpr static int value = 16;
  };
  #endif

  #if defined(KOKKOS_ENABLE_SERIAL)
  inline int get_team_work_size( const int team_size, Kokkos::Serial )
  {
    return 16;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::Serial, TTeamSize >
  {
    constexpr static int value = 16;
  };
  #endif

  #if defined(KOKKOS_ENABLE_THREADS)
  inline int get_team_work_size( const int team_size, Kokkos::Threads )
  {
    return 16;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::Threads, TTeamSize >
  {
    constexpr static int value = 16;
  };
  #endif
}  // namespace psi
#endif  // PSI_RANGE_SPARSE_BASE_HPP_
