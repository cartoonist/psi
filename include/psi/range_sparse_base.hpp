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

#include <Kokkos_Core.hpp>


namespace psi {
  // Accumulator tag
  template< typename TSpec >
  struct Accumulator {
    using spec_type = TSpec;
  };

  // Accumulator specialisation tag
  struct BTreeTag {};
  struct HBitVectorTag {};

  using NoAccumulator = Accumulator< void >;
  using BTreeAccumulator = Accumulator< BTreeTag >;
  using HBitVectorAccumulator = Accumulator< HBitVectorTag >;

  // Partition tag
  template< typename TSpec >
  struct ExecPartition {
    using type = TSpec;
  };

  // Partition specialisation tag
  struct ThreadRangeTag {};
  struct ThreadSequentialTag {};
  struct TeamSequentialTag {};
  // Unused specialisation tags yet
  //struct ThreadParallelTag {};
  //struct TeamFlatParallelTag {};

  using ThreadRangePartition = ExecPartition< ThreadRangeTag >;
  using ThreadSequentialPartition = ExecPartition< ThreadSequentialTag >;
  using TeamSequentialPartition = ExecPartition< TeamSequentialTag >;
  //using ThreadParallelPartition = ExecPartition< ThreadParallelTag >;
  //using TeamFlatParallelPartition = ExecPartition< TeamFlatParallelTag >;

  // Supported execution space by accumulators
  template< typename TAccumulator >
  struct AccumulatorExecSpace {
    using execution_space = Kokkos::DefaultExecutionSpace;      // By default run on device
  };

  template<>
  struct AccumulatorExecSpace< BTreeAccumulator > {
    using execution_space = Kokkos::DefaultHostExecutionSpace;  // Can only be executed on host
  };

  // Configuration tag
  template< typename TPartition, typename TAccumulator >
  struct SparseConfig {
    using partition_type = TPartition;
    using accumulator_type = TAccumulator;
    using execution_space = typename AccumulatorExecSpace< accumulator_type >::execution_space;
  };

  using DefaultSparseConfiguration = SparseConfig< ThreadRangePartition, BTreeAccumulator >;

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
