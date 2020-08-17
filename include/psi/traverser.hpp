/**
 *    @file  traverser.hpp
 *   @brief  Traverser main header file.
 *
 *  The main header file for Traversers.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Sep 04, 2017  03:29
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef PSI_TRAVERSER_HPP__
#define PSI_TRAVERSER_HPP__

#include "traverser_bfs.hpp"
#include "traverser_dfs.hpp"

namespace psi {
  template< typename TGraph,
    typename TIndex,
    typename TStrategy,
    template<typename, typename> class TMatchingTraits,
    typename TStatsSpec = WithStats >
    class Traverser;

  template< typename TGraph, typename TIndex, typename TStatsSpec >
    class Traverser< TGraph, TIndex, BFS, ExactMatching, TStatsSpec > {
      public:
        typedef TraverserBFS< TGraph, TIndex, ExactMatching, TStatsSpec > Type;
    };  /* ----------  end of template class Traverser  ---------- */

  template< typename TGraph, typename TIndex, typename TStatsSpec >
    class Traverser< TGraph, TIndex, DFS, ExactMatching, TStatsSpec > {
      public:
        typedef TraverserDFS< TGraph, TIndex, ExactMatching, TStatsSpec > Type;
    };  /* ----------  end of template class Traverser  ---------- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_TRAVERSER_HPP__ --- */
