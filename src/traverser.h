/**
 *    @file  traverser.h
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

#ifndef TRAVERSER_H__
#define TRAVERSER_H__

#include "traverser_bfs.h"
#include "traverser_dfs.h"

namespace grem {
  template< typename TIndexSpec,
    typename TStrategy,
    template<typename> class TMatchingTraits,
    typename TStatTag = void >
    class Traverser;

  template< typename IndexSpec, typename TStatTag >
    class Traverser< IndexSpec, BFS, ExactMatching, TStatTag > {
      public:
        typedef TraverserBFS< IndexSpec, ExactMatching, TStatTag > Type;
    };  /* ----------  end of template class Traverser  ---------- */
}  /* -----  end of namespace grem  ----- */

#endif  // end of TRAVERSER_H__
