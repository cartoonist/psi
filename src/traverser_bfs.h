/**
 *    @file  traverser_bfs.h
 *   @brief  BFS Traverser template class specialization.
 *
 *  The `Traverser` template class specialization for BFS traversal strategy.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Sep 04, 2017  03:24
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef TRAVERSER_BFS_H__
#define TRAVERSER_BFS_H__

#include "traverser_base.h"

namespace grem {
  /**
   *  @brief  BFS Traverser.
   *
   *  Traverser a graph from a starting point and find seed hits using the read index.
   */
  template< typename TIndex,
    template<typename> class TMatchingTraits,
    typename TStatSpec >
    class TraverserBFS;

  template< typename TIndex, typename TStatSpec >
    class TraverserBFS< TIndex, ExactMatching, TStatSpec >
    : public TraverserBase< TIndex, BFS, ExactMatching, TStatSpec >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TraverserBase< TIndex, BFS, ExactMatching, TStatSpec > TBase;
        typedef typename TBase::output_type output_type;
        typedef typename TBase::index_type index_type;
        typedef typename TBase::indexspec_type indexspec_type;
        typedef typename TBase::stringset_type stringset_type;
        typedef typename TBase::text_type text_type;
        typedef typename TBase::records_type records_type;
        typedef typename TBase::iterspec_type iterspec_type;
        typedef typename TBase::iterator_type iterator_type;
        typedef typename TBase::traits_type traits_type;
        typedef typename TBase::TSAValue TSAValue;
        typedef typename TBase::stats_type stats_type;
        /* ====================  LIFECYCLE     ======================================= */
        TraverserBFS( const VarGraph* graph, const records_type* r, TIndex* index,
            unsigned int len )
          : TBase( graph, r, index, len )
        { }

        TraverserBFS( const VarGraph* graph, unsigned int len )
          : TBase( graph, len )
        { }
        /* ====================  METHODS       ======================================= */
          inline void
        run( std::function< void( output_type const& ) > callback )
        {
          bool tie;
          do {
            std::size_t nofstates = this->states.size();
            tie = true;
            for( std::size_t idx = 0; idx < nofstates; ++idx ) {
              if ( this->states[ idx ].mismatches == 0 ) continue;
              filter( this->states[ idx ], callback );
              advance( this->states[ idx ] );
              if ( compute( this->states[ idx ] ) ) tie = false;
            }
          } while ( !tie );

          this->states.clear();
        }

          inline void
        filter( typename traits_type::TState& state,
            std::function< void( output_type const& ) > callback )
        {
          if ( state.mismatches != 0 && state.depth == this->seed_len ) {
            // Cross out the state.
            state.mismatches = 0;
            // Process the seed hit.
            seqan::String< TSAValue > saPositions = getOccurrences( state.iter.get_iter_() );
            typename seqan::Size< decltype( saPositions ) >::Type i;
            for ( i = 0; i < length( saPositions ); ++i )
            {
              output_type hit;
              hit.node_id = state.spos.node_id();
              hit.node_offset = state.spos.offset();
              hit.read_id = position_to_id( *(this->reads), saPositions[i].i1 );  // Read ID.
              hit.read_offset = position_to_offset( *(this->reads), saPositions[i] );  // Position in the read.
              callback( hit );
            }
          }
        }

          inline bool
        compute( typename traits_type::TState& state )
        {
          if ( state.mismatches == 0 ) return false;

          const auto& sequence = this->vargraph->node_sequence( state.cpos.node_id() );
          assert( state.depth < this->seed_len );
          std::make_unsigned_t< VarGraph::offset_type > end_idx =
            state.cpos.offset() + this->seed_len - state.depth;
          std::make_unsigned_t< VarGraph::offset_type > i;
          for ( i = state.cpos.offset(); i < end_idx && i < sequence.size(); ++i ) {
            if ( sequence[i] == 'N' || !go_down( state.iter, sequence[i] ) ) {
              state.mismatches--;
              break;
            }
            ++state.depth;
            stats_type::inc_total_nof_godowns();
          }

          state.cpos.set_offset( i );
          if ( i == sequence.size() ) state.end = true;
          return true;
        }

          inline void
        advance( typename traits_type::TState& state )
        {
          if ( state.mismatches == 0 || !state.end ) return;
          const auto& edges = this->vargraph->edges_from( state.cpos.node_id() );
          auto it = edges.begin();
          if ( it == edges.end() ) {
            state.mismatches = 0;
            return;
          }
          state.cpos.set_node_id( (*it).to() );
          state.cpos.set_offset( 0 );
          state.end = false;
          ++it;

          for ( ; it != edges.end(); ++it )
          {
            this->states.emplace_back( state );
            this->states.back().cpos.set_node_id( (*it).to() );
            this->states.back().cpos.set_offset( 0 );
          }
        }
    };  /* ----------  end of template class TraverserBFS  ---------- */
}  /* -----  end of namespace grem  ----- */

#endif  // end of TRAVERSER_BFS_H__
