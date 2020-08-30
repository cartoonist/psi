/**
 *    @file  traverser_bfs.hpp
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

#ifndef PSI_TRAVERSER_BFS_HPP__
#define PSI_TRAVERSER_BFS_HPP__

#include "traverser_base.hpp"

namespace psi {
  /**
   *  @brief  BFS Traverser.
   *
   *  Traverser a graph from a starting point and find seed hits using the read index.
   */
  template< class TGraph,
    typename TIndex,
    template<typename, typename> class TMatchingTraits,
    typename TStatsSpec >
    class TraverserBFS;

  template< class TGraph, typename TIndex, typename TStatsSpec >
    class TraverserBFS< TGraph, TIndex, ExactMatching, TStatsSpec >
    : public TraverserBase< TGraph, TIndex, BFS, ExactMatching, TStatsSpec >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TraverserBase< TGraph, TIndex, BFS, ExactMatching, TStatsSpec > base_type;
        typedef typename base_type::graph_type graph_type;
        typedef typename graph_type::id_type id_type;
        typedef typename graph_type::offset_type offset_type;
        typedef typename graph_type::linktype_type linktype_type;
        typedef typename base_type::output_type output_type;
        typedef typename base_type::index_type index_type;
        typedef typename base_type::indexspec_type indexspec_type;
        typedef typename base_type::stringset_type stringset_type;
        typedef typename base_type::text_type text_type;
        typedef typename base_type::records_type records_type;
        typedef typename base_type::iterspec_type iterspec_type;
        typedef typename base_type::iterator_type iterator_type;
        typedef typename base_type::traits_type traits_type;
        typedef typename base_type::TSAValue TSAValue;
        typedef typename base_type::stats_type stats_type;
        /* ====================  LIFECYCLE     ======================================= */
        TraverserBFS( const graph_type* g, const records_type* r, TIndex* index,
            unsigned int len )
          : base_type( g, r, index, len )
        { }

        TraverserBFS( const graph_type* g, unsigned int len )
          : base_type( g, len )
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
            stats_type::inc_total_seeds_off_paths( length( saPositions ) );
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

          const auto& sequence = this->graph_ptr->node_sequence( state.cpos.node_id() );
          assert( state.depth < this->seed_len );
          offset_type end_idx = state.cpos.offset() + this->seed_len - state.depth;
          offset_type i;
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
          if ( !this->graph_ptr->has_edges_out( state.cpos.node_id() ) ) {
            state.mismatches = 0;
            return;
          }
          bool first = true;
          this->graph_ptr->for_each_edges_out(
              state.cpos.node_id(),
              [this, &state, &first]( id_type to, linktype_type ) {
                if ( first ) {
                  state.cpos.set_node_id( to );
                  state.cpos.set_offset( 0 );
                  state.end = false;
                  first = false;
                  return true;
                }
                this->states.push_back( state );
                this->states.back().cpos.set_node_id( to );
                this->states.back().cpos.set_offset( 0 );
                return true;
              } );
        }
    };  /* --- end of template class TraverserBFS --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_TRAVERSER_BFS_HPP__ --- */
