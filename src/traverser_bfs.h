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
            unsigned int len, vg::Position s )
          : TBase( graph, r, index, len, s )
        { }

        TraverserBFS( const VarGraph* graph, const records_type* r, TIndex* index,
            unsigned int len )
          : TBase( graph, r, index, len )
        { }
        /* ====================  METHODS       ======================================= */
          inline void
        run( std::function< void( output_type const& ) >& callback )
        {
          this->frontier_states.push_back( {
                iterator_type( *( this->reads_index ) ),
                1,                     /**< @brief One more than allowed mismatches. */
                this->start_locus } );

          while( true ) {
            std::size_t current_size = this->frontier_states.size();
            for( std::size_t i = 0; i < current_size; ++i ) {
              filter( i, callback );
            }
            for( std::size_t i = 0; i < current_size; ++i ) {
              advance( i );
            }
            bool tie = true;
            for( std::size_t i = 0; i < current_size; ++i ) {
              if ( compute( i ) ) {
                tie = false;
              }
            }
            if ( tie ) {
              break;
            }
          }

          this->frontier_states.clear();
        }

          inline void
        filter( std::size_t state_idx,
            std::function< void( output_type const& ) >& callback )
        {
          assert( state_idx < this->frontier_states.size() );
          auto& state = this->frontier_states[ state_idx ];
          if ( state.mismatches != 0 && rep_length( state.iter ) == this->seed_len ) {
            // Cross out the state.
            state.mismatches = 0;
            // Process the seed hit.
            seqan::String< TSAValue > saPositions = getOccurrences( state.iter.get_iter_() );
            for ( unsigned i = 0; i < length( saPositions ); ++i )
            {
              output_type hit;
              hit.node_id = this->start_locus.node_id();
              hit.node_offset = this->start_locus.offset();
              auto id = position_to_id( *(this->reads), saPositions[i].i1 );
              hit.read_id = id;                // Read ID.
              hit.read_offset = saPositions[i].i2;  // Position in the read.
              callback( hit );
            }
          }
        }

          inline bool
        compute( std::size_t state_idx )
        {
          assert( state_idx < this->frontier_states.size() );
          auto& state = this->frontier_states[ state_idx ];
          if ( state.mismatches == 0 )
          {
            return false;
          }
          const auto& sequence = this->vargraph->node_sequence( state.pos.node_id() );
          const auto& itr_replen = rep_length( state.iter );

          assert( itr_replen < this->seed_len );
          auto remainder = this->seed_len - itr_replen;
          TIterRawText< typename iterator_type::base_type > partseq;
          partseq = sequence.substr( state.pos.offset(), remainder );

          typename seqan::Size< decltype( partseq ) >::Type i = 0;
          for ( ; i < seqan::length( partseq ); ++i ) {
            if ( partseq[i] == 'N' || !go_down( state.iter, partseq[i] ) ) {
              state.mismatches--;
              break;
            }
            stats_type::inc_total_nof_godowns();
          }

          state.pos.set_offset( state.pos.offset() + i );
          return true;
        }

          inline void
        advance( std::size_t state_idx )
        {
          assert( state_idx < this->frontier_states.size() );
          auto& state = this->frontier_states[ state_idx ];
          VarGraph::offset_type seqlen =
            this->vargraph->node_length( state.pos.node_id() );
          if ( state.mismatches == 0 || state.pos.offset() != seqlen )
          {
            return;
          }
          auto edges = this->vargraph->edges_from( state.pos.node_id() );
          auto it = edges.begin();
          if ( it == edges.end() ) {
            state.mismatches = 0;
            return;
          }
          state.pos.set_node_id( (*it).to() );
          state.pos.set_offset( 0 );
          ++it;

          for ( ; it != edges.end(); ++it )
          {
            typename traits_type::TState new_state = state;
            new_state.pos.set_node_id( (*it).to() );
            new_state.pos.set_offset( 0 );
            this->frontier_states.push_back( new_state );
          }
        }
    };  /* ----------  end of template class TraverserBFS  ---------- */
}  /* -----  end of namespace grem  ----- */

#endif  // end of TRAVERSER_BFS_H__
