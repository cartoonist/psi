/**
 *    @file  traverser_dfs.h
 *   @brief  DFS Traverser template class specialization.
 *
 *  The `Traverser` template class specialization for DFS traversal strategy.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Mon Sep 04, 2017  03:41
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef TRAVERSER_DFS_H__
#define TRAVERSER_DFS_H__

#include "traverser_base.h"

namespace grem {
  /**
   *  @brief  DFS Traverser.
   *
   *  Traverser a graph from a starting point and find seed hits using the read index.
   */
  template< typename TIndex,
    template<typename> class TMatchingTraits,
    typename TStatSpec >
      class TraverserDFS;

  template< typename TIndex, typename TStatSpec >
    class TraverserDFS< TIndex, ExactMatching, TStatSpec >
    : public TraverserBase< TIndex, DFS, ExactMatching, TStatSpec >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TraverserBase< TIndex, DFS, ExactMatching, TStatSpec > TBase;
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
        TraverserDFS( const VarGraph* graph, const records_type* r, TIndex* index,
            unsigned int len, vg::Position s )
          : TBase( graph, r, index, len, s ), cstate( {
              iterator_type( *( index ) ),
              0,
              this->start_locus } ) /**< @brief Should be initialised in `run`. */
        { }

        TraverserDFS( const VarGraph* graph, const records_type* r, TIndex* index,
            unsigned int len )
          : TBase( graph, r, index, len ), cstate( {
              iterator_type( *( index ) ),
              0,
              this->start_locus } ) /**< @brief Should be initialised in `run`. */
        { }
        /* ====================  METHODS       ======================================= */
          inline void
        run( std::function< void( output_type const& ) >& callback )
        {
          /* `cstate` should be updated in each `run` call with new `start_locus` value */
          this->frontier_states.push_back( {
              iterator_type( *( this->reads_index ) ),
              1,                     /**< @brief One more than allowed mismatches. */
              this->start_locus } );

          while( ! this->frontier_states.empty() )
          {
            filter( callback );
            advance( );
            compute( );
          }
        }

          inline void
        filter( std::function< void( output_type const& ) >& callback )
        {
          if ( cstate.mismatches != 0 && rep_length( cstate.iter ) == this->seed_len ) {
            // Cross out the cstate.
            cstate.mismatches = 0;
            // Process the seed hit.
            seqan::String< TSAValue > saPositions = getOccurrences( cstate.iter.get_iter_() );
            typename seqan::Size< decltype( saPositions ) >::Type i;
            for ( i = 0; i < length( saPositions ); ++i )
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
        compute( )
        {
          if ( cstate.mismatches == 0 ) return false;

          const auto& sequence = this->vargraph->node_sequence( cstate.pos.node_id() );
          const auto& itr_replen = rep_length( cstate.iter );
          assert( itr_replen < this->seed_len );
          VarGraph::offset_type end_idx =
            std::min( cstate.pos.offset() + this->seed_len - itr_replen, sequence.size() );
          VarGraph::offset_type i;
          for ( i = cstate.pos.offset(); i < end_idx; ++i ) {
            if ( sequence[i] == 'N' || !go_down( cstate.iter, sequence[i] ) ) {
              cstate.mismatches--;
              break;
            }
            stats_type::inc_total_nof_godowns();
          }

          cstate.pos.set_offset( i );
          return true;
        }

          inline void
        advance( )
        {
          if ( cstate.mismatches == 0 && ( !this->frontier_states.empty() ) )
          {
            this->cstate = this->frontier_states.back();
            this->frontier_states.pop_back();
            return;
          }

          VarGraph::offset_type seqlen =
            this->vargraph->node_length( cstate.pos.node_id() );
          if ( cstate.mismatches == 0 || cstate.pos.offset() != seqlen ) return;

          auto edges = this->vargraph->edges_from( this->cstate.pos.node_id() );
          auto it = edges.begin();
          if ( it == edges.end() ) {
            cstate.mismatches = 0;
            return;
          }
          cstate.pos.set_node_id( (*it).to() );
          cstate.pos.set_offset( 0 );
          ++it;

          for ( ; it != edges.end(); ++it )
          {
            typename traits_type::TState new_state = cstate;
            new_state.pos.set_node_id( (*it).to() );
            new_state.pos.set_offset( 0 );
            this->frontier_states.push_back( std::move( new_state ) );
          }
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        typename traits_type::TState cstate;
    };  /* ----------  end of template class TraverserDFS  ---------- */
}  /* -----  end of namespace grem  ----- */

#endif  // end of TRAVERSER_DFS_H__
