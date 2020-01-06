/**
 *    @file  traverser_dfs.hpp
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

#ifndef PSI_TRAVERSER_DFS_HPP__
#define PSI_TRAVERSER_DFS_HPP__

#include "traverser_base.hpp"

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
            unsigned int len )
          : TBase( graph, r, index, len ), cstate( index, 0, 0, 0, 0 )
        { }

        TraverserDFS( const VarGraph* graph, unsigned int len )
          : TBase( graph, len ), cstate( nullptr, 0, 0, 0, 0 )
        { }
        /* ====================  METHODS       ======================================= */
          inline void
        run( std::function< void( output_type const& ) > callback )
        {
          while( ! this->states.empty() )
          {
            filter( callback );
            advance( );
            compute( );
          }
        }

          inline void
        filter( std::function< void( output_type const& ) > callback )
        {
          if ( cstate.mismatches != 0 && cstate.depth == this->seed_len ) {
            // Cross out the cstate.
            cstate.mismatches = 0;
            // Process the seed hit.
            seqan::String< TSAValue > saPositions = getOccurrences( cstate.iter.get_iter_() );
            typename seqan::Size< decltype( saPositions ) >::Type i;
            for ( i = 0; i < length( saPositions ); ++i )
            {
              output_type hit;
              hit.node_id = cstate.spos.node_id();
              hit.node_offset = cstate.spos.offset();
              hit.read_id = position_to_id( *(this->reads), saPositions[i].i1 );  // Read ID.
              hit.read_offset = position_to_offset( *(this->reads), saPositions[i] );  // Position in the read.
              callback( hit );
            }
          }
        }

          inline bool
        compute( )
        {
          if ( cstate.mismatches == 0 ) return false;

          const auto& sequence = this->vargraph->node_sequence( cstate.cpos.node_id() );
          assert( cstate.depth < this->seed_len );
          std::make_unsigned_t< VarGraph::offset_type > end_idx =
            cstate.cpos.offset() + this->seed_len - cstate.depth;
          std::make_unsigned_t< VarGraph::offset_type > i;
          for ( i = cstate.cpos.offset(); i < end_idx && i < sequence.size(); ++i ) {
            if ( sequence[i] == 'N' || !go_down( cstate.iter, sequence[i] ) ) {
              cstate.mismatches--;
              break;
            }
            ++cstate.depth;
            stats_type::inc_total_nof_godowns();
          }

          cstate.cpos.set_offset( i );
          if ( i == sequence.size() ) cstate.end = true;
          return true;
        }

          inline void
        advance( )
        {
          if ( cstate.mismatches == 0 && ( !this->states.empty() ) )
          {
            this->cstate = this->states.back();
            this->states.pop_back();
            return;
          }

          if ( cstate.mismatches == 0 || !cstate.end ) return;

          const auto& edges = this->vargraph->edges_from( this->cstate.cpos.node_id() );
          auto it = edges.begin();
          if ( it == edges.end() ) {
            cstate.mismatches = 0;
            return;
          }
          cstate.cpos.set_node_id( (*it).to() );
          cstate.cpos.set_offset( 0 );
          cstate.end = false;
          ++it;

          for ( ; it != edges.end(); ++it )
          {
            this->states.emplace_back( cstate );
            this->states.back().cpos.set_node_id( (*it).to() );
            this->states.back().cpos.set_offset( 0 );
          }
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        typename traits_type::TState cstate;
    };  /* ----------  end of template class TraverserDFS  ---------- */
}  /* -----  end of namespace grem  ----- */

#endif  /* --- #ifndef PSI_TRAVERSER_DFS_HPP__ --- */
