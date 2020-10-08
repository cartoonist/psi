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

namespace psi {
  /**
   *  @brief  DFS Traverser.
   *
   *  Traverser a graph from a starting point and find seed hits using the read index.
   */
  template< class TGraph,
    typename TIndex,
    template<typename, typename> class TMatchingTraits,
    typename TStatsSpec >
      class TraverserDFS;

  template< class TGraph, typename TIndex, typename TStatsSpec >
    class TraverserDFS< TGraph, TIndex, ExactMatching, TStatsSpec >
    : public TraverserBase< TGraph, TIndex, DFS, ExactMatching, TStatsSpec >
    {
      public:
        /* ====================  TYPEDEFS      ======================================= */
        typedef TraverserBase< TGraph, TIndex, DFS, ExactMatching, TStatsSpec > base_type;
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
        TraverserDFS( const graph_type* g, const records_type* r, TIndex* index,
            unsigned int len )
          : base_type( g, r, index, len ), cstate( index, 0, 0, 0, 0 )
        { }

        TraverserDFS( const graph_type* g, unsigned int len )
          : base_type( g, len ), cstate( nullptr, 0, 0, 0, 0 )
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
            stats_type::inc_total_seeds_off_paths( length( saPositions ) );
            for ( i = 0; i < length( saPositions ); ++i )
            {
              output_type hit;
              hit.node_id = cstate.spos.node_id();
              hit.node_offset = cstate.spos.offset();
              hit.read_id = position_to_id( *(this->reads), saPositions[i].i1 );  // Read ID.
              hit.read_offset = position_to_offset( *(this->reads), saPositions[i] );  // Position in the read.
              hit.match_len = this->seed_len;
              hit.gocc = length( saPositions );
              callback( hit );
            }
          }
        }

          inline bool
        compute( )
        {
          if ( cstate.mismatches == 0 ) return false;

          const auto& sequence = this->graph_ptr->node_sequence( cstate.cpos.node_id() );
          assert( cstate.depth < this->seed_len );
          offset_type end_idx = cstate.cpos.offset() + this->seed_len - cstate.depth;
          offset_type i;
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

          if ( !this->graph_ptr->has_edges_out( this->cstate.cpos.node_id() ) ) {
            cstate.mismatches = 0;
            return;
          }
          bool first = true;
          this->graph_ptr->for_each_edges_out(
              this->cstate.cpos.node_id(),
              [this, &first]( id_type to, linktype_type ) {
                if ( first ) {
                  this->cstate.cpos.set_node_id( to );
                  this->cstate.cpos.set_offset( 0 );
                  this->cstate.end = false;
                  first = false;
                  return true;
                }
                this->states.emplace_back( this->cstate );
                this->states.back().cpos.set_node_id( to );
                this->states.back().cpos.set_offset( 0 );
                return true;
              } );
        }
      private:
        /* ====================  DATA MEMBERS  ======================================= */
        typename traits_type::TState cstate;
    };  /* --- end of template class TraverserDFS --- */
}  /* --- end of namespace psi --- */

#endif  /* --- #ifndef PSI_TRAVERSER_DFS_HPP__ --- */
