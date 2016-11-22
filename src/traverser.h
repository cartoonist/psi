/*
 * =====================================================================================
 *
 * Filename: traverser.h
 *
 * Created: Mon Nov 14, 2016  01:11
 * Last modified: Wed Nov 23, 2016  00:10
 *
 * Description: Traversers class definitions.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef TRAVERSER_H__
#define TRAVERSER_H__

#include <vector>
#include <functional>

#include "vg.pb.h"
#include "vargraph.h"
#include "types.h"

// TODO: refactor: types (const, * and &).

namespace grem
{
  class PathTraverser
  {
    public:
      // Member typedefs and classes
      typedef vg::Alignment Output;
      // Traverse parameters
      class Param
      {
        friend class PathTraverser;

        public:
          // Constructors
          Param(const ReadsChunk &reads_, unsigned int seed_len_) :
            reads(reads_), reads_index(reads.seqs), seed_len(seed_len_)
          {}

          // Attributes getters and setters
          inline const ReadsChunk     &get_reads()
          { return this->reads; }

          inline const DnaSeqSetIndex &get_reads_index()
          { return this->reads_index; }

          inline unsigned int          get_seed_len()
          { return this->seed_len; }

        private:
          ReadsChunk     reads;
          DnaSeqSetIndex reads_index;
          unsigned int   seed_len;
      };

      // Constructors
      PathTraverser(const VarGraph *graph, PathTraverser::Param *trav_params, vg::Position start);
      PathTraverser(const VarGraph &graph, PathTraverser::Param &trav_params, vg::Position start);
      PathTraverser(const PathTraverser &other, vg::Position new_locus);
      // TODO: Move/copy constructors.
      // TODO: Move/copy assignment operators.

      // Traverse interface functions (are friends!)
      friend void move_forward(PathTraverser &ptrav, std::vector< PathTraverser > &new_ptravs);
      friend bool is_finished(PathTraverser &ptrav);
      friend bool is_valid(PathTraverser &ptrav);
      friend void get_results(PathTraverser &ptrav, std::vector< PathTraverser::Output > &results);

      // Attributes getters and setters
      inline const VarGraph *             get_vargraph()
      { return this->vargraph; }

      inline const PathTraverser::Param * get_paramters()
      { return this->parameters; }

      inline vg::Position                 get_s_locus()
      { return this->s_locus; }

      inline vg::Position                 get_c_locus()
      { return this->c_locus; }

      inline const vg::Path &             get_path()
      { return this->path; }

    private:
      // Internal typedefs and classes
      typedef struct {
        DnaSSIndexIter iter;
        unsigned int   boffset;
      } IterState;

      // Attributes
      const VarGraph *         vargraph;        // pointer to variation graph.
      PathTraverser::Param *   parameters;      // pointer to params (shared between traversers).
      vg::Position             s_locus;         // starting locus
      vg::Position             c_locus;         // current locus
      std::vector< IterState > iters_state;
      vg::Path                 path;
      bool                     finished;

      // Internal methods
      bool is_seed_hit();
      bool go_down(IterState &its, seqan::Value<DnaSeq>::Type c);
      void go_down_all(seqan::Value<DnaSeq>::Type c);
      void extend_path(unsigned int visit_len);
      void one_node_forward();
      void get_results(std::vector< PathTraverser::Output > &results);
  };

  template <class TPathTraverser>
    class GraphTraverser
    {
      public:
        // Constructors
        GraphTraverser(const VarGraph *graph,
                       const std::vector< vg::Position > *start_loci);
        GraphTraverser(const VarGraph &graph,
                       const std::vector< vg::Position > *start_loci);
        GraphTraverser(const VarGraph *graph);
        GraphTraverser(const VarGraph &graph);
        // TODO: Move/copy constructors.
        // TODO: Move/copy assignment operators.
        // Public methods
        inline void add_start(vg::Position &locus)
        { this->starting_points.push_back(locus); }
        void traverse(typename TPathTraverser::Param trav_params,
                      std::function< void(typename TPathTraverser::Output &) > callback);
      private:
        // Attributes
        const VarGraph *             vargraph;
        std::vector< vg::Position >  starting_points;
    };
}

#endif
