/*
 * =====================================================================================
 *
 * Filename: vargraph_iter.h
 *
 * Created: Wed Jan 18, 2017  14:40
 * Last modified: Fri Jan 20, 2017  01:16
 *
 * Description: Variation graph iterator for traversing the graph.
 *
 * Copyright (c) 2017, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari (cartoonist), ali.ghaffaari@mpi-inf.mpg.de
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef VARGRAPH_ITER_H__
#define VARGRAPH_ITER_H__

#include <unordered_set>
#include <deque>

#include "vargraph.h"

namespace grem
{
  template<typename TGraph, typename TSpec>
    class Iterator;
  template<typename TGraph, typename TSpec>
    bool at_end(Iterator<TGraph, TSpec> it);
  template<typename TGraph, typename TSpec>
    Iterator<TGraph, TSpec> begin(const TGraph & g);

  template<typename TGraph, typename TSpec>
    class Iterator
    {
      friend bool at_end<TGraph, TSpec>(Iterator<TGraph, TSpec> it);
      friend Iterator<TGraph, TSpec> begin<TGraph, TSpec>(const TGraph & g);

      public:
        Iterator(const TGraph & vargraph) : Iterator(&vargraph) {}
        Iterator(const TGraph * _vargraph_ptr)
        {
          *this = begin<TGraph, TSpec>(*_vargraph_ptr);  // use move assign.
        }

        typename TSpec::Value operator*() { return this->itr_value; }
        typename TSpec::Value level() { return this->itr_level; }
        Iterator & operator++();
      private:
        const TGraph * vargraph_ptr;
        typename TSpec::Value itr_value;  // e.g. current node id pointed by iterator.
        typename TSpec::Value itr_level;
        std::unordered_set<typename TSpec::Value> seen;
        typename TSpec::TContainer visiting_buffer;
        typename TSpec::TContainer visiting_levels;

        Iterator() : vargraph_ptr(nullptr) {}
    };

  template<typename TSpec = void>
    struct BfsIterator
    {
      typedef unsigned long int Value;
      typedef std::deque<Value> TContainer;
    };

  template<>
    Iterator<VarGraph, BfsIterator<>> &
    Iterator<VarGraph, BfsIterator<>>::operator++()
    {
      if (this->vargraph_ptr->has_fwd_edge(this->itr_value))
      {
        auto edges = this->vargraph_ptr->fwd_edges(this->itr_value);
        for (auto it = edges.begin(); it != edges.end(); ++it)
        {
          auto adj_node = (*it)->to();
          if (seen.find(adj_node) == seen.end())
          {
            this->visiting_buffer.push_back(adj_node);
            this->visiting_levels.push_back(this->itr_level + 1);
            this->seen.insert(adj_node);
          }
        }
      }
      this->visiting_buffer.pop_front();
      this->visiting_levels.pop_front();
      this->itr_value = this->visiting_buffer.front();
      this->itr_level = this->visiting_levels.front();

      return *this;
    }

  template<>
  bool at_end(Iterator<VarGraph, BfsIterator<>> it)
  {
    return it.visiting_buffer.empty();
  }

  template<>
    Iterator<VarGraph, BfsIterator<>>
    begin(const VarGraph & g)
    {
      Iterator<VarGraph, BfsIterator<>> begin_it;
      begin_it.vargraph_ptr = &g;
      begin_it.seen.reserve(begin_it.vargraph_ptr->nodes_size());
      begin_it.visiting_buffer.push_back(g.node_at(0).id());
      begin_it.visiting_levels.push_back(0);
      begin_it.itr_value = begin_it.visiting_buffer.front();
      begin_it.itr_level = begin_it.visiting_levels.front();
      begin_it.seen.insert(begin_it.itr_value);

      return begin_it;
    }
}  // namespace grem

#endif  // VARGRAPH_ITER_H__
