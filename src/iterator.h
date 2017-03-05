/**
 *    @file  iterator.h
 *   @brief  Generic iterator template class.
 *
 *  Generic iterator template class header file mostly containing class documentations.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Thu Mar 02, 2017  18:00
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef  ITERATOR_H__
#define  ITERATOR_H__

namespace grem {

  /**
   *  @class  Iterator
   *  @brief  Generic iterator template class.
   *
   *  @tparam  TContainer The container.
   *  @tparam  TSpec The type for further specialization of the Iterator type.
   *
   *  The member typedef `Type` will be the proper iterator type of the `TContainer`.
   */
  template<typename TContainer, typename TSpec>
    class Iterator;

}  /* -----  end of namespace grem  ----- */

#endif  /* ----- #ifndef ITERATOR_H__  ----- */
