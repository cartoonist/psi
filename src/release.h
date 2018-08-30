/**
 *    @file  release.h
 *   @brief  Release information.
 *
 *  This header file contains release information, such as version, short and long
 *  description, and other constant values.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 *
 *  @internal
 *       Created:  Sun Nov 13, 2016  00:08
 *  Organization:  Max-Planck-Institut fuer Informatik
 *     Copyright:  Copyright (c) 2017, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef RELEASE_H__
#define RELEASE_H__

namespace grem
{
  const char * const VERSION = "0.0.2";                /**< @brief Version number. */
  const char * const PACKAGE = "grem";                 /**< @brief Package name. */
  const char * const SHORT_DESC = "Graph REad Mapper"; /**< @brief Short description. */
  /** @brief Long description. */
  const char * const LONG_DESC = "Map DNA reads to a reference graph.";
  /** @brief Banner to be printed on help and version messages. */
  const char * const BANNER =
    "\n"
    " ┌─┐┬─┐┌─┐┌┬┐\n"
    " │ ┬├┬┘├┤ │││\n"
    " └─┘┴└─└─┘┴ ┴\n";
}

#endif  // RELEASE_H__
