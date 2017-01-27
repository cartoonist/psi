/*
 * =====================================================================================
 *
 * Filename: release.h
 *
 * Created: Sun Nov 13, 2016  00:08
 * Last modified: Fri Jan 27, 2017  02:43
 *
 * Description: Release information.
 *
 * Copyright (c) 2016, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari, <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#ifndef RELEASE_H__
#define RELEASE_H__

namespace grem
{
  const char * const VERSION = "0.0.1.dev";
  const char * const PACKAGE = "grem";
  const char * const SHORT_DESC = "Graph REad Mapper";
  const char * const LONG_DESC = "Map DNA reads to a reference graph.";
  const char * const BANNER =
    "\n"
    " ┌─┐┬─┐┌─┐┌┬┐\n"
    " │ ┬├┬┘├┤ │││\n"
    " └─┘┴└─└─┘┴ ┴\n";
}

#endif  // RELEASE_H__
