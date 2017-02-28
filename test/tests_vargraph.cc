/*
 * =====================================================================================
 *
 * Filename: tests-vargraph.cc
 *
 * Created: Mon Feb 13, 2017  22:24
 * Last modified: Tue Feb 28, 2017  10:27
 *
 * Description: Test VarGraph class.
 *
 * Copyright (c) 2017, Ali Ghaffaari
 *
 * Author: Ali Ghaffaari (cartoonist), <ali.ghaffaari@mpi-inf.mpg.de>
 * Organization: Max-Planck-Institut fuer Informatik
 *
 * =====================================================================================
 */

#include <random>
#include <chrono>
#include <libgen.h>
#include <string>

#include <catch/catch.hpp>

#include "tests_base.h"
#include "vargraph.h"
#include "logger.h"

INITIALIZE_EASYLOGGINGPP

using namespace grem;


// Test if the parameter is a constant reference.
template<typename T>
bool is_const_ref(const T&) { return true; }

// Test if the parameter is a non-const pointer.
template<typename T>
bool is_mutable_ptr(T*) { return true; }

std::string testdir(TESTDIR);


SCENARIO( "Loading variation graph from a vg file", "[input]" )
{
  GIVEN( "A small graph without graph name" )
  {
    std::string vgpath = testdir + "/data/small/x.vg";
    VarGraph vargraph(vgpath.c_str());

    REQUIRE( vargraph.nodes_size() == 210 );

    grem::VarGraph::NodeID node_id = 0;
    REQUIRE( !vargraph.has_node(node_id) );

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<unsigned int> distribution(0, vargraph.nodes_size()-1);
    unsigned int random_idx = distribution(generator);

    auto& a_node = vargraph.node_at(random_idx);
    REQUIRE( vargraph.has_node(a_node.id()) );
    REQUIRE( vargraph.has_node(&a_node) );
    REQUIRE( vargraph.node_at(random_idx).id() == random_idx+1 );
    REQUIRE( vargraph.node_by(a_node.id()).sequence() == a_node.sequence() );

    REQUIRE( is_const_ref(vargraph.node_at(random_idx)) );
    REQUIRE( is_mutable_ptr(vargraph.mutable_node_at(random_idx)) );
    REQUIRE( is_const_ref(vargraph.node_by(a_node.id())) );
    REQUIRE( is_mutable_ptr(vargraph.mutable_node_by(a_node.id())) );

    REQUIRE( vargraph.edges_size() == 291 );

    random_idx = distribution(generator);

    auto& an_edge = vargraph.edge_at(random_idx);
  }
}
