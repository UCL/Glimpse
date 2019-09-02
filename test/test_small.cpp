/*
 * test_small.cpp
 *
 *  Created on: 21 Aug 2019
 *      Author: Timothy Spain, t.spain@ucl.ac.uk
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>

#include "../src/glimpse.h"

#ifdef DEBUG_FITS
  bool dbf = true;
#else
  bool dbf = false;
#endif

  void small_test() {

  std::cout << "DEBUG_FITS = " << dbf << std::endl;
}

TEST_CASE( "Test verbosity", "[verbosity]") {
    small_test();
    CHECK(!dbf);
}

TEST_CASE( "Reduced resolution", "[example]") {
    int argc = 6;
    char* argv[] = {"glimpse",
            "-g",
            "0",
            "../example/reduced_resolution/config3d.ini",
            "../example/reduced_resolution/cat_3_1.fits",
            "delta.fits"};

    boost::property_tree::ptree pt;
    po::variables_map vm;

    REQUIRE(create_config(argc, argv, pt, vm) == 0);


}
