/*
 * test_small.cpp
 *
 *  Created on: 21 Aug 2019
 *      Author: Timothy Spain, t.spain@ucl.ac.uk
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <memory>

#include <CCfits/CCfits>

#include "../src/glimpse.h"

#ifdef DEBUG_FITS
  bool dbf = true;
#else
  bool dbf = false;
#endif

void read_fits_result(char* filename, std::valarray<unsigned long>& contents);

void small_test() {

    std::cout << "DEBUG_FITS = " << dbf << std::endl;
}

TEST_CASE( "Test verbosity", "[verbosity]") {
    small_test();
    CHECK(!dbf);
}

TEST_CASE( "Reduced resolution", "[example]") {
    int argc = 6;

    char* out_file = "delta_test.fits";
    char* ref_file = "data/delta_ref.fits";
    char* argv[] = {"glimpse",
            "-g",
            "0",
            "data/config3d.ini",
            "data/cat_3_1.fits",
            out_file};

    boost::property_tree::ptree pt;
    po::variables_map vm;

    REQUIRE(create_config(argc, argv, pt, vm) == config_ok_go);

    REQUIRE(configure_and_run(pt, vm) == return_ok);

    // Read the output and reference FITS file
    std::valarray<unsigned long> output_contents;
    read_fits_result(out_file, output_contents);
    std::valarray<unsigned long> reference_contents;
    read_fits_result(ref_file, reference_contents);

    REQUIRE(reference_contents[0] == output_contents[0]);
}

void read_fits_result(char* filename, std::valarray<unsigned long>& contents) {
    std::unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(filename, CCfits::Read, true));
    CCfits::PHDU& image = pInfile->pHDU();
    image.readAllKeys();
    image.read(contents);
    pInfile->destroy();

}
