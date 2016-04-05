#!/bin/bash

pushd .
cd deps
git clone --depth 1 --singlebranch -b master https://github.com/OPM/opm-data
cd opm-data/spe1
$WORKSPACE/serial/build-opm-autodiff/bin/flow deck_filename=SPE1CASE2.DATA
cd ..
cd spe3
$WORKSPACE/serial/build-opm-autodiff/bin/flow max_iter=50 deck_filename=SPE3CASE1.DATA
cd ..
cd spe9
$WORKSPACE/serial/build-opm-autodiff/bin/flow max_iter=50 deck_filename=SPE9_CP.DATA

# Compare OPM with eclipse reference
PYTHONPATH=/usr/local/python python output_comparator/src/compare_eclipse.py spe1/eclipse-simulation/ spe1/ SPE1CASE2 0.01 0.01
PYTHONPATH=/usr/local/python python output_comparator/src/compare_eclipse.py spe3/eclipse-simulation/ spe3/ SPE3CASE1 0.02 0.02
PYTHONPATH=/usr/local/python python output_comparator/src/compare_eclipse.py spe9/eclipse-simulation/ spe9/ SPE9_CP 0.002 0.001

# Compare OPM with OPM reference
PYTHONPATH=/usr/local/python python output_comparator/src/compare_eclipse.py spe1/opm-simulation-reference/ spe1/ SPE1CASE2 0.001 0.001
PYTHONPATH=/usr/local/python python output_comparator/src/compare_eclipse.py spe3/opm-simulation-reference/ spe3/ SPE3CASE1 0.001 0.001
PYTHONPATH=/usr/local/python python output_comparator/src/compare_eclipse.py spe3/opm-simulation-reference/ spe9/ SPE9_CP 0.002 0.007
