#!/bin/bash

cd ../..
./PENTrack 0 test/IntegrationTest/config.in test/IntegrationTest
cd test/IntegrationTest
root -l 000000000000.root -c showintegrationresult.cxx
