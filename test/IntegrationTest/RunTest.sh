#!/bin/bash

cd ../..
./PENTrack 0 test/IntegrationTest/config.in test/IntegrationTest
cd test/IntegrationTest
root -l -q -c ../../out/merge_all.c
rm 000000000000neutronend.out 000000000000neutronhit.out
root -l out.root -c showintegrationresult.cxx
