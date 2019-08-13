#!/bin/bash

cd ../..
./PENTrack 0 test/HarmonicFieldTest/config.in test/HarmonicFieldTest/
cd test/HarmonicFieldTest
root -l -q -c '../../out/TREEendlog.c("BFCut.out")'
rm BFCut.out
root -l BFCut.out.root -c showfieldmap.cxx
