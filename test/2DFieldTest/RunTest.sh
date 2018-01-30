#!/bin/bash

cd ../..
./PENTrack 0 test/2DFieldTest/config.in test/2DFieldTest/
cd test/2DFieldTest
root -l -q -c '../../out/TREEendlog.c("BFCut.out")'
rm BFCut.out
root -l BFCut.out.root -c showfieldmap.cxx
