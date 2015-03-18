#!/bin/bash

cd ../..
./PENTrack 0 test/FieldTest/ test/FieldTest
cd test/FieldTest
root -l -q -c '../../out/TREEendlog.c("BFCut.out")'
rm BFCut.out
root -l BFCut.out.root -c showfieldmap.cxx
