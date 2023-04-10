#!/bin/bash

cd ../..
./PENTrack 0 test/HitTest/config.in test/HitTest
cd test/HitTest
root -l 000000000000.root showscatterdist.cxx
