#!/bin/bash

cd ../..
./PENTrack 0 test/comsolFieldTest/config.in test/comsolFieldTest
cd test/comsolFieldTest
python3 graphstuff.py
