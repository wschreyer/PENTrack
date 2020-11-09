#!/bin/bash

runTest () {
  cd ../..
  sed s/WPFREQ/${1}/g test/Ramsey/Ramsey_up.in > test/Ramsey/${1}.in
  nice -n 19 ./PENTrack $1 test/Ramsey/${1}.in test/Ramsey/
  rm test/Ramsey/${1}.in
  cd -
}

for f in {183200..183300}; do
  runTest $f #&
done
