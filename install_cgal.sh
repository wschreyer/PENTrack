#!/bin/sh

echo "I'm going to download and compile CGAL for you. This will take up to 1GB of disk space. Continue?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;;
    esac
done

git submodule init
git submodule update
mkdir cgal/build
cd cgal/build
cmake ../src
make -j8
cd ../..
