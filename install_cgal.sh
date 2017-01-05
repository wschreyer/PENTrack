#!/bin/sh

echo "I'm going to download and compile CGAL for you."
echo "This will take about 100MB of disk space."
echo "Press ENTER to continue, Ctrl+C to cancel."
read -p "" key

mkdir cgal/
cd cgal/
wget -nc https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.9/CGAL-4.9.tar.xz
tar -xf CGAL-4.9.tar.xz
rm CGAL-4.9.tar.xz
cd ..

mkdir cgal/build
cd cgal/build
cmake ../CGAL-4.9/ -DWITH_CGAL_Qt5=OFF -DWITH_CGAL_ImageIO=OFF -DCGAL_HEADER_ONLY=ON
cd ../..
