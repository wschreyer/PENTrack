#!/bin/sh

CGAL_VERSION=CGAL-4.11

echo "I'm going to download and compile ${CGAL_VERSION} for you."
echo "This will take about 100MB of disk space."
echo "Press ENTER to continue, Ctrl+C to cancel."
read -p "" key

mkdir cgal/
cd cgal/
wget -nc -O ${CGAL_VERSION}.tar.xz https://github.com/CGAL/cgal/releases/download/releases%2F${CGAL_VERSION}/${CGAL_VERSION}.tar.xz
tar -xf ${CGAL_VERSION}.tar.xz
rm ${CGAL_VERSION}.tar.xz
cd ..

mkdir cgal/build
cd cgal/build
cmake ../${CGAL_VERSION}/ -DWITH_CGAL_Qt5=OFF -DWITH_CGAL_ImageIO=OFF -DCGAL_HEADER_ONLY=ON
cd ../..
