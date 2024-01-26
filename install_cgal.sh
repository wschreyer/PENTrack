#!/bin/sh

CGAL_VERSION=5.6

echo "I'm going to download and compile ${CGAL_VERSION} for you."
echo "This will take about 50MB of disk space."
echo "Press ENTER to continue, Ctrl+C to cancel."
read -p "" key

mkdir cgal/
wget -nc https://github.com/CGAL/cgal/releases/download/v${CGAL_VERSION}/CGAL-${CGAL_VERSION}-library.tar.xz
tar -xf CGAL-${CGAL_VERSION}-library.tar.xz
rm CGAL-${CGAL_VERSION}-library.tar.xz
mv CGAL-${CGAL_VERSION}/ cgal
