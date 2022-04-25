#download hopscotch_map
cd extern/DYNAMIC
rm -rf  build
cd build
cmake ..
make
cd ../..

rm -rf build
mkdir build 
cd build
cmake ..
make
