lz-rlbwt: Run-Length Compressed Burrows-Wheeler transform with LZ77 suffix array sampling
===============
Author: Nicola Prezza (nicolapr@gmail.com)

From the paper: "Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza, and Mathieu Raffinot. Composite repetition-aware data structures"

### Brief description

The lz-rlbwt data structure is a fast index that takes advantage of text repetitions in order to reduce its memory footprint. Highly repetitive inputs produce longer BWT runs (RLE compressed using an Elias-Fano-encoded RLBWT data structure) and a smaller LZ77 parse. The lz-rlbwt exploits this fact and substitutes SA samples with components from a LZ77 index (basically, a sampled suffix tree and geometric range search data structures). 

The user can choose among 3 implementations:

1. Full index: forward RLBWT + reverse RLBWT + 2-sided and 4-sided range search. Linear-time search: O(m log n + occ log n)
2. Bidirectional index: forward RLBWT + 2-sided and 4-sided range search (reverse RLBWT is simulated using the forward RLBWT). Quadratic-time search: O(m^2 log n + occ log n)
3. Light index (default): reverse RLBWT + 2-sided range search. Linear-time search: O((m+1) occ log n)

### Download

To clone the repository, use the --recursive option:

> git clone --recursive http://github.com/nicolaprezza/lz-rlbwt

### Compile

The library has been tested under linux using gcc 4.9.2. You need the SDSL library installed on your system (https://github.com/simongog/sdsl-lite).

lz-rlbwt uses cmake to generate the Makefile. Create a build folder in the main lz-rlbwt folder:

> mkdir build

execute cmake:

> cd build; cmake ..

and compile:

> make

### Run

After compiling, run 

>  lz-rlbwt-build input.txt

This command will create the LZ-RLCSA index of the text file "input.txt" and will store it using as prefix "input.txt". Use option -o to specify a different basename for the index files. At the moment, building the index takes a lot of time since we use a slow but memory-efficient LZ77 parser.

Run

> lz-rlbwt-search index_basename patterns_file

to align a set of patterns in pizza&chilli format (http://pizzachili.dcc.uchile.cl/experiments.html) using the lz-rlcsa index. At the moment, this executable does not store the output anywhere (we use it only for benchmarking purposes)
