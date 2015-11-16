/*
 *  This file is part of BWTIL.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>
 *
 *   BWTIL is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   BWTIL is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

//============================================================================
// Name        : sparse_vector.h
// Author      : Nicola Prezza
// Version     : 1.0
// Description :

 /*
  * sparse_sd_vector: an wrapper on sd_vector of the sdsl library
  */

//============================================================================


#ifndef INTERNAL_SPARSE_SD_VECTOR_H_
#define INTERNAL_SPARSE_SD_VECTOR_H_

#include <vector>

using namespace std;
using namespace sdsl;

#ifndef ulint
typedef uint64_t ulint;
#endif

#ifndef uint
typedef uint32_t uint;
#endif

namespace lzrlbwt{

template<typename = ulint> //this is not used: legacy template parameter
class sparse_sd_vector{

public:

	/*
	 * empty constructor. Initialize bitvector with length 0.
	 */
	sparse_sd_vector(){}

	/*
	 * constructor. build bitvector given a vector of bools
	 */
	sparse_sd_vector(vector<bool> &b){

		if(b.size()==0) return;

		u = b.size();

		bit_vector bv(b.size());

		for(uint i=0;i<b.size();++i)
			bv[i] = b[i];

		sdv = sd_vector<>(bv);
		rank1 = sd_vector<>::rank_1_type(&sdv);
		select1 = sd_vector<>::select_1_type(&sdv);

	}

	sparse_sd_vector & operator= (const sparse_sd_vector & other) {

		u = other.sdv.size();
		sdv = sd_vector<>(other.sdv);
		rank1 = sd_vector<>::rank_1_type(&sdv);
		select1 = sd_vector<>::select_1_type(&sdv);

	    return *this;
	}

	/*
	 * not implemented
	 * argument: a boolean b
	 * behavior: append b at the end of the bitvector.
	 */
	//void push_back(bool b){}

	/*
	 * argument: position i in the bitvector
	 * returns: bit in position i
	 * only access! the bitvector is static.
	 */
	bool operator[](ulint i){

		assert(i<size());
		return sdv[i];

	}

	bool at(ulint i){
		return operator[](i);
	}

	/*
	 * argument: position i in the bitvector
	 * returns: number of bits equal to 1 before position i excluded
	 */
	ulint rank(ulint i){

		assert(i<=size());
		return rank1(i);

	}

	/*
	 * retrieve length of the i-th gap (i>=0). gap length includes the leading 1
	 * \param i<number_of_1()
	 *
	 */
	ulint gapAt(ulint i){

		assert(i<number_of_1());

		if(i==0) return select(0)+1;

		return select(i)-select(i-1);

	}

	/*
	 * argument: integer i>0
	 * returns: position of the i-th one in the bitvector. i starts from 0!
	 */
	ulint select(ulint i){

		assert(i<number_of_1());
		return select1(i+1);//in sd_vector, i starts from 1

	}

	/*
	* returns: size of the bitvector
	*/
	ulint size(){return u;}

	/*
	 * returns: number of 1s in the bitvector
	 */
	ulint number_of_1(){return rank1(size()); }

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint w_bytes = 0;

		out.write((char*)&u, sizeof(u));

		w_bytes += sizeof(u);

		if(u==0) return w_bytes;

		w_bytes += sdv.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&u, sizeof(u));

		if(u==0) return;

		sdv.load(in);
		rank1 = sd_vector<>::rank_1_type(&sdv);
		select1 = sd_vector<>::select_1_type(&sdv);

	}

private:

	//bitvector length
	ulint u = 0;

	sd_vector<> sdv;
	sd_vector<>::rank_1_type rank1;
	sd_vector<>::select_1_type select1;

};

}


#endif /* INTERNAL_SPARSE_SD_VECTOR_H_ */
