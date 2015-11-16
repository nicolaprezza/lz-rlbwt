/*
 *  This file is part of lz-rlbwt.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>,
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, and Mathieu Raffinot
 *
 *   lz-rlbwt is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   lz-rlbwt is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * st_subset.h
 *
 *  Created on: Mar 26, 2015
 *      Author: nicola
 *
 *  This class is a more efficient implementation of function I(W,V') described in Lemma 2
 *  in the paper "Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza, and Mathieu Raffinot.
 *  Composite repetition-aware data structures".
 *  This implementation is faster and takes half of the space with respect to the one described
 *  in the paper, and consists basically in the list of ST nodes (represented as BWT intervals) ordered
 *  according to the lexicographic comparator  of the strings they represent. An interval query then takes
 *  just 2 binary searches on the list.
 *
 *  Input: a subset of nodes V'={v_1, ..., v_k} of the Suffix Tree of a text T, represented as
 *  intervals on the BWT of T, plus the BWT length
 *
 *  The structure takes 2k words of space and permits to answer in O(log k) time queries of this kind:
 *
 *  Given the (closed) interval of a word W in the BWT of T, return the (closed) interval of W
 *  among the lexicographically sorted list of node labels l(v_1), ..., l(v_k)
 *
 *	The class is a template on the word type: to save space, choose the smallest word length
 *	compatible with the integers used (by default, 64 bits)
 *
 *	Indexes start from 0!
 *
 */

#ifndef ST_SUBSET_H_
#define ST_SUBSET_H_

#include <packed_view.h>
#include "../internal/sparse_sd_vector.h"

using namespace std;
using namespace bwtil;

namespace lzrlbwt{

class st_subset{

public:

	/* Default empty constructor
	 */
	st_subset(){}

	/* Constructor with bwt length and intervals. Builds the entire structure.
	 * Note that indexes start from 0!
	 *
	 * \param bwt_length	length of the BWT (i.e. 1 + text length)
	 * \param st_nodes		set of suffix tree nodes represented as (closed) intervals on the BWT
	 */
	st_subset(ulint bwt_length, set<pair<ulint,ulint> > st_n){

		this->bwt_length=bwt_length;

		auto st_nodes_vec = vector<pair<ulint,ulint> >(st_n.size());
		std::copy(st_n.begin(),st_n.end(),st_nodes_vec.begin());

		auto comp = [](const pair<ulint,ulint> &A, const pair<ulint,ulint> &B){

			//A precedes B
			bool preceded = A.second < B.first;

			//B is contained in A and the two intervals are distinct
			bool contained = 	(B.first >= A.first and B.second < A.second) or
								(B.second <= A.second and B.first > A.first);

			return preceded or contained;

		};

		std::sort(st_nodes_vec.begin(),st_nodes_vec.end(),comp);

		st_nodes = packed_interval_vector(st_nodes_vec);

	}

	/* Given the BWT interval of a string W, returns the (closed) interval of W
	 * among the lexicographically sorted list of node labels l(v_1), ..., l(v_k)
	 *
	 */
	pair<ulint,ulint> interval(pair<ulint,ulint> bwt_interval){

		//check if the object has been constructed properly
		assert(bwt_length>0);

		//init result with an empty interval
		pair<ulint,ulint> result = {1,0};

		/* this is the lexicographic comparator, i.e. the comparator between intervals that
		 * reflects the total order of the input suffix tree nodes:
		 * two distinct intervals A and B are such that A<B iif:
		 *
		 * 		1) A.second < B.first	(A precedes B)	or
		 * 		2) (B.first >= A.first and B.second < A.second) or	(B.second <= A.second and B.first > A.first)	(B is contained in A)
		 *
		 */
		auto comp1 = [](const pair<ulint,ulint> &A, const pair<ulint,ulint> &B){

			//A precedes B
			bool preceded = A.second < B.first;

			//B is contained in A and the two intervals are distinct
			bool contained = 	(B.first >= A.first and B.second < A.second) or
								(B.second <= A.second and B.first > A.first);

			return preceded or contained;

		};

		auto l = std::lower_bound( st_nodes.begin(),st_nodes.end(), bwt_interval, comp1);

		//all intervals are < bwt_interval
		if(l==st_nodes.end())
			return result;

		//if first interval found is disjoint with bwt_interval, then return empty interval
		if((*l).first > bwt_interval.second)
			return result;

		//first coordinate of the result: first interval that is >= bwt_interval
		result.first = l - st_nodes.begin();

		//l is the last interval : return <l,l>
		if((l+1)==st_nodes.end()){

			result.second = result.first;
			return result;

		}

		//if first interval is exactly bwt_interval, then move forward the iterator to the first interval > bwt_interval
		if(*l == bwt_interval)
			l++;


		//now search from l onwards with a new comparator:
		//A < B iif A precedes B OR A is contained in B
		auto comp2 = [](const pair<ulint,ulint> &A, const pair<ulint,ulint> &B){

			bool precedes = A.second < B.first;
			bool contained = A.first >= B.first and A.second <= B.second;

			return precedes or contained;

		};

		//r: pointer to first interval X > bwt_interval.
		auto r = std::lower_bound( l,st_nodes.end(), bwt_interval, comp2 );

		assert(r>st_nodes.begin());

		//last interval < bwt_interval
		r--;

		result.second = r - st_nodes.begin();

		return result;

	}

	// save the structure to file
	/* \param output_basename	basename of the output file. Extension '.st_subset' will be automatically added to the input path
	 */
	void save_to_file(string filename){

		std::ofstream out (filename,std::ofstream::binary);
		serialize(out);
		out.close();

	}

	// load index from file
	/* \param basename_path		basename of the index files
	 */
	void load_from_file(string filename){

		std::ifstream in (filename,std::ifstream::binary);
		load(in);
		in.close();

	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint w_bytes = 0;

		out.write((char *)&bwt_length, sizeof(bwt_length));

		w_bytes += sizeof(ulint);

		w_bytes += st_nodes.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char *)&bwt_length, sizeof(bwt_length));
		st_nodes.load(in);

	}

private:

	/*
	 * this class implements a vector of BWT intervals where first coordinates
	 * are sorted in increasing order. It uses Elias-Fano for the 1st coordinates, and
	 * a packed vector for the 2nd coordinates
	 */
	class packed_interval_vector{

	   class piv_iterator  :
			    public std::iterator<random_access_iterator_tag,pair<ulint,ulint> >{
			friend class packed_interval_vector;

			packed_interval_vector *_piv = nullptr;
			ulint _index = 0;

			piv_iterator(packed_interval_vector *v, ulint index)
				: _piv(v), _index(index) { }

		public:

            //using iterator_category = std::random_access_iterator_tag;

			piv_iterator() = default;
			piv_iterator(piv_iterator const&) = default;

			piv_iterator &operator=(piv_iterator const&) = default;

			// Iterator
			pair<ulint,ulint> operator*() const {
				return (*_piv)[_index];
			}

			piv_iterator &operator++() {
				++_index;
				return *this;
			}

			// EqualityComparable
			bool operator==(piv_iterator it) const {
				return _index == it._index;
			}

			// ForwardIterator
			bool operator!=(piv_iterator it) const {
				return _index != it._index;
			}

			piv_iterator operator++(int) {
				piv_iterator it(*this);
				++_index;
				return it;
			}

			// BidirectionalIterator
			piv_iterator &operator--() {
				--_index;
				return *this;
			}

			piv_iterator operator--(int) {
				piv_iterator it(*this);
				--_index;
				return it;
			}

			// RandomAccessIterator
			piv_iterator &operator+=(ulint n) {
				_index += n;
				return *this;
			}

			piv_iterator operator+(ulint n) const {
				piv_iterator it(*this);

				it += n;

				return it;
			}

			friend piv_iterator operator+(ulint n, piv_iterator it)
			{
				return it + n;
			}

			piv_iterator &operator-=(ulint n) {
				_index -= n;
				return *this;
			}

			piv_iterator operator-(ulint n) const {
				piv_iterator it(*this);

				it -= n;

				return it;
			}

			friend piv_iterator operator-(ulint n, piv_iterator it) {
				return it - n;
			}

			ulint operator-(piv_iterator it) {
				return ulint(_index) - ulint(it._index);
			}

			pair<ulint,ulint>  operator[](ulint i) const {
				return (*_piv)[_index + i];
			}

			bool operator<(piv_iterator it) const {
				return _index < it._index;
			}

			bool operator<=(piv_iterator it) const {
				return _index <= it._index;
			}

			bool operator>(piv_iterator it) const {
				return _index > it._index;
			}

			bool operator>=(piv_iterator it) const {
				return _index >= it._index;
			}

		};

	public:

		packed_interval_vector(){}

		packed_interval_vector(vector<pair<ulint,ulint> > pairs){

			n = pairs.size();

			{
				vector<ulint> x_vec;
				for(auto p : pairs)
					x_vec.push_back(p.first);

				x = {x_vec};
			}

			ulint max_y=0;
			for(auto p : pairs)
				if(p.second>max_y)
					max_y=p.second;

			bitsize_y = intlog2(max_y);

			y = packed_view<vector>(bitsize_y,n);
			for(ulint i=0;i<n;i++)
				y[i] = pairs[i].second;

		}

		ulint size(){return n;}

		pair<ulint,ulint> operator[](ulint i){

			assert(i<n);

			return {x[i],y[i]};

		}

		ulint serialize(ostream &out){

			ulint w_bytes = 0;

			out.write((char*)&n,sizeof(n));
			w_bytes += sizeof(n);

			if(n==0) return w_bytes;

			ulint y_container_size = y.container().size();
			out.write((char*)&y_container_size,sizeof(y_container_size));
			w_bytes += sizeof(y_container_size);

			out.write((char*)&bitsize_y,sizeof(bitsize_y));
			w_bytes += sizeof(bitsize_y);

			assert(bitsize_y>0);

			w_bytes += x.serialize(out);

			out.write((char*)y.container().data(),y_container_size*sizeof(ulint));
			w_bytes += y_container_size*sizeof(ulint);

			return w_bytes;

		}

		void load(istream &in){

			in.read((char*)&n,sizeof(n));

			if(n==0) return;

			ulint y_container_size;
			in.read((char*)&y_container_size,sizeof(y_container_size));

			in.read((char*)&bitsize_y,sizeof(bitsize_y));

			assert(bitsize_y>0);

			x.load(in);

			y = packed_view<vector>(bitsize_y,n);
			in.read((char*)y.container().data(),y_container_size*sizeof(ulint));

		}

		piv_iterator begin() { return piv_iterator(this, 0); }
		piv_iterator end()   { return piv_iterator(this, n); }

	private:

		class sd_nondecreasing_vec{

		public:

			sd_nondecreasing_vec(){}

			sd_nondecreasing_vec(vector<ulint> v){

				n = v.size();

				if(v.size()==0) return; //nothing to do

				{
					vector<bool> vb;

					ulint last_el = 0;
					for(ulint i=0;i<v.size();++i){

						assert(v[i]>=last_el && "error: vector is not nondecreasing");

						for(ulint j=0;j<v[i]-last_el;++j)
							vb.push_back(false);

						vb.push_back(true);

						last_el = v[i];

					}

					sdv = sparse_sd_vector<>(vb);

				}

			}

			ulint size(){return n;}

			ulint operator[](ulint i){

				assert(i<n);

				//position of i-th one
				auto pos = sdv.select(i);

				return pos - sdv.rank(pos);

			}

			ulint serialize(ostream& out){

				ulint w_bytes = 0;

				out.write((char*)&n,sizeof(n));
				w_bytes += sizeof(n);

				w_bytes += sdv.serialize(out);

				return w_bytes;

			}

			void load(istream& in){

				in.read((char*)&n,sizeof(n));
				sdv.load(in);

			}

		private:

			ulint n=0;
			sparse_sd_vector<> sdv;

		};

		sd_nondecreasing_vec x;
		packed_view<vector> y;

		ulint bitsize_y = 0;
		ulint n = 0;

	};

	ulint bwt_length=0;

	packed_interval_vector st_nodes;

};

}

#endif /* ST_SUBSET_H_ */
