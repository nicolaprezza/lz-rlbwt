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
 * lz_rlbwt.h
 *
 *  Created on: Mar 18, 2015
 *      Author: nicola
 *
 *  The lz-rlbwt data strucures combines the rlbwt and components from a LZ77 index. Space is
 *  O(r+z) words, r being the number of runs in the BWT of the text and z being the number of phrases
 *  in the LZ77 parse of the text.
 *
 *  This class permits to answer count and locate queries. Letting n be the text length,
 *  m the pattern length, and occ the number of pattern occurrences, count takes O(m log log n) time,
 *  and locate takes (depending on the version - see constructor):
 *
 *  -full: O( (m+occ) * log n) time
 *  -bidirectional: O( (m^2+occ) * log n ) time
 *  -light: O( m * log n * (occ+1) ) time
 *
 *	bidirectional and light indexes have several optimizations implemented, so the worst cases
 *	rarely should occur.
 *
 *  The class is a template on the word type and on its building blocks (FM index, range search,
 *  sparse bitvector.)
 *  Choosing a smaller ulint (default=32 bit) reduces the size of the structure. Even though it is possible
 *  to specify a 64-bit integer word type, using texts larger than 4GB may cause problems with the underlying
 *  rlbwt data structure (which uses 32-bits integers)
 *
 *  For more info, read the paper: "Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza, and Mathieu
 *  Raffinot. Composite repetition-aware data structures"
 *
 */

#ifndef LZrlbwt_H_
#define LZrlbwt_H_

#include <cstdlib>
#include <iostream>
#include <vector>

#include <sdsl/wavelet_trees.hpp>
#include <rlbwt.h>
#include <data_structures/lz77_parser.h>
#include <sparse_sd_vector.h>
#include <range_search_2D.h>
#include <st_subset.h>
#include <packed_view.h>
#include <dynamic.hpp>

using namespace std;
using namespace sdsl;
using namespace bwtil;
using namespace dyn;

namespace lzrlbwt{

typedef lz77_parser<gap_bv,packed_spsi> lz77_parser_t;
//typedef lz77_parser<> lz77_parser_t;

//type of the run-length BWT: Elias-Fano RLBWT / Huffman RLBWT
enum bwt_type {ef_rlbwt,h_rlbwt};

/*
 * Components stored in the index:
 *
 * - bidirectional: fwd BWT, 4-sided range, 2-sided range, st_subset, bitvector last, SA samples
 * - light: rev BWT, 2-sided range, bitvector last, SA samples
 * - full: fwd BWT, rev BWT, 4-sided range, 2-sided range, st_subset, bitvector last, SA samples
 *
 */
enum lz_rlbwt_type {bidirectional,light,full};

enum lzrlbwt_options{

	verbose_out,		//print verbose ouput
	bidirectional_bwt,	//if true, reverse index is not created / saved / loaded
	light_index,		//if true, both fwd index and 4-sided structure are not
						//created/saved/loaded

};

template<
	class fm_index_t = rlbwt_sd, 		// RLBWT using Elias-Fano compression on the gap lengths.
	class sparse_bitvector_t = sparse_sd_vector<>,		// elias-fano sparse bitvector
	class range_search_t = range_search_2D<> 	// class implementing 2-sided and 4-sided range search.
>
class lz_rlbwt{

public:

	~lz_rlbwt(){

		if(fm_index_fwd != NULL)
			delete fm_index_fwd;

		if(fm_index_rev != NULL)
			delete fm_index_rev;

		if(lz_4_sided_range_search != NULL)
			delete lz_4_sided_range_search;

		if(lz_2_sided_range_search != NULL)
			delete lz_2_sided_range_search;

		if(st_sub != NULL)
			delete st_sub;

		if(last != NULL)
			delete last;

	}

	//default empty constructor
	lz_rlbwt(){};

	// build the lz_rlbwt of the input string
	/* \param input		string to be indexed
	 * \param options	list of options: by default, empty (no verbose and no bidirectional bwt.)
	 * NOTE: we assume that all letters in the input text are right-maximal.
	 */
	lz_rlbwt(string &input, vector<lzrlbwt_options> opt = {}){

		bool verbose=false;
		type = full;

		for(auto o:opt){

			if(o==verbose_out)
				verbose=true;

			if(o==bidirectional_bwt)
				type = bidirectional;

			if(o==light_index)
				type = light;

		}

		build(input, verbose);

	}

	// return number of occurrences of input pattern
	/* \param pattern	string to be searched
	 * \return number of occurrences of the input pattern
	 */
	ulint count(string &pattern){

		ulint occ = 0;
		range_t range;

		//if type == light, we have only reverse index
		if(type==light){

			range = fm_index_rev->count(pattern,true);

		}else{

			range = fm_index_fwd->count(pattern);

		}

		if(range.second >= range.first)
			occ = (range.second - range.first) + 1;

		return occ;

	}

	// locate all occurrences of a pattern in the indexed text
	/* \param pattern	string to be searched
	 * \return all occurrences (i.e. text positions) of the input pattern
	 */
	vector<ulint> locate(string &pattern, bool optimize = true){

		//if index is light, then locate is particular
		if(type == light){

			return locate_light(pattern, optimize);

		}

		//else: locate with full or bidirectional index
		return locate_bidir_full(pattern, optimize);

	}

	// save the index to file
	/* \param 	basename	basename of the output index files. Extensions will be automatically
	 * 			added to the input path
	 * \param verbose	verbose output?
	 */
	void save_to_file(string base_name, bool verbose = false){

		string type_basename = string(base_name).append(".lzrlbwt.index_type");
		string fm_index_fwd_basename = string(base_name).append(".lzrlbwt.fmi_fwd");
		string fm_index_rev_basename = string(base_name).append(".lzrlbwt.fmi_rev");
		string two_sided_filename = string(base_name).append(".lzrlbwt.2_sided_range");
		string four_sided_filename = string(base_name).append(".lzrlbwt.4_sided_range");
		string st_subset_filename = string(base_name).append(".lzrlbwt.st_subset");
		string last_filename = string(base_name).append(".lzrlbwt.end_phrases");
		string start_filename = string(base_name).append(".lzrlbwt.beg_phrases");
		string sa_samples_filename = string(base_name).append(".lzrlbwt.sa_samples");

		//write lz-rlbwt type
		std::ofstream type_ofs (type_basename,std::ofstream::binary);
		type_ofs.write((char*)&type,sizeof(type));
		type_ofs.close();


		{

			if(verbose)
				cout << "Saving start of phrases (text positions) ... " << flush;

			std::ofstream out (start_filename,std::ofstream::binary);
			begin_of_phrase.serialize(out);

			if(verbose)
				cout << "done!" << endl;

		}

		if(type != light){

			if(verbose)
				cout << "Saving forward RLBWT ... " << flush;

			assert(fm_index_fwd!=0);

			fm_index_fwd->save_to_disk(fm_index_fwd_basename);

			if(verbose)
				cout << "done!" << endl;

		}

		if(type!=bidirectional){

			if(verbose)
				cout << "Saving reverse RLBWT ... " << flush;

			assert(fm_index_rev!=0);

			fm_index_rev->save_to_disk(fm_index_rev_basename);

			if(verbose)
				cout << "done!" << endl;

		}

		if(verbose)
			cout << "Saving two-sided range data structure ... " << flush;

		lz_2_sided_range_search->save_to_file(two_sided_filename);

		if(verbose)
			cout << "done!" << endl;

		if(type != light){

			if(verbose)
				cout << "Saving four-sided range data structure ... " << flush;

			lz_4_sided_range_search->save_to_file(four_sided_filename);

			if(verbose)
				cout << "done!" << endl;

			if(verbose)
				cout << "Saving interval data structure on the subset of ST nodes... " << flush;

			st_sub->save_to_file(st_subset_filename);

			if(verbose)
				cout << "done!" << endl;

		}

		if(verbose)
			cout << "Saving end of phrases (BWT positions) ... " << flush;

		std::ofstream out (last_filename,std::ofstream::binary);
		last->serialize(out);

		//save also global variables in this file
		out.write((char*)&terminator_pos,sizeof(ulint));
		out.write((char*)&bwt_length,sizeof(ulint));
		uint16_t a_size = alphabet.size();
		out.write((char*)&a_size,sizeof(uint16_t));
		out.write((char*)alphabet.data(),a_size*sizeof(uchar));

		out.close();

		if(verbose)
			cout << "done!" << endl;

		if(type==light){

			if(verbose)
				cout << "Saving SA samples ... " << flush;

			std::ofstream out (sa_samples_filename,std::ofstream::binary);

			ulint sa_samples_container_size = SA_samples.container().size();
			ulint bitlength = SA_samples.width();

			out.write((char*)&sa_samples_container_size,sizeof(sa_samples_container_size));
			out.write((char*)&bitlength,sizeof(bitlength));

			out.write((char*)SA_samples.container().data(),sa_samples_container_size*sizeof(ulint));
			out.close();

			if(verbose)
				cout << "done!" << endl;

		}

	}

	// index from file
	/* \param basename_path		basename of the index files
	 * \param opt	options: verbose output / load bidirectional index (i.e. do not load reverse bwt)
	 */
	void load_from_file(string base_name, vector<lzrlbwt_options> opt = {}){

		bool verbose=false;

		for(auto o:opt){

			if(o==verbose_out)
				verbose=true;

		}

		string type_basename = string(base_name).append(".lzrlbwt.index_type");
		string fm_index_fwd_basename = string(base_name).append(".lzrlbwt.fmi_fwd");
		string fm_index_rev_basename = string(base_name).append(".lzrlbwt.fmi_rev");
		string two_sided_filename = string(base_name).append(".lzrlbwt.2_sided_range");
		string four_sided_filename = string(base_name).append(".lzrlbwt.4_sided_range");
		string st_subset_filename = string(base_name).append(".lzrlbwt.st_subset");
		string last_filename = string(base_name).append(".lzrlbwt.end_phrases");
		string start_filename = string(base_name).append(".lzrlbwt.beg_phrases");
		string sa_samples_filename = string(base_name).append(".lzrlbwt.sa_samples");

		//read lz-rlbwt type
		std::ifstream type_ofs (type_basename,std::ifstream::binary);
		type_ofs.read((char*)&type,sizeof(type));
		type_ofs.close();

		{

			if(verbose)
				cout << "Loading sparse bitvector data structure (begin of phrases on text) ... " << flush;

			std::ifstream in (start_filename,std::ifstream::binary);
			begin_of_phrase.load(in);

			in.close();

			if(verbose)
				cout << "done!" << endl;

		}

		if(type != light){

			if(verbose)
				cout << "Loading forward FM index ... " << flush;

			fm_index_fwd = new fm_index_t();
			fm_index_fwd->load_from_disk(fm_index_fwd_basename);

			if(verbose)
				cout << "done!" << endl;

		}

		if(type != bidirectional){

			if(verbose)
				cout << "Loading reverse FM index ... " << flush;

			fm_index_rev = new fm_index_t();
			fm_index_rev->load_from_disk(fm_index_rev_basename);

			if(verbose)
				cout << "done!" << endl;

		}

		if(verbose)
			cout << "Loading two-sided range data structure ... " << flush;

		lz_2_sided_range_search = new range_search_t();
		lz_2_sided_range_search->load_from_file(two_sided_filename);

		if(verbose)
			cout << "done!" << endl;

		if(type != light){

			if(verbose)
				cout << "Loading four-sided range data structure ... " << flush;

			lz_4_sided_range_search = new range_search_t();
			lz_4_sided_range_search->load_from_file(four_sided_filename);

			if(verbose)
				cout << "done!" << endl;

			if(verbose)
				cout << "Loading interval data structure on the subset of ST nodes... " << flush;

			st_sub = new st_subset();
			st_sub->load_from_file(st_subset_filename);

			if(verbose)
				cout << "done!" << endl;

		}

		if(verbose)
			cout << "Loading sparse bitvector data structure (end of phrases on BWT) ... " << flush;

		std::ifstream in (last_filename,std::ifstream::binary);
		last = new sparse_bitvector_t();
		last->load(in);

		//load also global variables from this file
		in.read((char*)&terminator_pos,sizeof(ulint));
		in.read((char*)&bwt_length,sizeof(ulint));
		uint16_t a_size;
		in.read((char*)&a_size,sizeof(uint16_t));
		alphabet = vector<uchar>(a_size);
		in.read((char*)alphabet.data(),a_size*sizeof(uchar));

		in.close();

		if(verbose)
			cout << "done!" << endl;

		if(type==light){

			if(verbose)
				cout << "Loading SA samples ... " << flush;

			std::ifstream in (sa_samples_filename,std::ifstream::binary);

			ulint sa_samples_container_size;
			ulint bitlength;

			in.read((char*)&sa_samples_container_size,sizeof(sa_samples_container_size));
			in.read((char*)&bitlength,sizeof(bitlength));

			SA_samples = packed_view<vector>(bitlength,last->rank(last->size()));

			in.read((char*)SA_samples.container().data(),sa_samples_container_size*sizeof(ulint));
			in.close();

			if(verbose)
				cout << "done!" << endl;

		}

	}

private:

	vector<ulint> locate_bidir_full(string &pattern, bool optimize = true){

		//vector containing all occurrences of the pattern
		vector<ulint> occ;

		ulint n_occ = count(pattern);

		if(n_occ==0)
			return occ;

		ulint m = pattern.size();

		//interval of all proper reverse pattern suffixes
		//for convenience, we treat occurrences that start at the
		//beginning/end of a phrase as internal occurrences.
		//length of this vector is m-1
		vector<pair<ulint,ulint> > rev_ranges;

		//range in reverse bwt
		pair<ulint, ulint> rev_range = fm_index_fwd->get_char_range(pattern[0]);

		//if we have the reverse BWT, then we compute now the reverse ranges.
		if(type != bidirectional){

			for(ulint i=1;i<m;++i){

				//first_one = rank of the first one in this range, if any

				assert(rev_range.first<=last->size());
				ulint first_one = last->rank(rev_range.first);
				//ones = number of ones in last bitvector in this interval
				assert(rev_range.second+1<=last->size());
				ulint ones = last->rank(rev_range.second+1) - first_one;

				//if there are phrase ends, then insert a proper range. Otherwise, insert an empty range
				if(ones>0){

					rev_ranges.push_back({first_one,first_one+ones-1});

				}else{

					//insert an empty range
					rev_ranges.push_back({1,0});

				}

				rev_range = fm_index_rev->LF(rev_range,pattern[i]);

			}

		}

		//now, from end to begin of pattern, search fwd ranges in the st_sub structure

		//start from last character
		pair<ulint, ulint> fwd_range = fm_index_fwd->get_char_range(pattern[m-1]);

		for(ulint i=1; i<m and occ.size() < n_occ; ++i){

			pair<ulint, ulint> range_fwd = st_sub->interval(fwd_range);

			//if ranges are non-empty
			//if optimize = false, execute this in any case
			if((not optimize) or (range_fwd.second>=range_fwd.first)){

				if(type == bidirectional){

					//compute rev_range using bidirectional BWT
					rev_range = bidirectional_LF(pattern.substr(0,m-i));

					//first_one = rank of the first one in this range, if any
					assert(rev_range.first<=last->size());
					ulint first_one = last->rank(rev_range.first);
					//ones = number of ones in last bitvector in this interval
					assert(rev_range.second+1<=last->size());
					ulint ones = last->rank(rev_range.second+1) - first_one;

					if(ones>0){

						rev_range = {first_one,first_one+ones-1};

					}else{

						rev_range = {1,0};

					}

				}else{

					//we already computed rev_range
					rev_range = rev_ranges[m-i-1];

				}

				if( rev_range.second>=rev_range.first ){

					//rev ranges start from 1. remap so that they start from 0

					assert(rev_range.first>0);
					assert(rev_range.second>0);

					rev_range.first--;
					rev_range.second--;

					//query 4-sided range data structure
					auto four_sided_result =
							lz_4_sided_range_search->four_sided_range_search(rev_range,range_fwd);

					//cout << "OCC : " << four_sided_result.size() << endl;

					if(four_sided_result.size()>0){

						//store here positions on text found with 4-sided range search
						vector<ulint> occ_temp;

						for(auto o : four_sided_result){

							//o.second is the phrase number (of right factor)

							//convert from phrase rank to text position with a select query

							assert(o.second < begin_of_phrase.rank(begin_of_phrase.size()));
							ulint start_position = begin_of_phrase.select(o.second);

							//compute correct coordinate (left-shift)
							assert(start_position >= m-i);
							occ_temp.push_back( start_position - (m-i) );

						}

						for(auto o : occ_temp)
							occ.push_back(o);

						//call recursive 2-range search
						for(auto o : occ_temp){

							//perform 2-sided range search only if there still are occurrences to be found
							//if optimize = false, do this in any case
							if((not optimize) or occ.size() < n_occ){

								//recursively compute occurrences of phrases that copy T[o,...,o+m-1]
								auto occ_2_sided = query_2_sided_recursive(o,o+m-1);

								//add all 2 sided occurrences to the result
								for(auto h : occ_2_sided)
									occ.push_back(h);

							}

						}

					}

				}

			}

			//extend search
			fwd_range = fm_index_fwd->LF(fwd_range,pattern[m-i-1]);

		}

		assert(occ.size() == count(pattern));

		return occ;

	}

	vector<ulint> locate_light(string &pattern, bool optimize = true){

		//first, find range of pattern (reversed because we use the reverse index)
		auto range = fm_index_rev->count(pattern,true);

		ulint n_occ = range.second < range.first ? 0 : range.second - range.first + 1;

		//we will extract text in range and put here all primary occurrences found
		set<ulint> primary_occurrences;
		extract_occ_fwd(range,primary_occurrences,pattern.size());

		vector<ulint> occ;

		for(auto o1:primary_occurrences) occ.push_back(o1);

		for(auto o1:primary_occurrences){

			if((not optimize) or occ.size() < n_occ){

				for(auto o2 : query_2_sided_recursive(o1,o1+pattern.size()-1))
					occ.push_back(o2);

			}

		}

		assert(n_occ == count(pattern));

		return occ;

	}

	/*
	 * input: F range containing a unique character, a set where to put found
	 * primary occurrences, number i of FL steps left to do
	 */
	void extract_occ_fwd(range_t rn, set<ulint>& occ, ulint i, bool ignore_interval=true){

		if(i==0) return;

		//search and extract SA samples of marked positions inside interval rn
		//we ignore SA samples at last pattern position (ignore_interval)
		//because we want to treat occurrences that are suffix of a LZ phrase
		//as secondary occurrences.
		if(not ignore_interval) find_occ(rn, occ, i);

		//map range from F to L
		vector<range_t> ranges_on_F = fm_index_rev->FL(rn);

		//recursion
		for(auto r : ranges_on_F){

			extract_occ_fwd(r,occ,i-1,false);

		}

	}

	void find_occ(range_t rn, set<ulint>& occ, ulint i){

		assert(i>0);
		assert(rn.first>0); //# must not be contained in the range

		assert(rn.first<=last->size());
		assert(rn.second+1<=last->size());
		assert(SA_samples.size() == last->rank(last->size()));
		assert(SA_samples.size() == begin_of_phrase.rank(begin_of_phrase.size()));

		for(ulint j = last->rank(rn.first);j<last->rank(rn.second+1);++j){

			assert(j<SA_samples.size());
			ulint text_position = begin_of_phrase.select(SA_samples[j]);

			//decrement: text_position is the first position of a phrase, but we
			//are interested in the last position of previous phrase
			assert(text_position>0);
			text_position--;

			assert(text_position >= (i - 1));
			occ.insert( text_position - (i - 1) );

		}

	}

	void build(string &input, bool verbose){

		if(verbose)
			cout << "Initializing structures for the LZ77 parser ... " << endl;

		set<pair<uchar,ulint> > alphabet_and_freq;

		{
			std::istringstream iss(input);
			alphabet_and_freq = lz77_parser_t::get_alphabet_and_frequencies(iss);
		}

		alphabet = vector<uchar>();

		for(auto c : alphabet_and_freq)
			alphabet.push_back(c.first);

		std::sort(alphabet.begin(),alphabet.end());

		std::istringstream istr(input);

		auto parser = lz77_parser_t(istr,alphabet_and_freq,128,true);

		//size of the text
		ulint text_length = input.size();

		bwt_length=text_length + 1;

		{

			//this vector will contain true in positions that are the first character of a phrase
			vector<bool> phrase_start_vec;

			ulint z=0;	//number of LZ phrases

			//this vector will contain the 2D points for range search, one per phrase
			vector<pair<point_2d_t, ulint> > two_sided_points;

			ulint phrase_start_position = 0; //start position of the current phrase

			//start parsing
			while(not parser.eof()){

				auto token = parser.get_token(); //get a factor

				assert((not token.start_position_is_defined) or token.start_position < phrase_start_position);

				//append a true: we are at the beginning of a phrase
				phrase_start_vec.push_back(true);
				//append false for the remaining phrase characters
				for(ulint i=1;i<token.phrase.size();++i)
					phrase_start_vec.push_back(false);

				//for all factors with a start position, insert the corresponding point in the
				//2-sided range search structure
				if(token.start_position_is_defined){

					//pair: <start_position, end_position> of the copied string
					point_2d_t point = {(ulint)token.start_position,
										(ulint)(token.start_position + token.phrase.size()-1)};

					//insert the point in the vector of points.
					//we associate phrase rank (z) to each point
					two_sided_points.push_back({point,z});

				}

				z++;	//increment number of phrases

				//compute start position of next phrase
				phrase_start_position += token.phrase.size();

			}

			if(verbose)
				cout << "Building two-sided range data structure ... " << flush;

			//build range search structure. Do not re-map y coordinates: saves z words of space
			lz_2_sided_range_search = new range_search_t(two_sided_points,true,false);

			if(verbose)
				cout << "done!" << endl;

			assert(phrase_start_vec.size()==text_length);

			begin_of_phrase = sparse_bitvector_t(phrase_start_vec);

		}

		// Now build FM indexes

		// reverse index:

		if(verbose)
			cout << "Building reverse RLBWT ... " << flush;

		{
			string rev_text;

			for(ulint i=0;i<text_length;++i)
				rev_text.push_back(input[text_length - i - 1]);

			fm_index_rev = new fm_index_t(rev_text);
		}

		if(verbose)
			cout << "done!" << endl;

		// build array 'last'

		{

			if(verbose)
				cout << "Building sparse bitvector marking end of phrases ... " << flush;

			//first, build a standard bitvector where phrase ends (on the F column of the BWT of the
			//reversed text) are marked with a 1. Then, convert it to a sparse bitvector.
			vector<bool> last_bv(bwt_length, false);

			//position on rev BWT: initially, position 0 == position 0 on text
			ulint bwt_pos = 0;

			//navigate the reverse BWT, reading the (forward) text from first to last character
			for(ulint i=0;i<text_length;++i){

				//if i is the first position of a phrase and i>0, then the current BWT position is last
				//position of a phrase on column F
				if(begin_of_phrase[i] and i>0){

					assert(bwt_pos < last_bv.size());

					last_bv[bwt_pos] = true;		//mark the position

				}

				bwt_pos = fm_index_rev->LF(bwt_pos);		//update position in the BWT

			}

			//now bwt_pos is the position on BWT of #. With another step we go back to the
			//position of # on the F column, which must be equal to 0 if everything is correct
			bwt_pos = fm_index_rev->LF(bwt_pos);
			assert(bwt_pos==0);
			last_bv[bwt_pos] = true;

			//now convert the bitvector to a sparse bitvector

			last = new sparse_bitvector_t(last_bv);

			if(verbose)
				cout << "done!" << endl;

		}

		//populate SA_samples
		if(type==light){

			if(verbose)
				cout << "Sampling suffix array on LZ phrase borders ... " << flush;

			assert(text_length>0);
			//in SA_samples we write phrase ranks (on the text)

			ulint z = begin_of_phrase.rank(begin_of_phrase.size());

			assert( z == last->rank(last->size()));

			ulint bitlength =  64 - __builtin_clzll(z);

			SA_samples = packed_view<vector>(bitlength, z);
			assert(SA_samples.width()==bitlength);

			//position on F column of rev BWT
			ulint F_pos = 0;

			//navigate the reverse BWT, reading the (forward) text from first to last character
			//i = text position corresponding to F_pos
			for(ulint i=0;i<text_length;++i){

				//if i is the first position of a phrase and i>0, then the current BWT position
				//is last position of a phrase
				if(begin_of_phrase[i] and i>0){

					assert(F_pos < last->size());
					assert(last->at(F_pos));

					SA_samples[last->rank(F_pos)] = begin_of_phrase.rank(i);

				}

				F_pos = fm_index_rev->LF(F_pos);		//update position in the BWT

			}


			if(verbose)
				cout << "done!" << endl;

		}

		//forward index

		if(verbose)
			cout << "Building forward RLBWT ... " << flush;

		fm_index_fwd = new fm_index_t(input);

		//detect position of terminator character in the BWT
		terminator_pos = fm_index_fwd->get_terminator_position();

		if(verbose)
			cout << "done!" << endl;

		//Now build the interval data structure on the suffix tree nodes corresponding to the LZ factors

		if(type!=light){

			if(verbose)
				cout << "Building interval data structure on the subset of ST nodes ... " << flush;

			vector<pair<ulint, ulint> > bwt_intervals;

			//the bwt intervals corresponding to the factors
			//we use an ordered set to keep only one occurrence of each interval
			//here, ordering is lexicographic (in the cooordinates of the pairs)
			set<pair<ulint, ulint> > bwt_intervals_set;

			//initial range: full
			pair<ulint, ulint> range = {0,bwt_length-1};

			//scan text from end to begin, except first character: in this way we do not include
			//the first phrase. This is correct since before the first phrase there are no characters, so
			//the 2D point is not well-defined.
			for(ulint i=0;i<text_length-1;++i){

				//extend range with current text character
				range = fm_index_fwd->LF(range, input[text_length-i-1]);

				//if this character is the begin of a phrase, then current range is a the range of a factor
				if(begin_of_phrase[text_length-i-1]){

					//insert range in the set
					bwt_intervals_set.insert(range);
					//reset interval to full
					range = {0,bwt_length-1};

				}

			}

			//build structure
			st_sub = new st_subset(bwt_length,bwt_intervals_set);

			//convert intervals set to a vector
			bwt_intervals = vector<pair<ulint, ulint> >(bwt_intervals_set.size());
			std::copy(bwt_intervals_set.begin(), bwt_intervals_set.end(), bwt_intervals.begin());

			auto comp = [](const pair<ulint,ulint> &A, const pair<ulint,ulint> &B){

				//A precedes B
				bool preceded = A.second < B.first;

				//B is contained in A and the two intervals are distinct
				bool contained = 	(B.first >= A.first and B.second < A.second) or
									(B.second <= A.second and B.first > A.first);

				return preceded or contained;

			};

			//sort the vector according to the lexicographic order of the factors
			std::sort(bwt_intervals.begin(), bwt_intervals.end(), comp);

			if(verbose)
				cout << "done!" << endl;

			if(verbose)
				cout << "Building four-sided range data structure ... " << flush;

			//scan forward parsed text from last to first character and store in a vector the rankings
			//of begin of phrases

			//one ranking for each phrase, except the first phrase
			vector<ulint> fwd_ranks;

			//initial range: full
			range = {0,bwt_length-1};

			//scan the parsed text from last to second character and search factors
			//we skip the first character: in this way we do not create a 2D point for the 1st factor
			//(which does not have characters on the left)
			for(ulint i=0;i<text_length-1;++i){

				//search character: extend current search
				range = fm_index_fwd->LF(range, input[text_length-i-1]);

				//if this character is the begin of a phrase, then current range is a the range of a factor
				if(begin_of_phrase[text_length-i-1]){

					auto comp = [](const pair<ulint,ulint> &A, const pair<ulint,ulint> &B){

						//A precedes B
						bool preceded = A.second < B.first;

						//B is contained in A and the two intervals are distinct
						bool contained = 	(B.first >= A.first and B.second < A.second) or
											(B.second <= A.second and B.first > A.first);

						return preceded or contained;

					};

					auto it = std::lower_bound(bwt_intervals.begin(), bwt_intervals.end(),range,comp);

					//the interval must be inside the set
					assert(it!=bwt_intervals.end());
					assert(*it==range);

					ulint rank_in_st_nodes =  it - bwt_intervals.begin();

					//insert lexicographic rank in the set
					fwd_ranks.push_back(rank_in_st_nodes);

					//reset interval to full
					range = {0,bwt_length-1};

				}

			}

			//now scan reverse text and compute 2D points for four-sided range search structure

			//the set of pairs <point, text position>. Test position is the first position of the beginning
			//of the right factor
			vector<pair<point_2d_t, ulint> > points;

			//current range in the BWT of the reversed text
			range = {0,0};

			//current position in fwd_ranks
			ulint fwd_ranks_ptr = 0;

			//phrase number 0 is not used
			ulint phrase_nr = 1;

			//scan forward text from first to last character
			for(ulint i=0;i<text_length;++i){

				//if current position is the begin of a phrase AND current position is not the first
				//then range is the position in rev. BWT of the end of the left factor
				if(begin_of_phrase[i] and i>0){

					//make sure that there is a 1 in this BWT position of last
					assert(range.first<=last->size());
					assert(last->at(range.first));
					assert(fwd_ranks_ptr<fwd_ranks.size());

					//the 2D point: rank on reverse, rank on forward
					assert(last->rank(range.first)>0);

					//since rev ranks start from 1 (we never use last phrase as left phrase
					//in a 2D point AND rank of last phrase is 0 in bitvector last), we
					//subtract 1 so that x coordinates start from 0
					point_2d_t p = {(ulint)last->rank(range.first)-1,
									(ulint)fwd_ranks[fwd_ranks.size()-fwd_ranks_ptr-1]};

					fwd_ranks_ptr++;

					//insert new point
					points.push_back( {p,phrase_nr++} );

				}

				//extend search with new character
				range = fm_index_rev->LF(range, input[i]);

			}

			//build 4 sided range structure. Do not re-map coordinates since they are already in rank-space.
			lz_4_sided_range_search = new range_search_t(points, false, false);

		}

		//if type == bidirectional, we don't need the reverse index anymore.
		//if type == light, we don't need the forward index anymore.
		//if type == full, we keep both indexes
		if(type == bidirectional){

			delete fm_index_rev;
			fm_index_rev = NULL;

		}else if(type == light){

			delete fm_index_fwd;
			fm_index_fwd = NULL;

		}

		if(verbose)
			cout << "done!" << endl;

	}

	/*
	 * rangeCount: count number of characters equal to c in the bwt range
	 */
	ulint rangeCount(pair<ulint,ulint> bwt_range, uchar c){

		auto new_range = fm_index_fwd->LF(bwt_range,c);

		if(new_range.second<new_range.first)
			return 0;

		return (new_range.second-new_range.first)+1;

	}

	/*
	 * Synchronize range on forward and backward BWTs using only the forward BWT
	 *
	 * \param fwd_bwt_range  range in the forward BWT of word W^rev
	 * \param rev_bwt_range  range in the reverse BWT of word W
	 * \param c          character
	 *
	 * \return  range in the reverse BWT of Wc
	 *
	 */
	pair<ulint,ulint> synchronize(	pair<ulint,ulint> fwd_bwt_range,
									pair<ulint,ulint> rev_bwt_range,
									uchar c){

		//we do not accept reserved 0x0 character
		assert(c>0);

		//number of c in the fwd BWT interval
		ulint occ_c = rangeCount(fwd_bwt_range, c);

		//no extension is possible: return empty interval
		if(occ_c==0)
			return {1,0};

		//left bound of result
		ulint l = rev_bwt_range.first;

		//check if terminator is inside this BWT interval
		if(terminator_pos >= fwd_bwt_range.first and terminator_pos <= fwd_bwt_range.second)
			l++;

		//compute left interval
		for(uint i=0;i<alphabet.size() and alphabet[i]<c;++i)
			l += rangeCount(fwd_bwt_range,alphabet[i]);

		ulint r = (l + occ_c) - 1;

		return {l,r};

	}

	/*
	 * simulate reverse BWT using forward BWT and return range of s^rev
	 *
	 * \param s input string
	 * \return range of s^rev in the reverse BWT
	 *
	 */
	pair<ulint,ulint> bidirectional_LF(string s){

		ulint m = s.size();
		auto fwd_range = fm_index_fwd->get_char_range((uchar)s[m-1]);
		pair<ulint,ulint> rev_range(fwd_range);

		for(ulint i=1;i<m;++i){

			auto c = s[m-i-1];

			rev_range = synchronize(fwd_range,rev_range,c);
			fwd_range = fm_index_fwd->LF(fwd_range,c);

			if(rev_range.second < rev_range.first)
				return rev_range;

		}

		return rev_range;

	}

	/*
	 * find secondary occurrences copied from T[begin,...,end].
	 * Calls recursively itself until no more occurrences can be found
	 */
	vector<ulint> query_2_sided_recursive(ulint begin, ulint end){

		vector<ulint> result;

		assert(end>=begin);

		//pattern length
		ulint m = end-begin+1;

		//launch a 2 sided range search
		//result: a vector of pairs < <b,e>, i >, meaning that T[b,e] is copied in T[i,...,i+(e-b)]

		assert(end<bwt_length-1);

		auto points = lz_2_sided_range_search->two_sided_range_search_ul({begin,end});

		for(auto p : points){

			ulint b = p.first.x;

			//p.second is a phrase rank. We need to re-map this with a select query
			assert(p.second < begin_of_phrase.rank(begin_of_phrase.size()));
			ulint i = begin_of_phrase.select(p.second);

			ulint shift = begin - b;
			ulint j = i+shift; //new occurrence of the pattern in position j=i+shift

			//call recursively on the occurrence found
			assert(j+m-1<bwt_length-1);
			assert(j>begin);
			auto recursive_occ = query_2_sided_recursive(j,j+m-1);

			//append all found occurrences to the final result
			result.push_back(j);

			for(auto x : recursive_occ)
				result.push_back(x);

		}

		return result;

	}

	ulint bwt_length = 0;		//length of the BWT (i.e. 1+text length)
	ulint terminator_pos = 0;	//position of terminator character in the forward BWT
	vector<uchar> alphabet;		//alphabet characters sorted in lexicographic order

	range_search_t * lz_4_sided_range_search = NULL; 	// points representing (lexicographic ordering of)
														//adjacent lz factors
	range_search_t * lz_2_sided_range_search = NULL; 	// points with occurring positions of factors

	fm_index_t * fm_index_fwd = NULL;	// O(R) words. FM index for the forward text, where factors
										//are separated with a SEPARATOR character
	fm_index_t * fm_index_rev = NULL;	// O(R) words. FM index for the reversed text

	// z words. marks with a 1 positions in the F column of rev BWT
	//that correspond to end of phrases. Equivalently, marks with a 1 positions
	//on the L column of rev BWT that correspond to begin of phrases (except text position
	//0 == BWT position 0 which is not marked
	sparse_bitvector_t * last  = NULL;

	//this sparse bitvector marks with a 1 the first position of each phrase
	sparse_bitvector_t begin_of_phrase;

	st_subset * st_sub = NULL;	//interval data structure on the subset of nodes of the suffix tree
								//corresponding to the LZ phrases

	//used only if index is of light type (otherwise, samples
	//are stored inside lz_4_sided_range_search). Synchronized with bitvector last.
	//start of phrases on L column of rev BWT (except text position 0 == rev BWT position 0)
	packed_view<vector> SA_samples;

	lz_rlbwt_type type = full;

};//class lz_rlbwt

}//namespace lzrlbwt

#endif /* LZrlbwt_H_ */
