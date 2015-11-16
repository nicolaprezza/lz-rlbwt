/*
 * rlbwt.h
 *
 *  Created on: May 20, 2015
 *      Author: nicola
 *
 *  A run-length encoded BWT. This class implements LF function and count
 *  TODO: locate, display
 *
 *  Note: all ranges are inclusive
 *
 *  Note: caracters 0x0 and 0x1 are reserved. 0x0 is used as string terminator in the internal
 *  data structures (required by sdsl::wt_huff), while 0x1 is the BWT text terminator.
 *
 */

#ifndef RLBWT_H_
#define RLBWT_H_

#include <definitions.h>
#include <sparse_sd_vector.h>
#include <rle_string.h>

namespace lzrlbwt {

template<
	class rle_string_t	= rle_string<sparse_sd_vector<> >	//run-length encoded string
>
class rlbwt{

public:

	rlbwt(){}

	/*
	 * Build a run-length BWT from the input string
	 * input must not contain 0x0 and 0x1 characters
	 */
	rlbwt(string &input){

		assert(not contains_reserved_chars(input));
		string bwt_s = build_bwt(input);

		bwt = rle_string_t(bwt_s);

		//build F column
		F = vector<ulint>(256,0);
		for(uchar c : bwt_s)
			F[c]++;

		for(ulint i=255;i>0;--i)
			F[i] = F[i-1];

		F[0] = 0;

		for(ulint i=1;i<256;++i)
			F[i] += F[i-1];

		for(ulint i=0;i<bwt_s.size();++i)
			if(bwt_s[i]==TERMINATOR)
				terminator_position = i;

	}

	/*
	 * get full BWT range
	 */
	range_t full_range(){

		//inclusive range
		return {0,bwt_size()-1};

	}

	uchar operator[](ulint i ){
		return bwt[i];
	}

	/*
	 * \param r inclusive range of a string w
	 * \param c character
	 * \return inclusive range of cw
	 */
	range_t LF(range_t rn, uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		//number of c before the interval
		ulint c_before = bwt.rank(rn.first,c);

		//number of c inside the interval rn
		ulint c_inside = bwt.rank(rn.second+1,c) - c_before;

		//if there are no c in the interval, return empty range
		if(c_inside==0) return {1,0};

		ulint l = F[c] + c_before;

		return {l,l+c_inside-1};

	}

	/*
	 * \param rn inclusive F-range containing a single character (say, a)
	 * \return 	disjoint inclusive ranges on column L whose union are all a's contained
	 * 			in the input range
	 *
	 * 	complexity: |result| * FL_cost
	 * 				where result is the output vector and FL_cost
	 * 				is the cost of FL mapping (log n/R with elias-fano RLE string)
	 *
	 */
	vector<range_t> FL(range_t rn){

		//F-range rn must contain only one character
		assert(uniq_char(rn));

		//unique character in F-range rn
		uchar c = F_at(rn.first);
		auto l = rn.first;
		auto r = rn.second;

		//break range: given a range <l',r'> on BWT and a character c, this function
		//breaks <l',r'> in maximal sub-ranges containing character c.
		return bwt.break_range({FL(l,c),FL(r,c)},c);

	}

	//backward navigation of the BWT
	ulint LF(ulint  i){

		auto c = bwt[i];
		return F[c] + bwt.rank(i,c);

	}

	//forward navigation of the BWT
	ulint FL(ulint  i){

		//i-th character in first BWT column
		auto c = F_at(i);

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	//forward navigation of the BWT, where for efficiency we give c=F[i] as input
	ulint FL(ulint  i, uchar c){

		//i-th character in first BWT column
		assert(c == F_at(i));

		//this c is the j-th (counting from 0)
		ulint j = i - F[c];

		return bwt.select(j,uchar(c));

	}

	/*
	 * access column F at position i
	 */
	uchar F_at(ulint i){

		ulint c = (upper_bound(F.begin(),F.end(),i) - F.begin()) - 1;
		assert(c<256);
		assert(i>=F[c]);

		return uchar(c);

	}

	/*
	 * Return BWT range of character c
	 */
	range_t get_char_range(uchar c){

		//if character does not appear in the text, return empty pair
		if((c==255 and F[c]==bwt_size()) || F[c]>=F[c+1])
			return {1,0};

		ulint l = F[c];
		ulint r = bwt_size()-1;

		if(c<255)
			r = F[c+1]-1;

		return {l,r};

	}

	/*
	 * Return BWT range of pattern P. If reverse = true, search P reversed
	 */
	range_t count(string &P, bool reverse = false){

		auto range = full_range();
		ulint m = P.size();

		for(ulint i=0;i<m and range.second>=range.first;++i)
			range = LF(range,(reverse?P[i]:P[m-i-1]) );

		return range;

	}

	/*
	 * get number of runs in the BWT (terminator character included)
	 */
	ulint number_of_runs(){
		return bwt.number_of_runs();
	}

	/*
	 * get terminator (0x1) position in the BWT
	 */
	ulint get_terminator_position(){
		return terminator_position;
	}

	/*
	 * get BWT in string format
	 */
	string get_bwt(){
		return bwt.toString();
	}

	/* serialize the structure to the ostream
	 * \param out	 the ostream
	 */
	ulint serialize(std::ostream& out){

		ulint w_bytes = 0;

		assert(F.size()>0);
		assert(bwt.size()>0);

		out.write((char*)&terminator_position,sizeof(terminator_position));
		out.write((char*)F.data(),256*sizeof(ulint));

		w_bytes += sizeof(terminator_position) + 256*sizeof(ulint);

		w_bytes += bwt.serialize(out);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&terminator_position,sizeof(terminator_position));

		F = vector<ulint>(256);
		in.read((char*)F.data(),256*sizeof(ulint));

		bwt.load(in);

	}

	/*
	 * save the structure to the path specified.
	 * \param path_prefix prefix of the index files. suffix ".rlbwt" will be automatically added
	 */
	void save_to_disk(string path_prefix){

		string path = string(path_prefix).append(".rlbwt");

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	/*
	 * load the structure from the path specified.
	 * \param path_prefix prefix of the index files. suffix ".rlbwt" will be automatically added
	 */
	void load_from_disk(string path_prefix){

		string path = string(path_prefix).append(".rlbwt");
		std::ifstream in(path);
		load(in);
		in.close();

	}

	ulint text_size(){
		return bwt.size()-1;
	}

	ulint bwt_size(){
		return bwt.size();
	}



	uchar get_terminator(){
		return TERMINATOR;
	}

	string get_extension(){
		return string(".rlbwt");
	}

private:

	/*
	 * check if range rn on column F contains a
	 * single character
	 */
	bool uniq_char(range_t rn){

		for(ulint i=0;i<256;++i){

			ulint l = F[i];
			ulint r = (i==255?bwt_size():F[i+1]);

			if(rn.first >= l and rn.second < r ) return true;

		}

		return false;

	}

	/*
	 * builds BWT of input string using SE_SAIS algorithm
	 * uses 0x1 character as terminator
	 *
	 */
	static string build_bwt(string &s){

		string bwt_s;

	    cache_config cc;

	    int_vector<8> text(s.size());
	    assert(text.size()==s.size());

	    for(ulint i=0;i<s.size();++i)
	    	text[i] = (uchar)s[i];

	    assert(text.size()==s.size());

	    append_zero_symbol(text);

	    store_to_cache(text, conf::KEY_TEXT, cc);

	    construct_config::byte_algo_sa = SE_SAIS;
	    construct_sa<8>(cc);

	    //now build BWT from SA
	    int_vector_buffer<> sa(cache_file_name(conf::KEY_SA, cc));

	    {

	        for (ulint i=0; i<sa.size(); i++){
	            auto x = sa[i];

	            assert(x<=text.size());

	            if ( x > 0 )
	            	bwt_s.push_back((uchar)text[x-1]);
	            else
	            	bwt_s.push_back(TERMINATOR);

	        }

	    }

	    sdsl::remove(cache_file_name(conf::KEY_TEXT, cc));
	    sdsl::remove(cache_file_name(conf::KEY_SA, cc));

	    return bwt_s;

	}

	static bool contains_reserved_chars(string &s){

		for(auto c : s)
			if(c == 0 or c == 1)
				return true;

		return false;

	}

	static const uchar TERMINATOR = 1;

	rle_string_t bwt;

	ulint terminator_position = 0;

	//F column of the BWT
	vector<ulint> F;

};

typedef rlbwt<rle_string_sd> rlbwt_sd;		//Elias-Fano encoding of gap lengths (fast)

}

#endif /* RLBWT_H_ */
