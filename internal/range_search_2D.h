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
 * range_search_2D.h
 *
 *  Created on: Mar 24, 2015
 *      Author: Nicola Prezza
 *
 *  Description: a general 2D range search data structure based on wt_int.
 *  Differently than wt_int, this class permits the definition of points with same x coordinate,
 *  and even points having the same x and y coordinates.
 *  Moreover, the space of coordinates is bounded only by the memory word size: coordinates are re-mapped
 *  to rank space so that space of the structure is always O(z) words and operations cost O(log z),
 *  where z is the number of points.
 *
 *	The class permits to associate an object (of numeric type) with each 2D point. Range search functions
 *	return the objects lying in the range.
 *  The class is a template on the word type: to save space, use the smallest integer type.
 *
 */

#ifndef RANGE_SEARCH_2D_H_
#define RANGE_SEARCH_2D_H_

#include "sdsl/wavelet_trees.hpp"
#include <vector>
#include <set>
#include <packed_view.h>
#include "../internal/sparse_sd_vector.h"
#include "../internal/succinct_bit_vector.h"

using namespace std;
using namespace sdsl;
using namespace bv;

namespace lzrlbwt{

struct point_2d_t{

	ulint x;
	ulint y;

};

template<
			class sparse_bitvector_t = sparse_sd_vector<>,					// elias-fano sparse bitvector
			class succinct_bitvector_t = succinct_bit_vector				//succinct bitvector with constant-time rank and access
>
class range_search_2D{

public:

	/* Empty constructor
	 */
	range_search_2D(){};

	/* Constructor with 2D points associated to objects.
	 * The user can specify if the x and y coordinates should be remapped to
	 * rank space: this guarantees O(log z) operations for range search (z=number of points),
	 * but requires additional z words for each remapped coordinate (at most 2z words in total)
	 *
	 * Warning: if remap_x = false, then the x coordinates of the input points must
	 * be exactly 0, ..., |points|-1 (not necessarily in order). Otherwise, this constructor will produce a
	 * failed assertion.
	 *
	 * \param points	a vector containing pairs of point, object
	 * \param remap_x	remap x coordinate to rank space?
	 * \param remap_y	remap y coordinate to rank space?
	 */
	range_search_2D(vector<pair<point_2d_t, ulint> > points, bool remap_x=true, bool remap_y=true){

		using pair_t = pair<point_2d_t, ulint>;

		//sort the points with respect to their x coordinate
		std::sort(	points.begin(),
					points.end(),
					[](const pair_t &left, const pair_t &right) {
		    			return left.first.x < right.first.x;
					}
		);

		ulint max_val = 0;
		for(ulint i=0;i<points.size();++i)
			if(points[i].second>max_val)
				max_val = points[i].second;

		obj_bitlength = 64 - __builtin_clzll(max_val);
		assert(obj_bitlength>0);

		//compute max x and y values
		min_x = points[0].first.x;
		max_x = points[points.size()-1].first.x;
		max_y = 0;
		for(ulint i=0;i<points.size();++i)		//find max_y
			if(points[i].first.y>max_y)
				max_y = points[i].first.y;
		min_y = max_y;
		for(ulint i=0;i<points.size();++i)		//find max_y
			if(points[i].first.y<min_y)
				min_y = points[i].first.y;


		if(remap_x){

			{
				vector<bool> x_bools(max_x+1,false);

				for(auto p : points)
					x_bools[p.first.x] = true;

				remap_x_bv = sparse_bitvector_t(x_bools);
			}

			//there could be repeated x coordinates, and we must keep track of this.

			{
				vector<bool> dup_x;

				ulint prev_x=max_x+1;
				for(auto p : points){

					auto x = p.first.x;

					if(x==prev_x){

						dup_x.push_back(false);

					}else{

						dup_x.push_back(true);
						prev_x=x;

					}

				}

				duplicates_x = succinct_bitvector_t(dup_x);
				assert(duplicates_x.size()>0);

			}


		}else{

			//check if the x coordinates are exactly 0, ..., |points|-1
			bool correct_x_coordinates = points[0].first.x == 0;

			for(ulint i=1;i<points.size() and correct_x_coordinates;++i)
				correct_x_coordinates = correct_x_coordinates and
										points[i].first.x == (points[i-1].first.x + 1);

			assert(("Error: creating a 2D range structure with non-consecutive x coordinates. If in doubt, call the constructor with remap_x=true",correct_x_coordinates));

		}

		if(remap_y){

			//extract y coordinates, sort them and convert to sparse vector
			set<ulint> sorted_y_coord;
			for(ulint i=0;i<points.size();++i)
				sorted_y_coord.insert(points[i].first.y);

			vector<bool> y_bools(max_y+1,false);

			for(auto y : sorted_y_coord)
				y_bools[y] = true;

			remap_y_bv = sparse_bitvector_t(y_bools);

		}

		//objects are stored in sorted order with respect to x coordinate
		obj = packed_view<vector>(obj_bitlength,points.size());
		for(ulint i=0;i<points.size();++i)
			obj[i] = points[i].second;

		//compute maximum value that will be stored in the wavelet tree
		ulint max_y_in_wt = 0;

		if(remap_y)
			max_y_in_wt = remap_y_bv.number_of_1() - 1;		//rank space
		else
			max_y_in_wt = max_y;	//corresponds to max_y

		//compute depth of the WT (max number of bits of an integer)
		uint8_t wt_depth = 0;
		while(max_y_in_wt>0){
			max_y_in_wt = max_y_in_wt >> 1;
			wt_depth++;
		}

		if(wt_depth==0) wt_depth=1;	//degenerate case: all y coordinates are the same

		//create the integer vector that will contain point coordinates
		int_vector<> iv(points.size(),0,wt_depth);

		//insert (eventually remapped) y values
		for(ulint i=0;i<points.size();++i){

			ulint y; //the y coordinate

			if(remap_y){ //map to rank space

				y = remap_y_bv.rank(points[i].first.y);

				//y = std::lower_bound( y_coords.begin(),y_coords.end(),points[i].first.y) - y_coords.begin();

			}else{ //keep the coordinate as is

				y = points[i].first.y;

			}

			iv[i] = y;

		}

		//build wavelet tree
		construct_im(wt, iv);

		assert(remap_x_bv.size()==0 or duplicates_x.size()>0);

	}

	/* 4-sided range search
	 * \param x_bounds	x left and right boundaries (inclusive) of the range
	 * \param y_bounds	y left and right boundaries (inclusive) of the range
	 * \return	vector of objects (associated to points) lying in the range
	 */
	vector<pair<point_2d_t,ulint> > four_sided_range_search(pair<ulint,ulint> x_bounds, pair<ulint,ulint> y_bounds){

		assert(remap_x_bv.size()==0 or duplicates_x.size()>0);

		x_bounds.second = (x_bounds.second > max_x?max_x:x_bounds.second);
		y_bounds.second = (y_bounds.second > max_y?max_y:y_bounds.second);

		vector<pair<point_2d_t,ulint> > result;

		//the (eventually remapped) range
		pair<ulint,ulint> real_x_bounds = x_bounds;
		pair<ulint,ulint> real_y_bounds = y_bounds;

		//check if the range is empty
		if( x_bounds.first > max_x or
			x_bounds.second < min_x or
			y_bounds.first > max_y or
			y_bounds.second < min_y or
			x_bounds.first > x_bounds.second or
			y_bounds.first > y_bounds.second)
			return result;

		if(remap_x_bv.size()>0){	//if x coordinate has to be remapped to rank space

			assert(x_bounds.first<=remap_x_bv.size());
			assert(x_bounds.second+1<=remap_x_bv.size());

			ulint first_1_in_range = remap_x_bv.rank(x_bounds.first);
			ulint number_of_1_in_range = remap_x_bv.rank(x_bounds.second+1) - first_1_in_range;

			if(number_of_1_in_range==0){

				real_x_bounds = {1,0};

			}else{

				ulint last_1_in_range = (first_1_in_range + number_of_1_in_range) - 1;

				ulint low = duplicates_x.select(first_1_in_range);
				ulint up=0;

				if(last_1_in_range==duplicates_x.number_of_1()-1)
					up = duplicates_x.size()-1;
				else
					up = duplicates_x.select(last_1_in_range+1)-1;

				real_x_bounds = {low,up};

			}

		}

		if(remap_y_bv.size()>0){ //if y coordinate has to be remapped to rank space

			assert(y_bounds.first<=remap_y_bv.size());
			assert(y_bounds.second+1<=remap_y_bv.size());

			ulint first_1_in_range = remap_y_bv.rank(y_bounds.first);
			ulint number_of_1_in_range = remap_y_bv.rank(y_bounds.second+1) - first_1_in_range;

			if(number_of_1_in_range==0)
				real_y_bounds = {1,0};
			else
				real_y_bounds = {
					first_1_in_range,
					first_1_in_range + number_of_1_in_range -1
				};

		}

		//perform 2D range search on the wavelet tree
		auto wt_points = wt.range_search_2d(real_x_bounds.first,
											real_x_bounds.second,
											real_y_bounds.first,
											real_y_bounds.second).second;

		//the X coordinates of the wt_points are the indexes of the result objects

		for(ulint i=0;i<wt_points.size();++i){

			ulint x = wt_points[i].first;	//coordinates of the point as returned by the wt
			ulint y = wt_points[i].second;

			//if points in wt are in rank-space, then map back to coordinate space
			if(remap_x_bv.size()>0){
				if(duplicates_x[x]==0){
					x = remap_x_bv.select( duplicates_x.rank(x)-1 );
				}else{
					x = remap_x_bv.select( duplicates_x.rank(x) );
				}
			}

			if(remap_y_bv.size()>0)
				y = remap_y_bv.select(y);

			//push back a pair <point, object>
			result.push_back({
				{x,y},						//the point
				obj[ wt_points[i].first ]	//object
			});
		}

		return result;

	}

	/* 2-sided range search on the upper-left corner defined by the point
	 * \param point		the point that defines the 2 sided range (coordinates of the point included)
	 * \return	vector of objects (associated to points) lying in the range [0,point.x] x [point.y, max_y]
	 */
	vector<pair<point_2d_t,ulint> > two_sided_range_search_ul(point_2d_t point){
		return four_sided_range_search({0,point.x},{point.y,max_y});
	}

	/* 2-sided range search on the upper-right corner defined by the point
	 * \param point		the point that defines the 2 sided range (coordinates of the point included)
	 * \return	vector of objects (associated to points) lying in the range [point.x,max_x] x [point.y, max_y]
	 */
	vector<pair<point_2d_t,ulint> > two_sided_range_search_ur(point_2d_t point){
		return four_sided_range_search({point.x, max_x},{point.y,max_y});
	}

	/* 2-sided range search on the bottom-left corner defined by the point
	 * \param point		the point that defines the 2 sided range (coordinates of the point included)
	 * \return	vector of objects (associated to points) lying in the range [0, point.x] x [0,point.y]
	 */
	vector<pair<point_2d_t,ulint> > two_sided_range_search_bl(point_2d_t point){
		return four_sided_range_search({0, point.x},{0,point.y});
	}

	/* 2-sided range search on the bottom-right corner defined by the point
	 * \param point		the point that defines the 2 sided range (coordinates of the point included)
	 * \return	vector of objects (associated to points) lying in the range [point.x,max_x] x [0,point.y]
	 */
	vector<pair<point_2d_t,ulint> > two_sided_range_search_br(point_2d_t point){
		return four_sided_range_search({point.x,max_x},{0,point.y});
	}

	// save the structure to file
	/* \param filename	name of the output file.
	 */
	void save_to_file(string filename){

		std::ofstream out (filename,std::ofstream::binary);
		serialize(out);
		out.close();

	}

	// load index from file
	/* \param filename		filename of the index files
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

		ulint obj_length = obj.size();
		ulint obj_container_size = obj.container().size();

		assert(obj_bitlength>0);

		out.write((char *)&obj_length, sizeof(ulint));
		out.write((char *)&obj_container_size, sizeof(ulint));
		out.write((char *)&obj_bitlength, sizeof(obj_bitlength));

		w_bytes += sizeof(ulint);

		w_bytes += wt.serialize(out);

		out.write((char *)obj.container().data(), obj_container_size*sizeof(ulint));
		w_bytes += obj_container_size*sizeof(ulint);

		assert(remap_x_bv.size()==0 or duplicates_x.size()>0);

		w_bytes += remap_x_bv.serialize(out);
		w_bytes += duplicates_x.serialize(out);
		w_bytes += remap_y_bv.serialize(out);

		out.write((char *)&max_x, sizeof(ulint));
		out.write((char *)&max_y, sizeof(ulint));
		out.write((char *)&min_x, sizeof(ulint));
		out.write((char *)&min_y, sizeof(ulint));

		w_bytes += 4*sizeof(ulint);

		return w_bytes;

	}

	/* load the structure from the istream
	 * \param in the istream
	 */
	void load(std::istream& in) {

		ulint obj_length;
		ulint obj_container_size;

		in.read((char *)&obj_length, sizeof(ulint));
		in.read((char *)&obj_container_size, sizeof(ulint));
		in.read((char *)&obj_bitlength, sizeof(obj_bitlength));

		wt.load(in);

		obj = packed_view<vector>(obj_bitlength,obj_length);
		in.read((char *)obj.container().data(), obj_container_size*sizeof(ulint));

		remap_x_bv.load(in);
		duplicates_x.load(in);
		remap_y_bv.load(in);

		in.read((char *)&max_x, sizeof(ulint));
		in.read((char *)&max_y, sizeof(ulint));
		in.read((char *)&min_x, sizeof(ulint));
		in.read((char *)&min_y, sizeof(ulint));

		assert(remap_x_bv.size()==0 or duplicates_x.size()>0);

	}

	range_search_2D & operator= (const range_search_2D & other) {

		wt = wt_int<>(other.wt);
		obj = packed_view<vector>(other.obj);
		remap_x_bv = sparse_bitvector_t(other.remap_x_bv);
		duplicates_x = succinct_bitvector_t(other.duplicates_x);
		remap_y_bv = sparse_bitvector_t(other.remap_y_bv);

		max_x = other.max_x;
		max_y = other.max_y;
		min_x = other.min_x;
		min_y = other.min_y;
		obj_bitlength = other.obj_bitlength;

	    return *this;
	}

	uint64_t number_of_phrases(){
		return obj.size();
	}

private:

	wt_int<> wt;			//the underlying wavelet tree
	//vector<ulint> obj;

	uint obj_bitlength = 0;  //bitlength of objects associated with the points
	packed_view<vector> obj; //the vector of objects associated with the points

	//vector<ulint> x_coords;	//used to re-map x coordinates. Empty if remap_x=false
	//vector<ulint> y_coords; //used to re-map y coordinates. Empty if remap_y=false

	sparse_bitvector_t remap_x_bv; //used to re-map x coordinates. Empty if remap_x=false
	succinct_bitvector_t duplicates_x; //keep track of duplicate x coordinates

	sparse_bitvector_t remap_y_bv; //used to re-map y coordinates. Empty if remap_y=false

	ulint max_x=0;	//maximum x coordinate
	ulint max_y=0;	//maximum y coordinate
	ulint min_x=0;	//maximum x coordinate
	ulint min_y=0;	//maximum y coordinate

};

}

#endif /* RANGE_SEARCH_2D_H_ */
