/*
 *  This file is part of lz-rlcsa.
 *  Copyright (c) by
 *  Nicola Prezza <nicolapr@gmail.com>,
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, and Mathieu Raffinot
 *
 *   lz-rlcsa is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   lz-rlcsa is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * utils_lzrlcsa.h
 *
 *  Created on: Mar 31, 2015
 *      Author: nicola
 */

#include <iostream>

using namespace std;

#ifndef UTILS_LZRLCSA_H_
#define UTILS_LZRLCSA_H_

using ulint = uint64_t;

string get_time(uint64_t time){

	stringstream ss;

	if(time>=3600){

		uint64_t h = time/3600;
		uint64_t m = (time%3600)/60;
		uint64_t s = (time%3600)%60;

		ss  << time << " seconds. ("<< h << "h " << m << "m " << s << "s" << ")";

	}else if (time>=60){

		uint64_t m = time/60;
		uint64_t s = time%60;

		ss << time << " seconds. ("<< m << "m " << s << "s" << ")";

	}else{

		ss << time << " seconds.";

	}

	return ss.str();

}

//parse pizza&chilli patterns header:
void header_error(){
	cout << "Error: malformed header in patterns file" << endl;
	cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << endl;
	exit(0);
}

ulint get_number_of_patterns(string header){

	ulint start_pos = header.find("number=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

ulint get_patterns_length(string header){

	ulint start_pos = header.find("length=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}



#endif /* UTILS_LZRLCSA_H_ */
