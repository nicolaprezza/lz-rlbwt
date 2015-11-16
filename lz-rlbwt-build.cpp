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

#include <iostream>
#include <lz_rlbwt.h>
#include <utils.h>

using namespace lzrlbwt;
using namespace std;

string out_basename=string();
string input_file=string();

//default: light index
bool bidirectional_opt=false;
bool light_opt=true;

void help(){
	cout << "lz-rlbwt-build" << endl << endl;
	cout << "Usage: lz-rlbwt-build [options] <input_file_name>" << endl;
	cout << "   -o <basename>         use 'basename' as prefix for all index files. Default: basename is the specified input_file_name"<<endl<<endl;
	cout << "   -f, --full            build 2 RLBWT and 2 geometric range structures." << endl;
	cout << "   -b, --bidirectional   build 1 RLBWT and 2 geometric range structures." << endl;
	cout << "   -l, --light           build 1 RLBWT and 1 geometric range structure. DEFAULT." << endl <<endl;
	cout << "   <input_file_name>     input text file." << endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-o")==0){

		if(ptr>=argc-1){
			cout << "Error: missing parameter after -o option." << endl;
			help();
		}

		out_basename = string(argv[ptr]);
		ptr++;

	}else if(s.compare("--bidirectional")==0 or s.compare("-b")==0){

		bidirectional_opt = true;
		light_opt = false;

	}else if(s.compare("--light")==0 or s.compare("-l")==0){

		light_opt = true;
		bidirectional_opt = false;

	}else if(s.compare("--full")==0 or s.compare("-f")==0){

		bidirectional_opt = false;
		light_opt = false;

	}else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

void build_lz_rlbwt(string &input, vector<lzrlbwt::lzrlbwt_options> opt){

	auto idx = lz_rlbwt<>(input,opt);
	idx.save_to_file(out_basename, true);

}

int main(int argc, char** argv){

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    auto t1 = high_resolution_clock::now();

	//parse options

    bidirectional_opt=false;
    out_basename=string();
    input_file=string();
	int ptr = 1;

	if(argc<2) help();

	while(ptr<argc-1)
		parse_args(argv, argc, ptr);

	input_file = string(argv[ptr]);

	if(out_basename.compare("")==0)
		out_basename = string(input_file);

	cout << "Building lz-rlbwt of input file '" << input_file << "'" << endl;
	cout << "Prefix '" << out_basename << "' will be used for the index files." << endl;

	string input;

	{

		std::ifstream fs(input_file);
		std::stringstream buffer;
		buffer << fs.rdbuf();

		input = buffer.str();

	}

	vector<lzrlbwt::lzrlbwt_options> opt;

	if(bidirectional_opt)
		opt = {verbose_out, bidirectional_bwt};
	else if(light_opt)
		opt = {verbose_out, light_index};
	else
		opt = {verbose_out};

	auto idx = lz_rlbwt<>(input,opt);
	idx.save_to_file(out_basename, true);

	printRSSstat();

	auto t2 = high_resolution_clock::now();
	ulint total = duration_cast<duration<double, std::ratio<1>>>(t2 - t1).count();
	cout << "Build time : " << get_time(total) << endl;

}
