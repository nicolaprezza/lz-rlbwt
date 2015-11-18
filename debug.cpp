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

using namespace std;
using namespace lzrlbwt;

int main(int argc, char** argv){

	cout << "1" << endl;
	gap_bv g;
	cout << "2" << endl;

	g = gap_bv();
	cout << "3" << endl;

	g.insert(0,true);
	cout << "4" << endl;

}

