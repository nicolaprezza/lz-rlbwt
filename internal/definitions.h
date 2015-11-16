/*
 * definitions.h
 *
 *  Created on: Nov 16, 2015
 *      Author: nico
 */

#ifndef INCLUDE_DEFINITIONS_H_
#define INCLUDE_DEFINITIONS_H_

#include <string>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <istream>
#include <sys/stat.h>
#include <set>
#include "stdint.h"
#include <sstream>
#include <algorithm>
#include <sys/resource.h>
#include <sys/stat.h>
#include <fstream>
#include <assert.h>

using namespace std;
using namespace sdsl;

namespace lzrlbwt{

typedef uint64_t ulint;
typedef long int lint;
typedef unsigned int uint;
typedef unsigned short int t_char; ///< Type for char conversion
typedef unsigned short int t_errors; ///< Type for ERRORS

typedef unsigned char uchar;
typedef unsigned char symbol;
typedef unsigned char uint8;

typedef pair<ulint,ulint> range_t;

}

#endif /* INCLUDE_DEFINITIONS_H_ */
