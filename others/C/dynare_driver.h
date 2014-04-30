/*
 * Copyright (C) 2011-2012 Houtan Bastani, Daniel Waggoner, Tao Zha
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * More information is available at <http://www.gnu.org/licenses/>.
 */

#ifndef _DYNARE_C_DRIVER_H
#define _DYNARE_C_DRIVER_H

#include <cstdlib>
#include <vector>
#include <string>
#include <map>
#include <limits>

using namespace std;

struct aux_vars_t {
  int endo_index, type, orig_index, orig_lead_lag;
} ;

#endif // ! _DYNARE_C_DRIVER_H
