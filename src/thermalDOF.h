//<license>
/*
   This file is part of Behmeth.
   Author: Behram Kapadia, wiowou@hotmail.com

    Behmeth is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Behmeth is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Behmeth.  If not, see <http://www.gnu.org/licenses/>.
*/

//</license>

#ifndef _MATH_MODEL_FVM_THERMALDOF_h
#define _MATH_MODEL_FVM_THERMALDOF_h

#include "matrix/vect.h"
#include "matrix/matrix.h"

namespace math{
namespace model{
namespace fvm{

#ifdef MYDEBUG
  class ThermalDOFTest;
#endif //MYDEBUG
class ThermalDOF
{

public:
  //cell DOFs
  Vect<double> ts;
  //face DOFs
  // y = ax + b
  Vect<double> Q;
  Vect<double> Q_a;
  Vect<double> Q_b;
  Matrix<SparseVect> A_ts;
  
protected:

private:


#ifdef MYDEBUG
  friend class ThermalDOFTest;
#endif //MYDEBUG
};

}/*fvm*/ }/*model*/ }/*math*/ 

#endif /*_MATH_MODEL_FVM_THERMALDOF_h */
