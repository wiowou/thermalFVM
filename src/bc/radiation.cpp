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

#include "radiation.h"
#include <cmath>

namespace math{
namespace model{
namespace fvm{
namespace att{

const double Radiation::stephBoltz = 5.67e-8; // W/m2-K4

}/*att*/ }/*fvm*/ }/*model*/ }/*math*/

namespace math{
namespace model{
namespace fvm{
namespace bc{

void Radiation::Apply( Face* f1 )
{
  if (pID != f1->pID ) return;
  Domain* domain = static_cast<Domain*>(dom);
  ThermalDOF* dof = static_cast<ThermalDOF*>(dom);
  unsigned long long fidx = domain->Index(f1);
  unsigned long long cidx = domain->Index(f1->cell[0]);
  
  double tsurf = dof->ts[cidx]; //get surface temp from face's cell
  double tsurr;
  if (TSurrCell == nullptr) 
  {
    tsurr = Tsurr(f1->ctr[0], f1->ctr[1], f1->ctr[2], time);
  }
  else
  {
    unsigned long long bcidx = domain->Index(TSurrCell);
    tsurr = dof->ts[bcidx];
  }
  att::Radiation* radAttr = static_cast<att::Radiation*>(f1->att);
  double hA = f1->area * radAttr->emis * radAttr->shapeFact * att::Radiation::stephBoltz  
    * (pow(tsurr, 2.0) + pow(tsurf, 2.0) ) * (tsurr + tsurf);
  dof->Q_a[cidx] -= hA;
  dof->Q_b[cidx] += hA * tsurr;
}

}/*bc*/ }/*fvm*/ }/*model*/ }/*math*/ 

