//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
// Non--inline functions for BellBiharmonic elements
#include "foeppl_von_karman.h"

namespace oomph
{
//======================================================================
/// Set the data for the number of Variables at each node
//======================================================================
 template<>
 const unsigned FoepplVonKarmanC1CurvedBellElement<2,2,3>::Initial_Nvalue[3] = {8,8,8};

 template<>
 const unsigned FoepplVonKarmanC1CurvedBellElement<2,3,3>::Initial_Nvalue[6] = {8,8,8,2,2,2};
 
 template<>
 const unsigned FoepplVonKarmanC1CurvedBellElement<2,4,3>::Initial_Nvalue[10]= {8,8,8,2,2,2,2,2,2,2};

 template<>
 const unsigned FoepplVonKarmanC1CurvedBellElement<2,2,5>::Initial_Nvalue[3] = {8,8,8};

 template<>
 const unsigned FoepplVonKarmanC1CurvedBellElement<2,3,5>::Initial_Nvalue[6] = {8,8,8,2,2,2};
 
 template<>
 const unsigned FoepplVonKarmanC1CurvedBellElement<2,4,5>::Initial_Nvalue[10]= {8,8,8,2,2,2,2,2,2,2};

//=======================================================================
/// Shape function for specific TElement<DIM,NNODE,BOUNDARY_ORDER>
//=======================================================================
 template<unsigned DIM, unsigned NNODE, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE,BOUNDARY_ORDER>::shape_u(const Vector<double> &s, Shape &psi) const
   {
    // Use the base TElement version of shape
    this->shape(s,psi);
   }
//=======================================================================
/// Derivatives of shape functions for specific TElement<2,2,BOUNDARY_ORDER>
//=======================================================================
 template<unsigned DIM, unsigned NNODE, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE,BOUNDARY_ORDER>::dshape_u_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const
   {
    // Use the base TElement version of dshape_local
    this->dshape_local(s,psi,dpsids);
   }
//====================================================================
// Force build of templates
//====================================================================
template class FoepplVonKarmanC1CurvedBellElement<2,2,3>;

template class FoepplVonKarmanC1CurvedBellElement<2,3,3>;

template class FoepplVonKarmanC1CurvedBellElement<2,4,3>;

template class FoepplVonKarmanC1CurvedBellElement<2,2,5>;

template class FoepplVonKarmanC1CurvedBellElement<2,3,5>;

template class FoepplVonKarmanC1CurvedBellElement<2,4,5>;
}
