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
// template<unsigned DIM, unsigned NNODE_1D>
// const unsigned FoepplVonKarmanBellElement<DIM,NNODE_1D>::Initial_Nvalue[1] = 8;

 template<>
 const unsigned FoepplVonKarmanBellElement<2,2>::Initial_Nvalue[3] = {8,8,8};

 template<>
 const unsigned FoepplVonKarmanBellElement<2,3>::Initial_Nvalue[6] = {8,8,8,2,2,2};
 
 template<>
 const unsigned FoepplVonKarmanBellElement<2,4>::Initial_Nvalue[10]= {8,8,8,2,2,2,2,2,2,2};

//  template<>
// const unsigned FoepplVonKarmanBellElement<2,4>::Initial_Nvalue[9] = {8,8,8,2,2,2,2,2,2,};

//=======================================================================
/// Shape function for specific TElement<2,2>
//=======================================================================
 template<>
 void FoepplVonKarmanBellElement<2,2>::shape_u(const Vector<double> &s, Shape &psi) const
   {
    psi[0] = s[0];
    psi[1] = s[1];
    psi[2] = 1.0-s[0]-s[1];
   }
//=======================================================================
/// Derivatives of shape functions for specific TElement<2,2>
//=======================================================================
 template<>
 void FoepplVonKarmanBellElement<2,2>::dshape_u_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const
   {
    shape_u(s, psi);
    
    // Derivatives
    dpsids(0,0) = 1.0;
    dpsids(0,1) = 0.0;
    dpsids(1,0) = 0.0;
    dpsids(1,1) = 1.0;
    dpsids(2,0) = -1.0;
    dpsids(2,1) = -1.0;
   }

//=======================================================================
/// Shape function for specific TElement<2,3>
//=======================================================================
 template<>
 void FoepplVonKarmanBellElement<2,3>::shape_u(const Vector<double> &s, Shape &psi) const
{
 // Reconstruct the third area coordinate
 double s_2=1.0-s[0]-s[1];

 // note that s[2] needs replacing by s_2 everywhere since only 
 // two independent variables s[0],s[1] and s_2 is expressed in terms of those
 // later.
 psi[0] = 2.0*s[0]*(s[0]-0.5);
 psi[1] = 2.0*s[1]*(s[1]-0.5);
 psi[2] = 2.0*s_2 *(s_2 -0.5);
 psi[3] = 4.0*s[0]*s[1];
 psi[4] = 4.0*s[1]*s_2;
 psi[5] = 4.0*s_2*s[0];
}


//=======================================================================
/// Derivatives of shape functions for specific TElement<2,3>
//=======================================================================
 template<>
 void FoepplVonKarmanBellElement<2,3>::dshape_u_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const
   {
 //ALH: Don't know why object qualifier is needed
 shape_u(s, psi);

 dpsids(0,0) = 4.0*s[0]-1.0;
 dpsids(0,1) = 0.0;
 dpsids(1,0) = 0.0;
 dpsids(1,1) = 4.0*s[1]-1.0;
 dpsids(2,0) = 2.0*(2.0*s[0]-1.5+2.0*s[1]);
 dpsids(2,1) = 2.0*(2.0*s[0]-1.5+2.0*s[1]);
 dpsids(3,0) = 4.0*s[1];
 dpsids(3,1) = 4.0*s[0];
 dpsids(4,0) = -4.0*s[1];
 dpsids(4,1) = 4.0*(1.0-s[0]-2.0*s[1]);
 dpsids(5,0) = 4.0*(1.0-2.0*s[0]-s[1]);
 dpsids(5,1) = -4.0*s[0];
}

//=======================================================================
/// Shape function for specific TElement<2,4>
//=======================================================================
template<>
void FoepplVonKarmanBellElement<2,4>::shape_u(const Vector<double> &s, Shape &psi) const
{
 psi[0] = 0.5*s[0]*(3.0*s[0]-2.0)*(3.0*s[0]-1.0);
 psi[1] = 0.5*s[1]*(3.0*s[1]-2.0)*(3.0*s[1]-1.0);
 psi[2] = 0.5*(1.0-s[0]-s[1])*(1.0-3.0*s[0]-3.0*s[1])*(2.0-3.0*s[0]-3.0*s[1]);
 psi[3] = 4.5*s[0]*s[1]*(3.0*s[0]-1.0);
 psi[4] = 4.5*s[0]*s[1]*(3.0*s[1]-1.0);
 psi[5] = 4.5*s[1]*(1.0-s[0]-s[1])*(3.0*s[1]-1.0);
 psi[6] = 4.5*s[1]*(1.0-s[0]-s[1])*(3.0*(1.0-s[0]-s[1])-1.0);
 psi[7] = 4.5*s[0]*(1.0-s[0]-s[1])*(2.0-3*s[0]-3*s[1]);
 psi[8] = 4.5*s[0]*(1.0-s[0]-s[1])*(3.0*s[0]-1.0);
 psi[9] = 27.0*s[0]*s[1]*(1.0-s[0]-s[1]);
}

//=======================================================================
/// Derivatives of shape functions for specific TElement<2,4>
//=======================================================================
template<>
void FoepplVonKarmanBellElement<2,4>::dshape_u_local(const Vector<double> &s,
                   Shape &psi, DShape &dpsids) const
{
  
 //ALH: Don't know why object qualifier is needed
 shape_u(s, psi);
 
 dpsids(0,0) = 13.5*s[0]*s[0]-9.0*s[0]+1.0;
 dpsids(0,1) = 0.0;
 dpsids(1,0) = 0.0;
 dpsids(1,1) = 13.5*s[1]*s[1]-9.0*s[1]+1.0;
 dpsids(2,0) = 0.5*(36.0*s[0]+36.0*s[1]-27.0*s[0]*s[0]-
                    27.0*s[1]*s[1]-54.0*s[0]*s[1]-11.0);
 dpsids(2,1) = 0.5*(36.0*s[0]+36.0*s[1]-27.0*s[0]*s[0]-
                    27.0*s[1]*s[1]-54.0*s[0]*s[1]-11.0);
 dpsids(3,0) = 27.0*s[0]*s[1]-4.5*s[1];
 dpsids(3,1) = 4.5*s[0]*(3.0*s[0]-1.0);
 dpsids(4,0) = 4.5*s[1]*(3.0*s[1]-1.0);
 dpsids(4,1) = 27.0*s[0]*s[1]-4.5*s[0];
 dpsids(5,0) = 4.5*(s[1]-3.0*s[1]*s[1]);
 dpsids(5,1) = 4.5*(s[0]-6.0*s[0]*s[1]-9.0*s[1]*s[1]+8*s[1]-1.0);
 dpsids(6,0) = 4.5*(6.0*s[0]*s[1]-5.0*s[1]+6.0*s[1]*s[1]);
 dpsids(6,1) = 4.5*(2.0-5.0*s[0]+3.0*s[0]*s[0]+12.0*s[0]*s[1]-
                    10.0*s[1]+9.0*s[1]*s[1]);
 dpsids(7,0) = 4.5*(2.0-10.0*s[0]+9.0*s[0]*s[0]+12.0*s[0]*s[1]-
                    5.0*s[1]+3.0*s[1]*s[1]);
 dpsids(7,1) = 4.5*(6.0*s[0]*s[0]-5.0*s[0]+6.0*s[0]*s[1]);
 dpsids(8,0) = 4.5*(s[1]-6.0*s[0]*s[1]-9.0*s[0]*s[0]+8*s[0]-1.0);
 dpsids(8,1) = 4.5*(s[0]-3.0*s[0]*s[0]);
 dpsids(9,0) = 27.0*s[1]-54.0*s[0]*s[1]-27.0*s[1]*s[1];
 dpsids(9,1) = 27.0*s[0]-54.0*s[0]*s[1]-27.0*s[0]*s[0];
 
}
//====================================================================
// Force build of templates
//====================================================================
template class FoepplVonKarmanBellElement<2,2>;

template class FoepplVonKarmanBellElement<2,3>;

template class FoepplVonKarmanBellElement<2,4>;
}
