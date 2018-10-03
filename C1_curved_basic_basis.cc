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
//oomph-lib headers
#include "C1_curved_elements.h"

// The matrix that connverts from the set of monomials to a the shape functions
// defined on the basic element
// this output is directly from the mathematica script that generated the basic
// shape functions.

// We could possibly change this out for a static DenseMatrix which has the LU
// decomposition stored - might give smaller round off and be of comparably
// speed HERE

// We could definitely store this as a static 
namespace oomph {
namespace MyC1CurvedElements {

/// Precomputed basis polynomials for the 3rd order boundary representation
template <>
void  BernadouElementBasis<3>::full_basic_polynomials(const Vector<double>& s, Shape&
phi) const
 {
  // This will be replaced eventually - in favour of the explicit shape
  // functions
  // Get the number of basic basis function
  const unsigned nbasic=n_basic_basis_functions();
  // Initialise the p7 shape functions
  Shape p7(nbasic);
  full_basis_monomials(s,p7);
  // Initialise the basic association matrix
  DenseMatrix<double> a_matrix (nbasic,nbasic,0.0);
  monomial_to_basic_matrix(a_matrix);
  // Loop over the 36 Basis functions
  for(unsigned i=0;i<nbasic; ++i)
   {
    // Zero the shape functions
    phi(i) =0.0;
    // Zero nodal dofs
    for(unsigned j=0;j<nbasic; ++j)
     {phi(i)+=a_matrix(i,j)*p7(j);}
   }
 }

/// Precomputed dbasis polynomials for the 3rd order boundary representation
template <>
void  BernadouElementBasis<3>::dfull_basic_polynomials(const Vector<double>& s, DShape&
dphi) const
 {
  // This will be replaced eventually - in favour of the explicit shape
  // functions
  // Get the number of basic basis function
  const unsigned nbasic=n_basic_basis_functions();
  // Initialise the p7 shape functions
  DShape dp7(nbasic,2);
  dfull_basis_monomials(s,dp7);
  // Initialise the basic association matrix
  DenseMatrix<double> a_matrix (nbasic,nbasic,0.0);
  monomial_to_basic_matrix(a_matrix);
  // Loop over the 36 Basis functions
  for(unsigned i=0;i<nbasic; ++i)
   {
    for(unsigned alpha=0; alpha<2;++alpha)
     {
      // Zero the shape functions
      dphi(i,alpha) =0.0;
      // Zero nodal dofs
      for(unsigned j=0;j<nbasic; ++j)
       {dphi(i,alpha)+=a_matrix(i,j)*dp7(j,alpha);}
     }
   }
 }

/// Precomputed d2basis polynomials for the 3rd order boundary representation
template <>
void  BernadouElementBasis<3>::d2full_basic_polynomials(const Vector<double>& s, DShape&
d2phi) const
 {
  // This will be replaced eventually - in favour of the explicit shape
  // functions
  // Get the number of basic basis function
  const unsigned nbasic=n_basic_basis_functions();
  // Initialise the p7 shape functions
  DShape d2p7(nbasic,3);
  d2full_basis_monomials(s,d2p7);
  // Initialise the basic association matrix
  DenseMatrix<double> a_matrix (nbasic,nbasic,0.0);
  monomial_to_basic_matrix(a_matrix);
  // Loop over the 36 Basis functions
  for(unsigned i=0;i<nbasic; ++i)
   {
    for(unsigned alpha=0; alpha<3;++alpha)
     {
      // Zero the shape functions
      d2phi(i,alpha) =0.0;
      // Zero nodal dofs
      for(unsigned j=0;j<nbasic; ++j)
       {d2phi(i,alpha)+=a_matrix(i,j)*d2p7(j,alpha);}
     }
   }

 }

/// Precomputed basis polynomials for the 5th order boundary representation
template <>
void  BernadouElementBasis<5>::full_basic_polynomials(const Vector<double>& s_basic,
Shape& phi) const
 {
  // For convenience 
  const double s=s_basic[0], t=s_basic[1];
  // Now fill in (automatically generated code
  phi[0] = (pow(s, 2)*(598360*s - 8130732*pow(s, 2) + 40143300*pow(s, 3) -
94869252*pow(s, 4) + 116223876*pow(s, 5) - 71313696*pow(s, 6) + 17352144*pow(s,
7) - 10246440*pow(t, 2) + 126762915*s*pow(t, 2) - 538504848*pow(s, 2)*pow(t, 2)
+ 994544973*pow(s, 3)*pow(t, 2) - 823436604*pow(s, 4)*pow(t, 2) +
250880004*pow(s, 5)*pow(t, 2) + 22167939*pow(t, 3) - 247491396*s*pow(t, 3) +
924349995*pow(s, 2)*pow(t, 3) - 1297023192*pow(s, 3)*pow(t, 3) +
592038396*pow(s, 4)*pow(t, 3) - 12586266*pow(t, 4) + 92584809*s*pow(t, 4) -
342910908*pow(s, 2)*pow(t, 4) + 327894480*pow(s, 3)*pow(t, 4) + 5508243*pow(t,
5) + 42791328*s*pow(t, 5) - 33765336*pow(s, 2)*pow(t, 5) - 9799704*pow(t, 6) -
15545196*s*pow(t, 6) + 4956228*pow(t, 7)))/4000.;

  phi[1] = (pow(t, 2)*(598360*t - 10246440*pow(s, 2) + 126762915*t*pow(s, 2) +
22167939*pow(s, 3) - 247491396*t*pow(s, 3) - 12586266*pow(s, 4) +
92584809*t*pow(s, 4) + 5508243*pow(s, 5) + 42791328*t*pow(s, 5) - 9799704*pow(s,
6) - 15545196*t*pow(s, 6) + 4956228*pow(s, 7) - 8130732*pow(t, 2) -
538504848*pow(s, 2)*pow(t, 2) + 924349995*pow(s, 3)*pow(t, 2) - 342910908*pow(s,
4)*pow(t, 2) - 33765336*pow(s, 5)*pow(t, 2) + 40143300*pow(t, 3) +
994544973*pow(s, 2)*pow(t, 3) - 1297023192*pow(s, 3)*pow(t, 3) +
327894480*pow(s, 4)*pow(t, 3) - 94869252*pow(t, 4) - 823436604*pow(s, 2)*pow(t,
4) + 592038396*pow(s, 3)*pow(t, 4) + 116223876*pow(t, 5) + 250880004*pow(s,
2)*pow(t, 5) - 71313696*pow(t, 6) + 17352144*pow(t, 7)))/4000.;

  phi[2] = -((-2000 - 4000*s - 4000*t - 12000*s*t - 6000*pow(s, 2) -
24000*t*pow(s, 2) + 1978086*pow(s, 3) + 3932172*t*pow(s, 3) - 11934864*pow(s, 4)
- 19937556*t*pow(s, 4) + 26368362*pow(s, 5) + 32799168*t*pow(s, 5) -
25075656*pow(s, 6) - 17352144*t*pow(s, 6) + 8676072*pow(s, 7) - 6000*pow(t, 2) -
24000*s*pow(t, 2) - 16927218*pow(s, 2)*pow(t, 2) + 53140455*pow(s, 3)*pow(t, 2)
- 15361002*pow(s, 4)*pow(t, 2) - 31939920*pow(s, 5)*pow(t, 2) + 1978086*pow(t,
3) + 3932172*s*pow(t, 3) + 53140455*pow(s, 2)*pow(t, 3) - 197161452*pow(s,
3)*pow(t, 3) + 117457452*pow(s, 4)*pow(t, 3) - 11934864*pow(t, 4) -
19937556*s*pow(t, 4) - 15361002*pow(s, 2)*pow(t, 4) + 117457452*pow(s, 3)*pow(t,
4) + 26368362*pow(t, 5) + 32799168*s*pow(t, 5) - 31939920*pow(s, 2)*pow(t, 5) -
25075656*pow(t, 6) - 17352144*s*pow(t, 6) + 8676072*pow(t, 7))*pow(-1 + s + t,
2))/2000.;

  phi[3] = -(pow(s, 2)*(6280*s - 85636*pow(s, 2) + 425100*pow(s, 3) -
1011996*pow(s, 4) + 1250748*pow(s, 5) - 775008*pow(s, 6) + 190512*pow(s, 7) -
103738*pow(t, 2) + 1284763*s*pow(t, 2) - 5469804*pow(s, 2)*pow(t, 2) +
10136979*pow(s, 3)*pow(t, 2) - 8434692*pow(s, 4)*pow(t, 2) + 2586492*pow(s,
5)*pow(t, 2) + 222565*pow(t, 3) - 2481858*s*pow(t, 3) + 9277785*pow(s, 2)*pow(t,
3) - 13035816*pow(s, 3)*pow(t, 3) + 5967108*pow(s, 4)*pow(t, 3) - 127188*pow(t,
4) + 919827*s*pow(t, 4) - 3415284*pow(s, 2)*pow(t, 4) + 3259440*pow(s, 3)*pow(t,
4) + 59229*pow(t, 5) + 429624*s*pow(t, 5) - 331128*pow(s, 2)*pow(t, 5) -
102384*pow(t, 6) - 158436*s*pow(t, 6) + 51516*pow(t, 7)))/400.;

  phi[4] = -(t*pow(s, 2)*(6680 - 97796*s - 9488*t + 115031*s*t + 536936*pow(s,
2) - 418950*t*pow(s, 2) - 1445940*pow(s, 3) + 687627*t*pow(s, 3) +
2043000*pow(s, 4) - 601344*t*pow(s, 4) - 1450224*pow(s, 5) + 227124*t*pow(s, 5)
+ 406944*pow(s, 6) + 55323*pow(t, 2) - 638370*s*pow(t, 2) + 1640007*pow(s,
2)*pow(t, 2) - 1194588*pow(s, 3)*pow(t, 2) + 201204*pow(s, 4)*pow(t, 2) -
122310*pow(t, 3) + 1602207*s*pow(t, 3) - 3564648*pow(s, 2)*pow(t, 3) +
1718496*pow(s, 3)*pow(t, 3) + 18927*pow(t, 4) - 1133028*s*pow(t, 4) +
1737936*pow(s, 2)*pow(t, 4) + 102384*pow(t, 5) + 158436*s*pow(t, 5) -
51516*pow(t, 6)))/400.;

  phi[5] = (s*pow(t, 2)*(-6680 + 9488*s + 97796*t - 115031*s*t - 55323*pow(s, 2)
+ 638370*t*pow(s, 2) + 122310*pow(s, 3) - 1602207*t*pow(s, 3) - 18927*pow(s, 4)
+ 1133028*t*pow(s, 4) - 102384*pow(s, 5) - 158436*t*pow(s, 5) + 51516*pow(s, 6)
- 536936*pow(t, 2) + 418950*s*pow(t, 2) - 1640007*pow(s, 2)*pow(t, 2) +
3564648*pow(s, 3)*pow(t, 2) - 1737936*pow(s, 4)*pow(t, 2) + 1445940*pow(t, 3) -
687627*s*pow(t, 3) + 1194588*pow(s, 2)*pow(t, 3) - 1718496*pow(s, 3)*pow(t, 3) -
2043000*pow(t, 4) + 601344*s*pow(t, 4) - 201204*pow(s, 2)*pow(t, 4) +
1450224*pow(t, 5) - 227124*s*pow(t, 5) - 406944*pow(t, 6)))/400.;

  phi[6] = -(pow(t, 2)*(6280*t - 103738*pow(s, 2) + 1284763*t*pow(s, 2) +
222565*pow(s, 3) - 2481858*t*pow(s, 3) - 127188*pow(s, 4) + 919827*t*pow(s, 4) +
59229*pow(s, 5) + 429624*t*pow(s, 5) - 102384*pow(s, 6) - 158436*t*pow(s, 6) +
51516*pow(s, 7) - 85636*pow(t, 2) - 5469804*pow(s, 2)*pow(t, 2) + 9277785*pow(s,
3)*pow(t, 2) - 3415284*pow(s, 4)*pow(t, 2) - 331128*pow(s, 5)*pow(t, 2) +
425100*pow(t, 3) + 10136979*pow(s, 2)*pow(t, 3) - 13035816*pow(s, 3)*pow(t, 3) +
3259440*pow(s, 4)*pow(t, 3) - 1011996*pow(t, 4) - 8434692*pow(s, 2)*pow(t, 4) +
5967108*pow(s, 3)*pow(t, 4) + 1250748*pow(t, 5) + 2586492*pow(s, 2)*pow(t, 5) -
775008*pow(t, 6) + 190512*pow(t, 7)))/400.;

  phi[7] = -(s*(-200 - 400*s - 400*t - 1200*s*t + 26178*pow(s, 2) +
51156*t*pow(s, 2) - 142272*pow(s, 3) - 233388*t*pow(s, 3) + 300726*pow(s, 4) +
368064*t*pow(s, 4) - 279288*pow(s, 5) - 190512*t*pow(s, 5) + 95256*pow(s, 6) +
29518*pow(t, 2) - 218907*s*pow(t, 2) + 355743*pow(s, 2)*pow(t, 2) +
278154*pow(s, 3)*pow(t, 2) - 563760*pow(s, 4)*pow(t, 2) - 181350*pow(t, 3) +
1104822*s*pow(t, 3) - 1916298*pow(s, 2)*pow(t, 3) + 722196*pow(s, 3)*pow(t, 3) +
448020*pow(t, 4) - 1782000*s*pow(t, 4) + 1684800*pow(s, 2)*pow(t, 4) -
495720*pow(t, 5) + 939600*s*pow(t, 5) + 203472*pow(t, 6))*pow(-1 + s + t,
2))/200.;

  phi[8] = -(t*(-200 - 400*s - 400*t - 1200*s*t + 29518*pow(s, 2) -
218907*t*pow(s, 2) - 181350*pow(s, 3) + 1104822*t*pow(s, 3) + 448020*pow(s, 4) -
1782000*t*pow(s, 4) - 495720*pow(s, 5) + 939600*t*pow(s, 5) + 203472*pow(s, 6) +
26178*pow(t, 2) + 51156*s*pow(t, 2) + 355743*pow(s, 2)*pow(t, 2) -
1916298*pow(s, 3)*pow(t, 2) + 1684800*pow(s, 4)*pow(t, 2) - 142272*pow(t, 3) -
233388*s*pow(t, 3) + 278154*pow(s, 2)*pow(t, 3) + 722196*pow(s, 3)*pow(t, 3) +
300726*pow(t, 4) + 368064*s*pow(t, 4) - 563760*pow(s, 2)*pow(t, 4) -
279288*pow(t, 5) - 190512*s*pow(t, 5) + 95256*pow(t, 6))*pow(-1 + s + t,
2))/200.;

  phi[9] = (pow(s, 2)*(40*s - 548*pow(s, 2) + 2740*pow(s, 3) - 6588*pow(s, 4) +
8244*pow(s, 5) - 5184*pow(s, 6) + 1296*pow(s, 7) - 628*pow(t, 2) + 7783*s*pow(t,
2) - 33192*pow(s, 2)*pow(t, 2) + 61677*pow(s, 3)*pow(t, 2) - 51516*pow(s,
4)*pow(t, 2) + 15876*pow(s, 5)*pow(t, 2) + 1339*pow(t, 3) - 14904*s*pow(t, 3) +
55755*pow(s, 2)*pow(t, 3) - 78408*pow(s, 3)*pow(t, 3) + 35964*pow(s, 4)*pow(t,
3) - 774*pow(t, 4) + 5481*s*pow(t, 4) - 20412*pow(s, 2)*pow(t, 4) + 19440*pow(s,
3)*pow(t, 4) + 387*pow(t, 5) + 2592*s*pow(t, 5) - 1944*pow(s, 2)*pow(t, 5) -
648*pow(t, 6) - 972*s*pow(t, 6) + 324*pow(t, 7)))/80.;

  phi[10] = (t*pow(s, 2)*(40 - 588*s - 50*t + 599*s*t + 3248*pow(s, 2) -
2070*t*pow(s, 2) - 8820*pow(s, 3) + 3141*t*pow(s, 3) + 12600*pow(s, 4) -
2592*t*pow(s, 4) - 9072*pow(s, 5) + 972*t*pow(s, 5) + 2592*pow(s, 6) +
325*pow(t, 2) - 3780*s*pow(t, 2) + 9621*pow(s, 2)*pow(t, 2) - 6804*pow(s,
3)*pow(t, 2) + 972*pow(s, 4)*pow(t, 2) - 720*pow(t, 3) + 9621*s*pow(t, 3) -
21384*pow(s, 2)*pow(t, 3) + 10368*pow(s, 3)*pow(t, 3) + 81*pow(t, 4) -
6804*s*pow(t, 4) + 10368*pow(s, 2)*pow(t, 4) + 648*pow(t, 5) + 972*s*pow(t, 5) -
324*pow(t, 6)))/40.;

  phi[11] = (pow(s, 2)*pow(t, 2)*(10 - 51*s - 325*t + 3780*s*t - 630*pow(s, 2) -
9621*t*pow(s, 2) + 2979*pow(s, 3) + 6804*t*pow(s, 3) - 3888*pow(s, 4) -
972*t*pow(s, 4) + 1620*pow(s, 5) + 720*pow(t, 2) - 9621*s*pow(t, 2) +
21384*pow(s, 2)*pow(t, 2) - 10368*pow(s, 3)*pow(t, 2) - 81*pow(t, 3) +
6804*s*pow(t, 3) - 10368*pow(s, 2)*pow(t, 3) - 648*pow(t, 4) - 972*s*pow(t, 4) +
324*pow(t, 5)))/80.;

  phi[12] = (pow(s, 2)*pow(t, 2)*(10 - 325*s - 51*t + 3780*s*t + 720*pow(s, 2) -
9621*t*pow(s, 2) - 81*pow(s, 3) + 6804*t*pow(s, 3) - 648*pow(s, 4) -
972*t*pow(s, 4) + 324*pow(s, 5) - 630*pow(t, 2) - 9621*s*pow(t, 2) +
21384*pow(s, 2)*pow(t, 2) - 10368*pow(s, 3)*pow(t, 2) + 2979*pow(t, 3) +
6804*s*pow(t, 3) - 10368*pow(s, 2)*pow(t, 3) - 3888*pow(t, 4) - 972*s*pow(t, 4)
+ 1620*pow(t, 5)))/80.;

  phi[13] = -(s*pow(t, 2)*(-40 + 50*s + 588*t - 599*s*t - 325*pow(s, 2) +
3780*t*pow(s, 2) + 720*pow(s, 3) - 9621*t*pow(s, 3) - 81*pow(s, 4) +
6804*t*pow(s, 4) - 648*pow(s, 5) - 972*t*pow(s, 5) + 324*pow(s, 6) - 3248*pow(t,
2) + 2070*s*pow(t, 2) - 9621*pow(s, 2)*pow(t, 2) + 21384*pow(s, 3)*pow(t, 2) -
10368*pow(s, 4)*pow(t, 2) + 8820*pow(t, 3) - 3141*s*pow(t, 3) + 6804*pow(s,
2)*pow(t, 3) - 10368*pow(s, 3)*pow(t, 3) - 12600*pow(t, 4) + 2592*s*pow(t, 4) -
972*pow(s, 2)*pow(t, 4) + 9072*pow(t, 5) - 972*s*pow(t, 5) - 2592*pow(t,
6)))/40.;

  phi[14] = (pow(t, 2)*(40*t - 628*pow(s, 2) + 7783*t*pow(s, 2) + 1339*pow(s, 3)
- 14904*t*pow(s, 3) - 774*pow(s, 4) + 5481*t*pow(s, 4) + 387*pow(s, 5) +
2592*t*pow(s, 5) - 648*pow(s, 6) - 972*t*pow(s, 6) + 324*pow(s, 7) - 548*pow(t,
2) - 33192*pow(s, 2)*pow(t, 2) + 55755*pow(s, 3)*pow(t, 2) - 20412*pow(s,
4)*pow(t, 2) - 1944*pow(s, 5)*pow(t, 2) + 2740*pow(t, 3) + 61677*pow(s,
2)*pow(t, 3) - 78408*pow(s, 3)*pow(t, 3) + 19440*pow(s, 4)*pow(t, 3) -
6588*pow(t, 4) - 51516*pow(s, 2)*pow(t, 4) + 35964*pow(s, 3)*pow(t, 4) +
8244*pow(t, 5) + 15876*pow(s, 2)*pow(t, 5) - 5184*pow(t, 6) + 1296*pow(t,
7)))/80.;

  phi[15] = -((-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*(-20 + 74*s - 40*t + 108*s*t -
90*pow(s, 2) - 72*t*pow(s, 2) + 36*pow(s, 3) + 279*pow(t, 2) - 360*s*pow(t, 2) -
234*pow(t, 3))*pow(-1 + s + t, 2))/40.;

  phi[16] = -(s*t*(-1 + 2*s + 2*t)*(-2 + 3*s + 3*t)*(-1 + 3*s + 3*t)*(-5 + 6*s +
6*t)*(-1 + 6*s + 6*t)*pow(-1 + s + t, 2))/10.;

  phi[17] = -((-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(-20 - 40*s + 74*t + 108*s*t +
279*pow(s, 2) - 360*t*pow(s, 2) - 234*pow(s, 3) - 90*pow(t, 2) - 72*s*pow(t, 2)
+ 36*pow(t, 3))*pow(-1 + s + t, 2))/40.;

  phi[18] = -16*s*(-1 + 3*t)*(-2 + 3*s + 3*t)*(-1 + 6*t)*(-5 + 6*s + 6*t)*pow(t,
2)*pow(-1 + s + t, 2);

  phi[19] = -16*(-1 + 3*s)*(-1 + 6*s)*t*(-2 + 3*s + 3*t)*(-5 + 6*s + 6*t)*pow(s,
2)*pow(-1 + s + t, 2);

  phi[20] = 8*(-1 + 3*s)*(-1 + 6*s)*(-1 + s + t)*(-1 + 3*t)*(-1 + 6*t)*pow(s,
2)*pow(t, 2)*sqrt(2);

  phi[21] = (-1944*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(16*t + 32*s*t - 279*pow(s,
2) + 360*t*pow(s, 2) + 234*pow(s, 3) - 40*pow(t, 2) - 48*s*pow(t, 2) + 24*pow(t,
3))*pow(-1 + s + t, 2))/125.;

  phi[22] = (243*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(10*t + 20*s*t - 183*pow(s, 2)
+ 168*t*pow(s, 2) + 234*pow(s, 3) - 22*pow(t, 2) - 24*s*pow(t, 2) + 12*pow(t,
3))*pow(-1 + s + t, 2))/16.;

  phi[23] = (-81*(-1 + 6*t)*pow(t, 2)*(-60*t - 120*s*t + 1003*pow(s, 2) -
1809*t*pow(s, 2) - 2898*pow(s, 3) + 3078*t*pow(s, 3) + 2016*pow(s, 4) +
222*pow(t, 2) + 324*s*pow(t, 2) + 648*pow(s, 2)*pow(t, 2) - 270*pow(t, 3) -
216*s*pow(t, 3) + 108*pow(t, 4))*pow(-1 + s + t, 2))/16.;

  phi[24] = (324*pow(t, 2)*(480*t + 960*s*t - 5929*pow(s, 2) + 22458*t*pow(s, 2)
+ 25473*pow(s, 3) - 77868*t*pow(s, 3) - 35388*pow(s, 4) + 57888*t*pow(s, 4) +
16092*pow(s, 5) - 3216*pow(t, 2) - 5472*s*pow(t, 2) - 17172*pow(s, 2)*pow(t, 2)
+ 48600*pow(s, 3)*pow(t, 2) + 7488*pow(t, 3) + 9504*s*pow(t, 3) - 2592*pow(s,
2)*pow(t, 3) - 7344*pow(t, 4) - 5184*s*pow(t, 4) + 2592*pow(t, 5))*pow(-1 + s +
t, 2))/125.;

  phi[25] = (324*pow(s, 2)*(480*s + 960*s*t - 3216*pow(s, 2) - 5472*t*pow(s, 2)
+ 7488*pow(s, 3) + 9504*t*pow(s, 3) - 7344*pow(s, 4) - 5184*t*pow(s, 4) +
2592*pow(s, 5) - 5929*pow(t, 2) + 22458*s*pow(t, 2) - 17172*pow(s, 2)*pow(t, 2)
- 2592*pow(s, 3)*pow(t, 2) + 25473*pow(t, 3) - 77868*s*pow(t, 3) + 48600*pow(s,
2)*pow(t, 3) - 35388*pow(t, 4) + 57888*s*pow(t, 4) + 16092*pow(t, 5))*pow(-1 + s
+ t, 2))/125.;

  phi[26] = (-81*(-1 + 6*s)*pow(s, 2)*(-60*s - 120*s*t + 222*pow(s, 2) +
324*t*pow(s, 2) - 270*pow(s, 3) - 216*t*pow(s, 3) + 108*pow(s, 4) + 1003*pow(t,
2) - 1809*s*pow(t, 2) + 648*pow(s, 2)*pow(t, 2) - 2898*pow(t, 3) + 3078*s*pow(t,
3) + 2016*pow(t, 4))*pow(-1 + s + t, 2))/16.;

  phi[27] = (243*(-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*(10*s + 20*s*t - 22*pow(s, 2) -
24*t*pow(s, 2) + 12*pow(s, 3) - 183*pow(t, 2) + 168*s*pow(t, 2) + 234*pow(t,
3))*pow(-1 + s + t, 2))/16.;

  phi[28] = (-1944*(-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*(16*s + 32*s*t - 40*pow(s, 2)
- 48*t*pow(s, 2) + 24*pow(s, 3) - 279*pow(t, 2) + 360*s*pow(t, 2) + 234*pow(t,
3))*pow(-1 + s + t, 2))/125.;

  phi[29] = (-162*pow(s, 2)*pow(t, 2)*(-518 + 7483*s - 1313*t + 14000*s*t -
40700*pow(s, 2) - 21891*t*pow(s, 2) + 97059*pow(s, 3) - 20016*t*pow(s, 3) -
100800*pow(s, 4) + 29700*t*pow(s, 4) + 37476*pow(s, 5) + 4540*pow(t, 2) -
56751*s*pow(t, 2) + 126144*pow(s, 2)*pow(t, 2) - 58968*pow(s, 3)*pow(t, 2) -
1161*pow(t, 3) + 40464*s*pow(t, 3) - 63288*pow(s, 2)*pow(t, 3) - 3168*pow(t, 4)
- 5292*s*pow(t, 4) + 1620*pow(t, 5)))/125.;

  phi[30] = (-81*pow(s, 2)*pow(t, 2)*(-274 + 3167*s + 2021*t - 22688*s*t -
11236*pow(s, 2) + 72087*t*pow(s, 2) + 15867*pow(s, 3) - 82044*t*pow(s, 3) -
9144*pow(s, 4) + 30564*t*pow(s, 4) + 1620*pow(s, 5) - 3382*pow(t, 2) +
37767*s*pow(t, 2) - 94392*pow(s, 2)*pow(t, 2) + 58752*pow(s, 3)*pow(t, 2) +
1167*pow(t, 3) - 20268*s*pow(t, 3) + 33696*pow(s, 2)*pow(t, 3) + 1008*pow(t, 4)
+ 2052*s*pow(t, 4) - 540*pow(t, 5)))/32.;

  phi[31] = (81*pow(s, 2)*pow(t, 2)*(274 - 2021*s - 3167*t + 22688*s*t +
3382*pow(s, 2) - 37767*t*pow(s, 2) - 1167*pow(s, 3) + 20268*t*pow(s, 3) -
1008*pow(s, 4) - 2052*t*pow(s, 4) + 540*pow(s, 5) + 11236*pow(t, 2) -
72087*s*pow(t, 2) + 94392*pow(s, 2)*pow(t, 2) - 33696*pow(s, 3)*pow(t, 2) -
15867*pow(t, 3) + 82044*s*pow(t, 3) - 58752*pow(s, 2)*pow(t, 3) + 9144*pow(t, 4)
- 30564*s*pow(t, 4) - 1620*pow(t, 5)))/32.;

  phi[32] = (-162*pow(s, 2)*pow(t, 2)*(-518 - 1313*s + 7483*t + 14000*s*t +
4540*pow(s, 2) - 56751*t*pow(s, 2) - 1161*pow(s, 3) + 40464*t*pow(s, 3) -
3168*pow(s, 4) - 5292*t*pow(s, 4) + 1620*pow(s, 5) - 40700*pow(t, 2) -
21891*s*pow(t, 2) + 126144*pow(s, 2)*pow(t, 2) - 63288*pow(s, 3)*pow(t, 2) +
97059*pow(t, 3) - 20016*s*pow(t, 3) - 58968*pow(s, 2)*pow(t, 3) - 100800*pow(t,
4) + 29700*s*pow(t, 4) + 37476*pow(t, 5)))/125.;

  phi[33] = (-648*s*(-1 + 2*t)*(-2 + 3*t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*pow(-1
+ s + t, 2))/25.;

  phi[34] = (81*s*(-1 + 2*t)*(-1 + 3*t)*(-1 + 6*t)*(-5 + 6*s + 6*t)*pow(t,
2)*pow(-1 + s + t, 2))/4.;

  phi[35] = (81*s*(-1 + 2*s + 2*t)*(-2 + 3*s + 3*t)*(-1 + 6*t)*(-5 + 6*s +
6*t)*pow(t, 2)*pow(-1 + s + t, 2))/4.;

  phi[36] = (-648*s*(-1 + 2*s + 2*t)*(-2 + 3*s + 3*t)*(-1 + 3*s + 3*t)*(-5 + 6*s
+ 6*t)*pow(t, 2)*pow(-1 + s + t, 2))/25.;

  phi[37] = (-648*t*(-1 + 2*s + 2*t)*(-2 + 3*s + 3*t)*(-1 + 3*s + 3*t)*(-5 + 6*s
+ 6*t)*pow(s, 2)*pow(-1 + s + t, 2))/25.;

  phi[38] = (81*(-1 + 6*s)*t*(-1 + 2*s + 2*t)*(-2 + 3*s + 3*t)*(-5 + 6*s +
6*t)*pow(s, 2)*pow(-1 + s + t, 2))/4.;

  phi[39] = (81*(-1 + 2*s)*(-1 + 3*s)*(-1 + 6*s)*t*(-5 + 6*s + 6*t)*pow(s,
2)*pow(-1 + s + t, 2))/4.;

  phi[40] = (-648*(-1 + 2*s)*(-2 + 3*s)*(-1 + 3*s)*(-1 + 6*s)*t*pow(s, 2)*pow(-1
+ s + t, 2))/25.;

  phi[41] = (324*(-1 + 2*s)*(-2 + 3*s)*(-1 + 3*s)*(-1 + 6*s)*(-1 + s + t)*pow(s,
2)*pow(t, 2)*sqrt(2))/25.;

  phi[42] = (81*(-1 + 2*s)*(-1 + 3*s)*(-1 + 6*s)*(-1 + s + t)*(-1 + 6*t)*pow(s,
2)*pow(t, 2)*pow(sqrt(2), -1))/4.;

  phi[43] = (81*(-1 + 6*s)*(-1 + s + t)*(-1 + 2*t)*(-1 + 3*t)*(-1 + 6*t)*pow(s,
2)*pow(t, 2)*pow(sqrt(2), -1))/4.;

  phi[44] = (324*(-1 + s + t)*(-1 + 2*t)*(-2 + 3*t)*(-1 + 3*t)*(-1 + 6*t)*pow(s,
2)*pow(t, 2)*sqrt(2))/25.;

  phi[45] = 2916*(-1 + 2*t)*(-1 + 3*t)*(-1 + 6*t)*pow(s, 2)*pow(t, 2)*pow(-1 + s
+ t, 2);

  phi[46] = -1296*(-1 + 3*t)*(-1 + 6*t)*(-5 + 6*s + 6*t)*pow(s, 2)*pow(t,
2)*pow(-1 + s + t, 2);

  phi[47] = 1296*(-2 + 3*s + 3*t)*(-1 + 6*t)*(-5 + 6*s + 6*t)*pow(s, 2)*pow(t,
2)*pow(-1 + s + t, 2);

  phi[48] = -2916*(-1 + 2*s + 2*t)*(-2 + 3*s + 3*t)*(-5 + 6*s + 6*t)*pow(s,
2)*pow(t, 2)*pow(-1 + s + t, 2);

  phi[49] = 1296*(-1 + 6*s)*(-2 + 3*s + 3*t)*(-5 + 6*s + 6*t)*pow(s, 2)*pow(t,
2)*pow(-1 + s + t, 2);

  phi[50] = -1296*(-1 + 3*s)*(-1 + 6*s)*(-5 + 6*s + 6*t)*pow(s, 2)*pow(t,
2)*pow(-1 + s + t, 2);

  phi[51] = 2916*(-1 + 2*s)*(-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*pow(t, 2)*pow(-1 + s
+ t, 2);

  phi[52] = 1296*(-1 + 3*s)*(-1 + 6*s)*(-1 + 6*t)*pow(s, 2)*pow(t, 2)*pow(-1 + s + t, 2);

  phi[53] = 1296*(-1 + 6*s)*(-1 + 3*t)*(-1 + 6*t)*pow(s, 2)*pow(t, 2)*pow(-1 + s
+ t, 2);

  phi[54] = -729*(-1 + 6*s)*(-1 + 6*t)*(-5 + 6*s + 6*t)*pow(s, 2)*pow(t,
2)*pow(-1 + s + t, 2);

 } 

/// Precomputed dbasis polynomials for the 5th order boundary representation
template <> void  BernadouElementBasis<5>::dfull_basic_polynomials(const Vector<double>&
s_basic, DShape& dphi) const 
{ 
  // For convenience 
  const double s=s_basic[0], t=s_basic[1];
  // Now fill in (automatically generated code)
  dphi(0,0) = (3*s*(598360*s - 10840976*pow(s, 2) + 66905500*pow(s, 3) -
189738504*pow(s, 4) + 271189044*pow(s, 5) - 190169856*pow(s, 6) +
52056432*pow(s, 7) - 6830960*pow(t, 2) + 126762915*s*pow(t, 2) -
718006464*pow(s, 2)*pow(t, 2) + 1657574955*pow(s, 3)*pow(t, 2) -
1646873208*pow(s, 4)*pow(t, 2) + 585386676*pow(s, 5)*pow(t, 2) + 14778626*pow(t,
3) - 247491396*s*pow(t, 3) + 1232466660*pow(s, 2)*pow(t, 3) - 2161705320*pow(s,
3)*pow(t, 3) + 1184076792*pow(s, 4)*pow(t, 3) - 8390844*pow(t, 4) +
92584809*s*pow(t, 4) - 457214544*pow(s, 2)*pow(t, 4) + 546490800*pow(s,
3)*pow(t, 4) + 3672162*pow(t, 5) + 42791328*s*pow(t, 5) - 45020448*pow(s,
2)*pow(t, 5) - 6533136*pow(t, 6) - 15545196*s*pow(t, 6) + 3304152*pow(t,
7)))/4000.; 

  dphi(0,1) = (3*t*pow(s, 2)*(-6830960 + 84508610*s + 22167939*t -
247491396*s*t - 359003232*pow(s, 2) + 924349995*t*pow(s, 2) + 663029982*pow(s,
3) - 1297023192*t*pow(s, 3) - 548957736*pow(s, 4) + 592038396*t*pow(s, 4) +
167253336*pow(s, 5) - 16781688*pow(t, 2) + 123446412*s*pow(t, 2) -
457214544*pow(s, 2)*pow(t, 2) + 437192640*pow(s, 3)*pow(t, 2) + 9180405*pow(t,
3) + 71318880*s*pow(t, 3) - 56275560*pow(s, 2)*pow(t, 3) - 19599408*pow(t, 4) -
31090392*s*pow(t, 4) + 11564532*pow(t, 5)))/4000.; 

  dphi(1,0) = (3*s*pow(t,
2)*(-6830960 + 22167939*s + 84508610*t - 247491396*s*t - 16781688*pow(s, 2) +
123446412*t*pow(s, 2) + 9180405*pow(s, 3) + 71318880*t*pow(s, 3) -
19599408*pow(s, 4) - 31090392*t*pow(s, 4) + 11564532*pow(s, 5) -
359003232*pow(t, 2) + 924349995*s*pow(t, 2) - 457214544*pow(s, 2)*pow(t, 2) -
56275560*pow(s, 3)*pow(t, 2) + 663029982*pow(t, 3) - 1297023192*s*pow(t, 3) +
437192640*pow(s, 2)*pow(t, 3) - 548957736*pow(t, 4) + 592038396*s*pow(t, 4) +
167253336*pow(t, 5)))/4000.; 

  dphi(1,1) = (3*t*(598360*t - 6830960*pow(s, 2) +
126762915*t*pow(s, 2) + 14778626*pow(s, 3) - 247491396*t*pow(s, 3) -
8390844*pow(s, 4) + 92584809*t*pow(s, 4) + 3672162*pow(s, 5) + 42791328*t*pow(s,
5) - 6533136*pow(s, 6) - 15545196*t*pow(s, 6) + 3304152*pow(s, 7) -
10840976*pow(t, 2) - 718006464*pow(s, 2)*pow(t, 2) + 1232466660*pow(s, 3)*pow(t,
2) - 457214544*pow(s, 4)*pow(t, 2) - 45020448*pow(s, 5)*pow(t, 2) +
66905500*pow(t, 3) + 1657574955*pow(s, 2)*pow(t, 3) - 2161705320*pow(s,
3)*pow(t, 3) + 546490800*pow(s, 4)*pow(t, 3) - 189738504*pow(t, 4) -
1646873208*pow(s, 2)*pow(t, 4) + 1184076792*pow(s, 3)*pow(t, 4) +
271189044*pow(t, 5) + 585386676*pow(s, 2)*pow(t, 5) - 190169856*pow(t, 6) +
52056432*pow(t, 7)))/4000.; 

  dphi(2,0) = (-3*s*(-1 + s + t)*(-1986086*s -
1986086*s*t + 19209962*pow(s, 2) + 17223876*t*pow(s, 2) - 67816998*pow(s, 3) -
50593122*t*pow(s, 3) + 111677490*pow(s, 4) + 61084368*t*pow(s, 4) -
87112584*pow(s, 5) - 26028216*t*pow(s, 5) + 26028216*pow(s, 6) + 11244812*pow(t,
2) - 71777907*s*pow(t, 2) + 82465353*pow(s, 2)*pow(t, 2) + 77176476*pow(s,
3)*pow(t, 2) - 109230768*pow(s, 4)*pow(t, 2) - 42779610*pow(t, 3) +
321155847*s*pow(t, 3) - 505693692*pow(s, 2)*pow(t, 3) + 181681704*pow(s,
3)*pow(t, 3) + 25730082*pow(t, 4) - 335100240*s*pow(t, 4) + 352372356*pow(s,
2)*pow(t, 4) + 43851780*pow(t, 5) + 74870892*s*pow(t, 5) - 38645424*pow(t,
6)))/2000.; 

  dphi(2,1) = (3*t*(-1 + s + t)*(1986086*t + 1986086*s*t -
11244812*pow(s, 2) + 71777907*t*pow(s, 2) + 42779610*pow(s, 3) -
321155847*t*pow(s, 3) - 25730082*pow(s, 4) + 335100240*t*pow(s, 4) -
43851780*pow(s, 5) - 74870892*t*pow(s, 5) + 38645424*pow(s, 6) - 19209962*pow(t,
2) - 17223876*s*pow(t, 2) - 82465353*pow(s, 2)*pow(t, 2) + 505693692*pow(s,
3)*pow(t, 2) - 352372356*pow(s, 4)*pow(t, 2) + 67816998*pow(t, 3) +
50593122*s*pow(t, 3) - 77176476*pow(s, 2)*pow(t, 3) - 181681704*pow(s, 3)*pow(t,
3) - 111677490*pow(t, 4) - 61084368*s*pow(t, 4) + 109230768*pow(s, 2)*pow(t, 4)
+ 87112584*pow(t, 5) + 26028216*s*pow(t, 5) - 26028216*pow(t, 6)))/2000.;


  dphi(3,0) = -(s*(18840*s - 342544*pow(s, 2) + 2125500*pow(s, 3) - 6071976*pow(s,
4) + 8755236*pow(s, 5) - 6200064*pow(s, 6) + 1714608*pow(s, 7) - 207476*pow(t,
2) + 3854289*s*pow(t, 2) - 21879216*pow(s, 2)*pow(t, 2) + 50684895*pow(s,
3)*pow(t, 2) - 50608152*pow(s, 4)*pow(t, 2) + 18105444*pow(s, 5)*pow(t, 2) +
445130*pow(t, 3) - 7445574*s*pow(t, 3) + 37111140*pow(s, 2)*pow(t, 3) -
65179080*pow(s, 3)*pow(t, 3) + 35802648*pow(s, 4)*pow(t, 3) - 254376*pow(t, 4) +
2759481*s*pow(t, 4) - 13661136*pow(s, 2)*pow(t, 4) + 16297200*pow(s, 3)*pow(t,
4) + 118458*pow(t, 5) + 1288872*s*pow(t, 5) - 1324512*pow(s, 2)*pow(t, 5) -
204768*pow(t, 6) - 475308*s*pow(t, 6) + 103032*pow(t, 7)))/400.; 

  dphi(3,1) =
-(t*pow(s, 2)*(-207476 + 2569526*s + 667695*t - 7445574*s*t - 10939608*pow(s, 2)
+ 27833355*t*pow(s, 2) + 20273958*pow(s, 3) - 39107448*t*pow(s, 3) -
16869384*pow(s, 4) + 17901324*t*pow(s, 4) + 5172984*pow(s, 5) - 508752*pow(t, 2)
+ 3679308*s*pow(t, 2) - 13661136*pow(s, 2)*pow(t, 2) + 13037760*pow(s, 3)*pow(t,
2) + 296145*pow(t, 3) + 2148120*s*pow(t, 3) - 1655640*pow(s, 2)*pow(t, 3) -
614304*pow(t, 4) - 950616*s*pow(t, 4) + 360612*pow(t, 5)))/400.; 

  dphi(4,0) =
-(s*t*(13360 - 293388*s - 18976*t + 345093*s*t + 2147744*pow(s, 2) -
1675800*t*pow(s, 2) - 7229700*pow(s, 3) + 3438135*t*pow(s, 3) + 12258000*pow(s,
4) - 3608064*t*pow(s, 4) - 10151568*pow(s, 5) + 1589868*t*pow(s, 5) +
3255552*pow(s, 6) + 110646*pow(t, 2) - 1915110*s*pow(t, 2) + 6560028*pow(s,
2)*pow(t, 2) - 5972940*pow(s, 3)*pow(t, 2) + 1207224*pow(s, 4)*pow(t, 2) -
244620*pow(t, 3) + 4806621*s*pow(t, 3) - 14258592*pow(s, 2)*pow(t, 3) +
8592480*pow(s, 3)*pow(t, 3) + 37854*pow(t, 4) - 3399084*s*pow(t, 4) +
6951744*pow(s, 2)*pow(t, 4) + 204768*pow(t, 5) + 475308*s*pow(t, 5) -
103032*pow(t, 6)))/400.; 

  dphi(4,1) = -(pow(s, 2)*(6680 - 97796*s - 18976*t +
230062*s*t + 536936*pow(s, 2) - 837900*t*pow(s, 2) - 1445940*pow(s, 3) +
1375254*t*pow(s, 3) + 2043000*pow(s, 4) - 1202688*t*pow(s, 4) - 1450224*pow(s,
5) + 454248*t*pow(s, 5) + 406944*pow(s, 6) + 165969*pow(t, 2) - 1915110*s*pow(t,
2) + 4920021*pow(s, 2)*pow(t, 2) - 3583764*pow(s, 3)*pow(t, 2) + 603612*pow(s,
4)*pow(t, 2) - 489240*pow(t, 3) + 6408828*s*pow(t, 3) - 14258592*pow(s,
2)*pow(t, 3) + 6873984*pow(s, 3)*pow(t, 3) + 94635*pow(t, 4) - 5665140*s*pow(t,
4) + 8689680*pow(s, 2)*pow(t, 4) + 614304*pow(t, 5) + 950616*s*pow(t, 5) -
360612*pow(t, 6)))/400.; 

  dphi(5,0) = -(pow(t, 2)*(6680 - 18976*s - 97796*t +
230062*s*t + 165969*pow(s, 2) - 1915110*t*pow(s, 2) - 489240*pow(s, 3) +
6408828*t*pow(s, 3) + 94635*pow(s, 4) - 5665140*t*pow(s, 4) + 614304*pow(s, 5) +
950616*t*pow(s, 5) - 360612*pow(s, 6) + 536936*pow(t, 2) - 837900*s*pow(t, 2) +
4920021*pow(s, 2)*pow(t, 2) - 14258592*pow(s, 3)*pow(t, 2) + 8689680*pow(s,
4)*pow(t, 2) - 1445940*pow(t, 3) + 1375254*s*pow(t, 3) - 3583764*pow(s,
2)*pow(t, 3) + 6873984*pow(s, 3)*pow(t, 3) + 2043000*pow(t, 4) -
1202688*s*pow(t, 4) + 603612*pow(s, 2)*pow(t, 4) - 1450224*pow(t, 5) +
454248*s*pow(t, 5) + 406944*pow(t, 6)))/400.; 

  dphi(5,1) = (s*t*(-13360 + 18976*s
+ 293388*t - 345093*s*t - 110646*pow(s, 2) + 1915110*t*pow(s, 2) + 244620*pow(s,
3) - 4806621*t*pow(s, 3) - 37854*pow(s, 4) + 3399084*t*pow(s, 4) - 204768*pow(s,
5) - 475308*t*pow(s, 5) + 103032*pow(s, 6) - 2147744*pow(t, 2) +
1675800*s*pow(t, 2) - 6560028*pow(s, 2)*pow(t, 2) + 14258592*pow(s, 3)*pow(t, 2)
- 6951744*pow(s, 4)*pow(t, 2) + 7229700*pow(t, 3) - 3438135*s*pow(t, 3) +
5972940*pow(s, 2)*pow(t, 3) - 8592480*pow(s, 3)*pow(t, 3) - 12258000*pow(t, 4) +
3608064*s*pow(t, 4) - 1207224*pow(s, 2)*pow(t, 4) + 10151568*pow(t, 5) -
1589868*s*pow(t, 5) - 3255552*pow(t, 6)))/400.; 

  dphi(6,0) = -(s*pow(t,
2)*(-207476 + 667695*s + 2569526*t - 7445574*s*t - 508752*pow(s, 2) +
3679308*t*pow(s, 2) + 296145*pow(s, 3) + 2148120*t*pow(s, 3) - 614304*pow(s, 4)
- 950616*t*pow(s, 4) + 360612*pow(s, 5) - 10939608*pow(t, 2) + 27833355*s*pow(t,
2) - 13661136*pow(s, 2)*pow(t, 2) - 1655640*pow(s, 3)*pow(t, 2) +
20273958*pow(t, 3) - 39107448*s*pow(t, 3) + 13037760*pow(s, 2)*pow(t, 3) -
16869384*pow(t, 4) + 17901324*s*pow(t, 4) + 5172984*pow(t, 5)))/400.; 

  dphi(6,1)
= -(t*(18840*t - 207476*pow(s, 2) + 3854289*t*pow(s, 2) + 445130*pow(s, 3) -
7445574*t*pow(s, 3) - 254376*pow(s, 4) + 2759481*t*pow(s, 4) + 118458*pow(s, 5)
+ 1288872*t*pow(s, 5) - 204768*pow(s, 6) - 475308*t*pow(s, 6) + 103032*pow(s, 7)
- 342544*pow(t, 2) - 21879216*pow(s, 2)*pow(t, 2) + 37111140*pow(s, 3)*pow(t, 2)
- 13661136*pow(s, 4)*pow(t, 2) - 1324512*pow(s, 5)*pow(t, 2) + 2125500*pow(t, 3)
+ 50684895*pow(s, 2)*pow(t, 3) - 65179080*pow(s, 3)*pow(t, 3) + 16297200*pow(s,
4)*pow(t, 3) - 6071976*pow(t, 4) - 50608152*pow(s, 2)*pow(t, 4) +
35802648*pow(s, 3)*pow(t, 4) + 8755236*pow(t, 5) + 18105444*pow(s, 2)*pow(t, 5)
- 6200064*pow(t, 6) + 1714608*pow(t, 7)))/400.; 

  dphi(7,0) = -((-1 + s + t)*(200
  + 200*s + 200*t + 400*s*t - 80134*pow(s, 2) - 79734*t*pow(s, 2) +
699978*pow(s, 3) + 620244*t*pow(s, 3) - 2357262*pow(s, 4) - 1737018*t*pow(s, 4)
+ 3780810*pow(s, 5) + 2043792*t*pow(s, 5) - 2901096*pow(s, 6) - 857304*t*pow(s,
6) + 857304*pow(s, 7) - 29918*pow(t, 2) + 523968*s*pow(t, 2) - 1789389*pow(s,
2)*pow(t, 2) - 267453*pow(s, 3)*pow(t, 2) + 6328044*pow(s, 4)*pow(t, 2) -
5089392*pow(s, 5)*pow(t, 2) + 210868*pow(t, 3) - 3191508*s*pow(t, 3) +
11235411*pow(s, 2)*pow(t, 3) - 11357658*pow(s, 3)*pow(t, 3) + 1514376*pow(s,
4)*pow(t, 3) - 629370*pow(t, 4) + 7117704*s*pow(t, 4) - 17931294*pow(s,
2)*pow(t, 4) + 11312784*pow(s, 3)*pow(t, 4) + 943740*pow(t, 5) -
6930360*s*pow(t, 5) + 8812800*pow(s, 2)*pow(t, 5) - 699192*pow(t, 6) +
2489616*s*pow(t, 6) + 203472*pow(t, 7)))/200.; 

  dphi(7,1) = (s*t*(-1 + s +
t)*(30118 - 246625*s - 331061*t + 2367072*s*t + 497916*pow(s, 2) -
5243166*t*pow(s, 2) + 272493*pow(s, 3) + 3401433*t*pow(s, 3) - 1394010*pow(s, 4)
+ 44226*t*pow(s, 4) + 849528*pow(s, 5) + 1349415*pow(t, 2) - 7222095*s*pow(t, 2)
+ 11724345*pow(s, 2)*pow(t, 2) - 5175090*pow(s, 3)*pow(t, 2) - 2583360*pow(t, 3)
+ 8934300*s*pow(t, 3) - 7403400*pow(s, 2)*pow(t, 3) + 2345436*pow(t, 4) -
3899016*s*pow(t, 4) - 813888*pow(t, 5)))/100.; 

  dphi(8,0) = -(s*t*(-1 + s +
t)*(-30118 + 331061*s + 246625*t - 2367072*s*t - 1349415*pow(s, 2) +
7222095*t*pow(s, 2) + 2583360*pow(s, 3) - 8934300*t*pow(s, 3) - 2345436*pow(s,
4) + 3899016*t*pow(s, 4) + 813888*pow(s, 5) - 497916*pow(t, 2) +
5243166*s*pow(t, 2) - 11724345*pow(s, 2)*pow(t, 2) + 7403400*pow(s, 3)*pow(t, 2)
- 272493*pow(t, 3) - 3401433*s*pow(t, 3) + 5175090*pow(s, 2)*pow(t, 3) +
1394010*pow(t, 4) - 44226*s*pow(t, 4) - 849528*pow(t, 5)))/100.; 

  dphi(8,1) =
-((-1 + s + t)*(200 + 200*s + 200*t + 400*s*t - 29918*pow(s, 2) +
523968*t*pow(s, 2) + 210868*pow(s, 3) - 3191508*t*pow(s, 3) - 629370*pow(s, 4) +
7117704*t*pow(s, 4) + 943740*pow(s, 5) - 6930360*t*pow(s, 5) - 699192*pow(s, 6)
+ 2489616*t*pow(s, 6) + 203472*pow(s, 7) - 80134*pow(t, 2) - 79734*s*pow(t, 2) -
1789389*pow(s, 2)*pow(t, 2) + 11235411*pow(s, 3)*pow(t, 2) - 17931294*pow(s,
4)*pow(t, 2) + 8812800*pow(s, 5)*pow(t, 2) + 699978*pow(t, 3) + 620244*s*pow(t,
3) - 267453*pow(s, 2)*pow(t, 3) - 11357658*pow(s, 3)*pow(t, 3) + 11312784*pow(s,
4)*pow(t, 3) - 2357262*pow(t, 4) - 1737018*s*pow(t, 4) + 6328044*pow(s,
2)*pow(t, 4) + 1514376*pow(s, 3)*pow(t, 4) + 3780810*pow(t, 5) +
2043792*s*pow(t, 5) - 5089392*pow(s, 2)*pow(t, 5) - 2901096*pow(t, 6) -
857304*s*pow(t, 6) + 857304*pow(t, 7)))/200.; 

  dphi(9,0) = (s*(120*s -
2192*pow(s, 2) + 13700*pow(s, 3) - 39528*pow(s, 4) + 57708*pow(s, 5) -
41472*pow(s, 6) + 11664*pow(s, 7) - 1256*pow(t, 2) + 23349*s*pow(t, 2) -
132768*pow(s, 2)*pow(t, 2) + 308385*pow(s, 3)*pow(t, 2) - 309096*pow(s,
4)*pow(t, 2) + 111132*pow(s, 5)*pow(t, 2) + 2678*pow(t, 3) - 44712*s*pow(t, 3) +
223020*pow(s, 2)*pow(t, 3) - 392040*pow(s, 3)*pow(t, 3) + 215784*pow(s,
4)*pow(t, 3) - 1548*pow(t, 4) + 16443*s*pow(t, 4) - 81648*pow(s, 2)*pow(t, 4) +
97200*pow(s, 3)*pow(t, 4) + 774*pow(t, 5) + 7776*s*pow(t, 5) - 7776*pow(s,
2)*pow(t, 5) - 1296*pow(t, 6) - 2916*s*pow(t, 6) + 648*pow(t, 7)))/80.;


  dphi(9,1) = (t*pow(s, 2)*(-1256 + 15566*s + 4017*t - 44712*s*t - 66384*pow(s, 2)
+ 167265*t*pow(s, 2) + 123354*pow(s, 3) - 235224*t*pow(s, 3) - 103032*pow(s, 4)
+ 107892*t*pow(s, 4) + 31752*pow(s, 5) - 3096*pow(t, 2) + 21924*s*pow(t, 2) -
81648*pow(s, 2)*pow(t, 2) + 77760*pow(s, 3)*pow(t, 2) + 1935*pow(t, 3) +
12960*s*pow(t, 3) - 9720*pow(s, 2)*pow(t, 3) - 3888*pow(t, 4) - 5832*s*pow(t, 4)
+ 2268*pow(t, 5)))/80.; 

  dphi(10,0) = (s*t*(80 - 1764*s - 100*t + 1797*s*t +
12992*pow(s, 2) - 8280*t*pow(s, 2) - 44100*pow(s, 3) + 15705*t*pow(s, 3) +
75600*pow(s, 4) - 15552*t*pow(s, 4) - 63504*pow(s, 5) + 6804*t*pow(s, 5) +
20736*pow(s, 6) + 650*pow(t, 2) - 11340*s*pow(t, 2) + 38484*pow(s, 2)*pow(t, 2)
- 34020*pow(s, 3)*pow(t, 2) + 5832*pow(s, 4)*pow(t, 2) - 1440*pow(t, 3) +
28863*s*pow(t, 3) - 85536*pow(s, 2)*pow(t, 3) + 51840*pow(s, 3)*pow(t, 3) +
162*pow(t, 4) - 20412*s*pow(t, 4) + 41472*pow(s, 2)*pow(t, 4) + 1296*pow(t, 5) +
2916*s*pow(t, 5) - 648*pow(t, 6)))/40.; 

  dphi(10,1) = (pow(s, 2)*(40 - 588*s -
100*t + 1198*s*t + 3248*pow(s, 2) - 4140*t*pow(s, 2) - 8820*pow(s, 3) +
6282*t*pow(s, 3) + 12600*pow(s, 4) - 5184*t*pow(s, 4) - 9072*pow(s, 5) +
1944*t*pow(s, 5) + 2592*pow(s, 6) + 975*pow(t, 2) - 11340*s*pow(t, 2) +
28863*pow(s, 2)*pow(t, 2) - 20412*pow(s, 3)*pow(t, 2) + 2916*pow(s, 4)*pow(t, 2)
- 2880*pow(t, 3) + 38484*s*pow(t, 3) - 85536*pow(s, 2)*pow(t, 3) + 41472*pow(s,
3)*pow(t, 3) + 405*pow(t, 4) - 34020*s*pow(t, 4) + 51840*pow(s, 2)*pow(t, 4) +
3888*pow(t, 5) + 5832*s*pow(t, 5) - 2268*pow(t, 6)))/40.; 

  dphi(11,0) = (s*pow(t,
2)*(20 - 153*s - 650*t + 11340*s*t - 2520*pow(s, 2) - 38484*t*pow(s, 2) +
14895*pow(s, 3) + 34020*t*pow(s, 3) - 23328*pow(s, 4) - 5832*t*pow(s, 4) +
11340*pow(s, 5) + 1440*pow(t, 2) - 28863*s*pow(t, 2) + 85536*pow(s, 2)*pow(t, 2)
- 51840*pow(s, 3)*pow(t, 2) - 162*pow(t, 3) + 20412*s*pow(t, 3) - 41472*pow(s,
2)*pow(t, 3) - 1296*pow(t, 4) - 2916*s*pow(t, 4) + 648*pow(t, 5)))/80.;


  dphi(11,1) = (t*pow(s, 2)*(20 - 102*s - 975*t + 11340*s*t - 1260*pow(s, 2) -
28863*t*pow(s, 2) + 5958*pow(s, 3) + 20412*t*pow(s, 3) - 7776*pow(s, 4) -
2916*t*pow(s, 4) + 3240*pow(s, 5) + 2880*pow(t, 2) - 38484*s*pow(t, 2) +
85536*pow(s, 2)*pow(t, 2) - 41472*pow(s, 3)*pow(t, 2) - 405*pow(t, 3) +
34020*s*pow(t, 3) - 51840*pow(s, 2)*pow(t, 3) - 3888*pow(t, 4) - 5832*s*pow(t,
4) + 2268*pow(t, 5)))/80.; 

  dphi(12,0) = (s*pow(t, 2)*(20 - 975*s - 102*t +
11340*s*t + 2880*pow(s, 2) - 38484*t*pow(s, 2) - 405*pow(s, 3) + 34020*t*pow(s,
3) - 3888*pow(s, 4) - 5832*t*pow(s, 4) + 2268*pow(s, 5) - 1260*pow(t, 2) -
28863*s*pow(t, 2) + 85536*pow(s, 2)*pow(t, 2) - 51840*pow(s, 3)*pow(t, 2) +
5958*pow(t, 3) + 20412*s*pow(t, 3) - 41472*pow(s, 2)*pow(t, 3) - 7776*pow(t, 4)
- 2916*s*pow(t, 4) + 3240*pow(t, 5)))/80.; 

  dphi(12,1) = (t*pow(s, 2)*(20 - 650*s
  - 153*t + 11340*s*t + 1440*pow(s, 2) - 28863*t*pow(s, 2) - 162*pow(s, 3) +
    20412*t*pow(s, 3) - 1296*pow(s, 4) - 2916*t*pow(s, 4) + 648*pow(s, 5) -
2520*pow(t, 2) - 38484*s*pow(t, 2) + 85536*pow(s, 2)*pow(t, 2) - 41472*pow(s,
3)*pow(t, 2) + 14895*pow(t, 3) + 34020*s*pow(t, 3) - 51840*pow(s, 2)*pow(t, 3) -
23328*pow(t, 4) - 5832*s*pow(t, 4) + 11340*pow(t, 5)))/80.; 

  dphi(13,0) = (pow(t,
2)*(40 - 100*s - 588*t + 1198*s*t + 975*pow(s, 2) - 11340*t*pow(s, 2) -
2880*pow(s, 3) + 38484*t*pow(s, 3) + 405*pow(s, 4) - 34020*t*pow(s, 4) +
3888*pow(s, 5) + 5832*t*pow(s, 5) - 2268*pow(s, 6) + 3248*pow(t, 2) -
4140*s*pow(t, 2) + 28863*pow(s, 2)*pow(t, 2) - 85536*pow(s, 3)*pow(t, 2) +
51840*pow(s, 4)*pow(t, 2) - 8820*pow(t, 3) + 6282*s*pow(t, 3) - 20412*pow(s,
2)*pow(t, 3) + 41472*pow(s, 3)*pow(t, 3) + 12600*pow(t, 4) - 5184*s*pow(t, 4) +
2916*pow(s, 2)*pow(t, 4) - 9072*pow(t, 5) + 1944*s*pow(t, 5) + 2592*pow(t,
6)))/40.; 

  dphi(13,1) = -(s*t*(-80 + 100*s + 1764*t - 1797*s*t - 650*pow(s, 2) +
11340*t*pow(s, 2) + 1440*pow(s, 3) - 28863*t*pow(s, 3) - 162*pow(s, 4) +
20412*t*pow(s, 4) - 1296*pow(s, 5) - 2916*t*pow(s, 5) + 648*pow(s, 6) -
12992*pow(t, 2) + 8280*s*pow(t, 2) - 38484*pow(s, 2)*pow(t, 2) + 85536*pow(s,
3)*pow(t, 2) - 41472*pow(s, 4)*pow(t, 2) + 44100*pow(t, 3) - 15705*s*pow(t, 3) +
34020*pow(s, 2)*pow(t, 3) - 51840*pow(s, 3)*pow(t, 3) - 75600*pow(t, 4) +
15552*s*pow(t, 4) - 5832*pow(s, 2)*pow(t, 4) + 63504*pow(t, 5) - 6804*s*pow(t,
5) - 20736*pow(t, 6)))/40.; 

  dphi(14,0) = (s*pow(t, 2)*(-1256 + 4017*s + 15566*t
- 44712*s*t - 3096*pow(s, 2) + 21924*t*pow(s, 2) + 1935*pow(s, 3) +
12960*t*pow(s, 3) - 3888*pow(s, 4) - 5832*t*pow(s, 4) + 2268*pow(s, 5) -
66384*pow(t, 2) + 167265*s*pow(t, 2) - 81648*pow(s, 2)*pow(t, 2) - 9720*pow(s,
3)*pow(t, 2) + 123354*pow(t, 3) - 235224*s*pow(t, 3) + 77760*pow(s, 2)*pow(t, 3)
- 103032*pow(t, 4) + 107892*s*pow(t, 4) + 31752*pow(t, 5)))/80.; 

  dphi(14,1) =
  (t*(120*t - 1256*pow(s, 2) + 23349*t*pow(s, 2) + 2678*pow(s, 3) -
44712*t*pow(s, 3) - 1548*pow(s, 4) + 16443*t*pow(s, 4) + 774*pow(s, 5) +
7776*t*pow(s, 5) - 1296*pow(s, 6) - 2916*t*pow(s, 6) + 648*pow(s, 7) -
2192*pow(t, 2) - 132768*pow(s, 2)*pow(t, 2) + 223020*pow(s, 3)*pow(t, 2) -
81648*pow(s, 4)*pow(t, 2) - 7776*pow(s, 5)*pow(t, 2) + 13700*pow(t, 3) +
308385*pow(s, 2)*pow(t, 3) - 392040*pow(s, 3)*pow(t, 3) + 97200*pow(s, 4)*pow(t,
3) - 39528*pow(t, 4) - 309096*pow(s, 2)*pow(t, 4) + 215784*pow(s, 3)*pow(t, 4) +
57708*pow(t, 5) + 111132*pow(s, 2)*pow(t, 5) - 41472*pow(t, 6) + 11664*pow(t,
7)))/80.; 

  dphi(15,0) = -(s*(-1 + s + t)*(40 - 842*s + 40*t - 802*s*t +
5734*pow(s, 2) + 4932*t*pow(s, 2) - 17586*pow(s, 3) - 12654*t*pow(s, 3) +
26910*pow(s, 4) + 14256*t*pow(s, 4) - 20088*pow(s, 5) - 5832*t*pow(s, 5) +
5832*pow(s, 6) - 638*pow(t, 2) + 11133*s*pow(t, 2) - 54459*pow(s, 2)*pow(t, 2) +
94932*pow(s, 3)*pow(t, 2) - 53136*pow(s, 4)*pow(t, 2) + 1026*pow(t, 3) -
15867*s*pow(t, 3) + 60426*pow(s, 2)*pow(t, 3) - 57672*pow(s, 3)*pow(t, 3) -
468*pow(t, 4) + 6318*s*pow(t, 4) - 16848*pow(s, 2)*pow(t, 4)))/40.; 

  dphi(15,1) =
(3*(-1 + 3*s)*(-1 + 6*s)*t*(-1 + s + t)*pow(s, 2)*(113 - 267*s - 303*t + 357*s*t
+ 156*pow(s, 2) + 195*pow(t, 2)))/20.; 

  dphi(16,0) = -(t*(-1 + s + t)*(10 - 304*s
- 147*t + 3385*s*t + 2573*pow(s, 2) - 21195*t*pow(s, 2) - 9495*pow(s, 3) +
54990*t*pow(s, 3) + 17280*pow(s, 4) - 63180*t*pow(s, 4) - 15228*pow(s, 5) +
26568*t*pow(s, 5) + 5184*pow(s, 6) + 812*pow(t, 2) - 13905*s*pow(t, 2) +
61290*pow(s, 2)*pow(t, 2) - 100440*pow(s, 3)*pow(t, 2) + 55080*pow(s, 4)*pow(t,
2) - 2205*pow(t, 3) + 26730*s*pow(t, 3) - 74520*pow(s, 2)*pow(t, 3) +
58320*pow(s, 3)*pow(t, 3) + 3150*pow(t, 4) - 24300*s*pow(t, 4) + 32400*pow(s,
2)*pow(t, 4) - 2268*pow(t, 5) + 8424*s*pow(t, 5) + 648*pow(t, 6)))/10.;


  dphi(16,1) = -(s*(-1 + s + t)*(10 - 147*s - 304*t + 3385*s*t + 812*pow(s, 2) -
13905*t*pow(s, 2) - 2205*pow(s, 3) + 26730*t*pow(s, 3) + 3150*pow(s, 4) -
24300*t*pow(s, 4) - 2268*pow(s, 5) + 8424*t*pow(s, 5) + 648*pow(s, 6) +
2573*pow(t, 2) - 21195*s*pow(t, 2) + 61290*pow(s, 2)*pow(t, 2) - 74520*pow(s,
3)*pow(t, 2) + 32400*pow(s, 4)*pow(t, 2) - 9495*pow(t, 3) + 54990*s*pow(t, 3) -
100440*pow(s, 2)*pow(t, 3) + 58320*pow(s, 3)*pow(t, 3) + 17280*pow(t, 4) -
63180*s*pow(t, 4) + 55080*pow(s, 2)*pow(t, 4) - 15228*pow(t, 5) + 26568*s*pow(t,
5) + 5184*pow(t, 6)))/10.; 

  dphi(17,0) = (3*s*(-1 + s + t)*(-1 + 3*t)*(-1 +
6*t)*pow(t, 2)*(113 - 303*s - 267*t + 357*s*t + 195*pow(s, 2) + 156*pow(t,
2)))/20.; 

  dphi(17,1) = -(t*(-1 + s + t)*(40 + 40*s - 842*t - 802*s*t -
638*pow(s, 2) + 11133*t*pow(s, 2) + 1026*pow(s, 3) - 15867*t*pow(s, 3) -
468*pow(s, 4) + 6318*t*pow(s, 4) + 5734*pow(t, 2) + 4932*s*pow(t, 2) -
54459*pow(s, 2)*pow(t, 2) + 60426*pow(s, 3)*pow(t, 2) - 16848*pow(s, 4)*pow(t,
2) - 17586*pow(t, 3) - 12654*s*pow(t, 3) + 94932*pow(s, 2)*pow(t, 3) -
57672*pow(s, 3)*pow(t, 3) + 26910*pow(t, 4) + 14256*s*pow(t, 4) - 53136*pow(s,
2)*pow(t, 4) - 20088*pow(t, 5) - 5832*s*pow(t, 5) + 5832*pow(t, 6)))/40.;


  dphi(18,0) = -16*(-1 + s + t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(-10 + 84*s + 37*t
- 207*s*t - 162*pow(s, 2) + 198*t*pow(s, 2) + 90*pow(s, 3) - 45*pow(t, 2) +
126*s*pow(t, 2) + 18*pow(t, 3)); 

  dphi(18,1) = -16*s*t*(-1 + s + t)*(-1 + s +
2*t)*(20 - 54*s - 351*t + 837*s*t + 36*pow(s, 2) - 486*t*pow(s, 2) + 1647*pow(t,
2) - 3078*s*pow(t, 2) + 1296*pow(s, 2)*pow(t, 2) - 2592*pow(t, 3) +
2592*s*pow(t, 3) + 1296*pow(t, 4)); 

  dphi(19,0) = -16*s*t*(-1 + s + t)*(-1 + 2*s
+ t)*(20 - 351*s - 54*t + 837*s*t + 1647*pow(s, 2) - 3078*t*pow(s, 2) -
2592*pow(s, 3) + 2592*t*pow(s, 3) + 1296*pow(s, 4) + 36*pow(t, 2) - 486*s*pow(t,
2) + 1296*pow(s, 2)*pow(t, 2)); 

  dphi(19,1) = -16*(-1 + 3*s)*(-1 + 6*s)*(-1 + s +
t)*pow(s, 2)*(-10 + 37*s + 84*t - 207*s*t - 45*pow(s, 2) + 126*t*pow(s, 2) +
18*pow(s, 3) - 162*pow(t, 2) + 198*s*pow(t, 2) + 90*pow(t, 3)); 

  dphi(20,0) =
8*s*(-1 + 3*t)*(-1 + 6*t)*(-2 + 30*s + 2*t - 27*s*t - 108*pow(s, 2) +
72*t*pow(s, 2) + 90*pow(s, 3))*pow(t, 2)*sqrt(2); 

  dphi(20,1) = 8*(-1 + 3*s)*(-1
+ 6*s)*t*pow(s, 2)*(-2 + 2*s + 30*t - 27*s*t - 108*pow(t, 2) + 72*s*pow(t, 2) +
90*pow(t, 3))*sqrt(2); 

  dphi(21,0) = (-11664*s*(-1 + s + t)*(-1 + 3*t)*(-1 +
6*t)*pow(t, 2)*(93 - 303*s - 197*t + 357*s*t + 195*pow(s, 2) + 96*pow(t,
2)))/125.; 

  dphi(21,1) = (-5832*t*(-1 + s + t)*(-16*t - 16*s*t + 186*pow(s, 2) -
3211*t*pow(s, 2) - 342*pow(s, 3) + 5289*t*pow(s, 3) + 156*pow(s, 4) -
2106*t*pow(s, 4) + 272*pow(t, 2) + 256*s*pow(t, 2) + 15353*pow(s, 2)*pow(t, 2) -
20142*pow(s, 3)*pow(t, 2) + 5616*pow(s, 4)*pow(t, 2) - 1488*pow(t, 3) -
1232*s*pow(t, 3) - 25644*pow(s, 2)*pow(t, 3) + 19224*pow(s, 3)*pow(t, 3) +
3440*pow(t, 4) + 2208*s*pow(t, 4) + 13392*pow(s, 2)*pow(t, 4) - 3504*pow(t, 5) -
1296*s*pow(t, 5) + 1296*pow(t, 6)))/125.; 

  dphi(22,0) = (729*s*(-1 + s + t)*(-1 +
3*t)*(-1 + 6*t)*pow(t, 2)*(61 - 239*s - 107*t + 229*s*t + 195*pow(s, 2) +
44*pow(t, 2)))/8.; 

  dphi(22,1) = (729*t*(-1 + s + t)*(-10*t - 10*s*t + 122*pow(s,
2) - 2039*t*pow(s, 2) - 278*pow(s, 3) + 4233*t*pow(s, 3) + 156*pow(s, 4) -
2106*t*pow(s, 4) + 166*pow(t, 2) + 156*s*pow(t, 2) + 9161*pow(s, 2)*pow(t, 2) -
15534*pow(s, 3)*pow(t, 2) + 5616*pow(s, 4)*pow(t, 2) - 874*pow(t, 3) -
718*s*pow(t, 3) - 13692*pow(s, 2)*pow(t, 3) + 13464*pow(s, 3)*pow(t, 3) +
1918*pow(t, 4) + 1200*s*pow(t, 4) + 6192*pow(s, 2)*pow(t, 4) - 1848*pow(t, 5) -
648*s*pow(t, 5) + 648*pow(t, 6)))/16.; 

  dphi(23,0) = (-81*s*(-1 + s + t)*(-1 +
6*t)*pow(t, 2)*(-1003 + 6353*s + 2632*t - 12582*s*t - 11277*pow(s, 2) +
11727*t*pow(s, 2) + 6048*pow(s, 3) - 1971*pow(t, 2) + 5913*s*pow(t, 2) +
324*pow(t, 3)))/8.; 

  dphi(23,1) = (-81*t*(-1 + s + t)*(-180*t - 180*s*t +
2006*pow(s, 2) - 27133*t*pow(s, 2) - 7802*pow(s, 3) + 96471*t*pow(s, 3) +
9828*pow(s, 4) - 105750*t*pow(s, 4) - 4032*pow(s, 5) + 36288*t*pow(s, 5) +
2628*pow(t, 2) + 2448*s*pow(t, 2) + 80967*pow(s, 2)*pow(t, 2) - 222210*pow(s,
3)*pow(t, 2) + 134352*pow(s, 4)*pow(t, 2) - 11502*pow(t, 3) - 9054*s*pow(t, 3) -
77652*pow(s, 2)*pow(t, 3) + 130248*pow(s, 3)*pow(t, 3) + 21582*pow(t, 4) +
12528*s*pow(t, 4) + 19440*pow(s, 2)*pow(t, 4) - 18360*pow(t, 5) - 5832*s*pow(t,
5) + 5832*pow(t, 6)))/16.; 

  dphi(24,0) = (324*s*(-1 + s + t)*pow(t, 2)*(11858 -
100135*s - 53894*t + 399855*s*t + 268917*pow(s, 2) - 762444*t*pow(s, 2) -
292788*pow(s, 3) + 427788*t*pow(s, 3) + 112644*pow(s, 4) + 62844*pow(t, 2) -
448092*s*pow(t, 2) + 474552*pow(s, 2)*pow(t, 2) - 648*pow(t, 3) +
135432*s*pow(t, 3) - 20736*pow(t, 4)))/125.; 

  dphi(24,1) = (648*t*(-1 + s +
t)*(-720*t - 720*s*t + 5929*pow(s, 2) - 44105*t*pow(s, 2) - 31402*pow(s, 3) +
201435*t*pow(s, 3) + 60861*pow(s, 4) - 274410*t*pow(s, 4) - 51480*pow(s, 5) +
119016*t*pow(s, 5) + 16092*pow(s, 6) + 7632*pow(t, 2) + 6912*s*pow(t, 2) +
79545*pow(s, 2)*pow(t, 2) - 326214*pow(s, 3)*pow(t, 2) + 241920*pow(s, 4)*pow(t,
2) - 28368*pow(t, 3) - 21456*s*pow(t, 3) - 21276*pow(s, 2)*pow(t, 3) +
139320*pow(s, 3)*pow(t, 3) + 48240*pow(t, 4) + 26784*s*pow(t, 4) - 24624*pow(s,
2)*pow(t, 4) - 38448*pow(t, 5) - 11664*s*pow(t, 5) + 11664*pow(t, 6)))/125.;


  dphi(25,0) = (648*s*(-1 + s + t)*(-720*s - 720*s*t + 7632*pow(s, 2) +
6912*t*pow(s, 2) - 28368*pow(s, 3) - 21456*t*pow(s, 3) + 48240*pow(s, 4) +
26784*t*pow(s, 4) - 38448*pow(s, 5) - 11664*t*pow(s, 5) + 11664*pow(s, 6) +
5929*pow(t, 2) - 44105*s*pow(t, 2) + 79545*pow(s, 2)*pow(t, 2) - 21276*pow(s,
3)*pow(t, 2) - 24624*pow(s, 4)*pow(t, 2) - 31402*pow(t, 3) + 201435*s*pow(t, 3)
- 326214*pow(s, 2)*pow(t, 3) + 139320*pow(s, 3)*pow(t, 3) + 60861*pow(t, 4) -
274410*s*pow(t, 4) + 241920*pow(s, 2)*pow(t, 4) - 51480*pow(t, 5) +
119016*s*pow(t, 5) + 16092*pow(t, 6)))/125.; 

  dphi(25,1) = (-324*t*(-1 + s +
t)*pow(s, 2)*(-11858 + 53894*s + 100135*t - 399855*s*t - 62844*pow(s, 2) +
448092*t*pow(s, 2) + 648*pow(s, 3) - 135432*t*pow(s, 3) + 20736*pow(s, 4) -
268917*pow(t, 2) + 762444*s*pow(t, 2) - 474552*pow(s, 2)*pow(t, 2) +
292788*pow(t, 3) - 427788*s*pow(t, 3) - 112644*pow(t, 4)))/125.; 

  dphi(26,0) =
(-81*s*(-1 + s + t)*(-180*s - 180*s*t + 2628*pow(s, 2) + 2448*t*pow(s, 2) -
11502*pow(s, 3) - 9054*t*pow(s, 3) + 21582*pow(s, 4) + 12528*t*pow(s, 4) -
18360*pow(s, 5) - 5832*t*pow(s, 5) + 5832*pow(s, 6) + 2006*pow(t, 2) -
27133*s*pow(t, 2) + 80967*pow(s, 2)*pow(t, 2) - 77652*pow(s, 3)*pow(t, 2) +
19440*pow(s, 4)*pow(t, 2) - 7802*pow(t, 3) + 96471*s*pow(t, 3) - 222210*pow(s,
2)*pow(t, 3) + 130248*pow(s, 3)*pow(t, 3) + 9828*pow(t, 4) - 105750*s*pow(t, 4)
+ 134352*pow(s, 2)*pow(t, 4) - 4032*pow(t, 5) + 36288*s*pow(t, 5)))/16.;


  dphi(26,1) = (-81*(-1 + 6*s)*t*(-1 + s + t)*pow(s, 2)*(-1003 + 2632*s + 6353*t -
12582*s*t - 1971*pow(s, 2) + 5913*t*pow(s, 2) + 324*pow(s, 3) - 11277*pow(t, 2)
+ 11727*s*pow(t, 2) + 6048*pow(t, 3)))/8.; 

  dphi(27,0) = (729*s*(-1 + s +
t)*(-10*s - 10*s*t + 166*pow(s, 2) + 156*t*pow(s, 2) - 874*pow(s, 3) -
718*t*pow(s, 3) + 1918*pow(s, 4) + 1200*t*pow(s, 4) - 1848*pow(s, 5) -
648*t*pow(s, 5) + 648*pow(s, 6) + 122*pow(t, 2) - 2039*s*pow(t, 2) + 9161*pow(s,
2)*pow(t, 2) - 13692*pow(s, 3)*pow(t, 2) + 6192*pow(s, 4)*pow(t, 2) - 278*pow(t,
3) + 4233*s*pow(t, 3) - 15534*pow(s, 2)*pow(t, 3) + 13464*pow(s, 3)*pow(t, 3) +
156*pow(t, 4) - 2106*s*pow(t, 4) + 5616*pow(s, 2)*pow(t, 4)))/16.; 

  dphi(27,1) =
(729*(-1 + 3*s)*(-1 + 6*s)*t*(-1 + s + t)*pow(s, 2)*(61 - 107*s - 239*t +
229*s*t + 44*pow(s, 2) + 195*pow(t, 2)))/8.; 

  dphi(28,0) = (-5832*s*(-1 + s +
t)*(-16*s - 16*s*t + 272*pow(s, 2) + 256*t*pow(s, 2) - 1488*pow(s, 3) -
1232*t*pow(s, 3) + 3440*pow(s, 4) + 2208*t*pow(s, 4) - 3504*pow(s, 5) -
1296*t*pow(s, 5) + 1296*pow(s, 6) + 186*pow(t, 2) - 3211*s*pow(t, 2) +
15353*pow(s, 2)*pow(t, 2) - 25644*pow(s, 3)*pow(t, 2) + 13392*pow(s, 4)*pow(t,
2) - 342*pow(t, 3) + 5289*s*pow(t, 3) - 20142*pow(s, 2)*pow(t, 3) + 19224*pow(s,
3)*pow(t, 3) + 156*pow(t, 4) - 2106*s*pow(t, 4) + 5616*pow(s, 2)*pow(t,
4)))/125.; 

  dphi(28,1) = (-11664*(-1 + 3*s)*(-1 + 6*s)*t*(-1 + s + t)*pow(s,
2)*(93 - 197*s - 303*t + 357*s*t + 96*pow(s, 2) + 195*pow(t, 2)))/125.;


  dphi(29,0) = (-162*s*pow(t, 2)*(-1036 + 22449*s - 2626*t + 42000*s*t -
162800*pow(s, 2) - 87564*t*pow(s, 2) + 485295*pow(s, 3) - 100080*t*pow(s, 3) -
604800*pow(s, 4) + 178200*t*pow(s, 4) + 262332*pow(s, 5) + 9080*pow(t, 2) -
170253*s*pow(t, 2) + 504576*pow(s, 2)*pow(t, 2) - 294840*pow(s, 3)*pow(t, 2) -
2322*pow(t, 3) + 121392*s*pow(t, 3) - 253152*pow(s, 2)*pow(t, 3) - 6336*pow(t,
4) - 15876*s*pow(t, 4) + 3240*pow(t, 5)))/125.; 

  dphi(29,1) = (-162*t*pow(s,
2)*(-1036 + 14966*s - 3939*t + 42000*s*t - 81400*pow(s, 2) - 65673*t*pow(s, 2) +
194118*pow(s, 3) - 60048*t*pow(s, 3) - 201600*pow(s, 4) + 89100*t*pow(s, 4) +
74952*pow(s, 5) + 18160*pow(t, 2) - 227004*s*pow(t, 2) + 504576*pow(s, 2)*pow(t,
2) - 235872*pow(s, 3)*pow(t, 2) - 5805*pow(t, 3) + 202320*s*pow(t, 3) -
316440*pow(s, 2)*pow(t, 3) - 19008*pow(t, 4) - 31752*s*pow(t, 4) + 11340*pow(t,
5)))/125.; 

  dphi(30,0) = (-81*s*pow(t, 2)*(-548 + 9501*s + 4042*t - 68064*s*t -
44944*pow(s, 2) + 288348*t*pow(s, 2) + 79335*pow(s, 3) - 410220*t*pow(s, 3) -
54864*pow(s, 4) + 183384*t*pow(s, 4) + 11340*pow(s, 5) - 6764*pow(t, 2) +
113301*s*pow(t, 2) - 377568*pow(s, 2)*pow(t, 2) + 293760*pow(s, 3)*pow(t, 2) +
2334*pow(t, 3) - 60804*s*pow(t, 3) + 134784*pow(s, 2)*pow(t, 3) + 2016*pow(t, 4)
+ 6156*s*pow(t, 4) - 1080*pow(t, 5)))/32.; 

  dphi(30,1) = (-81*t*pow(s, 2)*(-548 +
6334*s + 6063*t - 68064*s*t - 22472*pow(s, 2) + 216261*t*pow(s, 2) +
31734*pow(s, 3) - 246132*t*pow(s, 3) - 18288*pow(s, 4) + 91692*t*pow(s, 4) +
3240*pow(s, 5) - 13528*pow(t, 2) + 151068*s*pow(t, 2) - 377568*pow(s, 2)*pow(t,
2) + 235008*pow(s, 3)*pow(t, 2) + 5835*pow(t, 3) - 101340*s*pow(t, 3) +
168480*pow(s, 2)*pow(t, 3) + 6048*pow(t, 4) + 12312*s*pow(t, 4) - 3780*pow(t,
5)))/32.; 

  dphi(31,0) = (81*s*pow(t, 2)*(548 - 6063*s - 6334*t + 68064*s*t +
13528*pow(s, 2) - 151068*t*pow(s, 2) - 5835*pow(s, 3) + 101340*t*pow(s, 3) -
6048*pow(s, 4) - 12312*t*pow(s, 4) + 3780*pow(s, 5) + 22472*pow(t, 2) -
216261*s*pow(t, 2) + 377568*pow(s, 2)*pow(t, 2) - 168480*pow(s, 3)*pow(t, 2) -
31734*pow(t, 3) + 246132*s*pow(t, 3) - 235008*pow(s, 2)*pow(t, 3) + 18288*pow(t,
4) - 91692*s*pow(t, 4) - 3240*pow(t, 5)))/32.; 

  dphi(31,1) = (81*t*pow(s, 2)*(548
- 4042*s - 9501*t + 68064*s*t + 6764*pow(s, 2) - 113301*t*pow(s, 2) -
2334*pow(s, 3) + 60804*t*pow(s, 3) - 2016*pow(s, 4) - 6156*t*pow(s, 4) +
1080*pow(s, 5) + 44944*pow(t, 2) - 288348*s*pow(t, 2) + 377568*pow(s, 2)*pow(t,
2) - 134784*pow(s, 3)*pow(t, 2) - 79335*pow(t, 3) + 410220*s*pow(t, 3) -
293760*pow(s, 2)*pow(t, 3) + 54864*pow(t, 4) - 183384*s*pow(t, 4) - 11340*pow(t,
5)))/32.; 

  dphi(32,0) = (-162*s*pow(t, 2)*(-1036 - 3939*s + 14966*t + 42000*s*t +
18160*pow(s, 2) - 227004*t*pow(s, 2) - 5805*pow(s, 3) + 202320*t*pow(s, 3) -
19008*pow(s, 4) - 31752*t*pow(s, 4) + 11340*pow(s, 5) - 81400*pow(t, 2) -
65673*s*pow(t, 2) + 504576*pow(s, 2)*pow(t, 2) - 316440*pow(s, 3)*pow(t, 2) +
194118*pow(t, 3) - 60048*s*pow(t, 3) - 235872*pow(s, 2)*pow(t, 3) -
201600*pow(t, 4) + 89100*s*pow(t, 4) + 74952*pow(t, 5)))/125.; 

  dphi(32,1) =
(-162*t*pow(s, 2)*(-1036 - 2626*s + 22449*t + 42000*s*t + 9080*pow(s, 2) -
170253*t*pow(s, 2) - 2322*pow(s, 3) + 121392*t*pow(s, 3) - 6336*pow(s, 4) -
15876*t*pow(s, 4) + 3240*pow(s, 5) - 162800*pow(t, 2) - 87564*s*pow(t, 2) +
504576*pow(s, 2)*pow(t, 2) - 253152*pow(s, 3)*pow(t, 2) + 485295*pow(t, 3) -
100080*s*pow(t, 3) - 294840*pow(s, 2)*pow(t, 3) - 604800*pow(t, 4) +
178200*s*pow(t, 4) + 262332*pow(t, 5)))/125.; 

  dphi(33,0) = (-648*(-1 + s +
t)*(-1 + 3*s + t)*(-1 + 2*t)*(-2 + 3*t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2))/25.;


  dphi(33,1) = (-648*s*t*(-1 + s + t)*(-4 + 4*s + 83*t - 75*s*t - 545*pow(t, 2) +
420*s*pow(t, 2) + 1530*pow(t, 3) - 900*s*pow(t, 3) - 1908*pow(t, 4) +
648*s*pow(t, 4) + 864*pow(t, 5)))/25.; 

  dphi(34,0) = (81*(-1 + s + t)*(-1 +
2*t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(5 - 27*s - 11*t + 30*s*t + 24*pow(s, 2) +
6*pow(t, 2)))/4.; 

  dphi(34,1) = (81*s*t*(-1 + s + t)*(-10 + 22*s + 203*t -
405*s*t - 12*pow(s, 2) + 198*t*pow(s, 2) - 1289*pow(t, 2) + 2178*s*pow(t, 2) -
864*pow(s, 2)*pow(t, 2) + 3456*pow(t, 3) - 4356*s*pow(t, 3) + 1080*pow(s,
2)*pow(t, 3) - 4068*pow(t, 4) + 2808*s*pow(t, 4) + 1728*pow(t, 5)))/4.;


  dphi(35,0) = (81*(-1 + s + t)*(-1 + 6*t)*pow(t, 2)*(10 - 124*s - 57*t + 523*s*t
+ 404*pow(s, 2) - 1116*t*pow(s, 2) - 504*pow(s, 3) + 684*t*pow(s, 3) +
216*pow(s, 4) + 119*pow(t, 2) - 720*s*pow(t, 2) + 756*pow(s, 2)*pow(t, 2) -
108*pow(t, 3) + 324*s*pow(t, 3) + 36*pow(t, 4)))/4.; 

  dphi(35,1) = (81*s*t*(-1 +
s + t)*(-20 + 114*s + 361*t - 1787*s*t - 238*pow(s, 2) + 3186*t*pow(s, 2) +
216*pow(s, 3) - 2412*t*pow(s, 3) - 72*pow(s, 4) + 648*t*pow(s, 4) - 1951*pow(t,
2) + 7434*s*pow(t, 2) - 9180*pow(s, 2)*pow(t, 2) + 3672*pow(s, 3)*pow(t, 2) +
4464*pow(t, 3) - 11412*s*pow(t, 3) + 7128*pow(s, 2)*pow(t, 3) - 4572*pow(t, 4) +
5832*s*pow(t, 4) + 1728*pow(t, 5)))/4.; 

  dphi(36,0) = (-648*(-1 + s + t)*pow(t,
2)*(-10 + 184*s + 87*t - 1237*s*t - 947*pow(s, 2) + 4611*t*pow(s, 2) +
2073*pow(s, 3) - 6516*t*pow(s, 3) - 2052*pow(s, 4) + 3132*t*pow(s, 4) +
756*pow(s, 5) - 290*pow(t, 2) + 3003*s*pow(t, 2) - 7236*pow(s, 2)*pow(t, 2) +
4968*pow(s, 3)*pow(t, 2) + 465*pow(t, 3) - 3132*s*pow(t, 3) + 3672*pow(s,
2)*pow(t, 3) - 360*pow(t, 4) + 1188*s*pow(t, 4) + 108*pow(t, 5)))/25.;


  dphi(36,1) = (-648*s*t*(-1 + s + t)*(-20 + 174*s + 271*t - 1817*s*t - 580*pow(s,
2) + 4398*t*pow(s, 2) + 930*pow(s, 3) - 4572*t*pow(s, 3) - 720*pow(s, 4) +
1728*t*pow(s, 4) + 216*pow(s, 5) - 1237*pow(t, 2) + 6006*s*pow(t, 2) -
9396*pow(s, 2)*pow(t, 2) + 4752*pow(s, 3)*pow(t, 2) + 2538*pow(t, 3) -
7956*s*pow(t, 3) + 6048*pow(s, 2)*pow(t, 3) - 2412*pow(t, 4) + 3672*s*pow(t, 4)
+ 864*pow(t, 5)))/25.; 

  dphi(37,0) = (-648*s*t*(-1 + s + t)*(-20 + 271*s + 174*t
- 1817*s*t - 1237*pow(s, 2) + 6006*t*pow(s, 2) + 2538*pow(s, 3) - 7956*t*pow(s,
3) - 2412*pow(s, 4) + 3672*t*pow(s, 4) + 864*pow(s, 5) - 580*pow(t, 2) +
4398*s*pow(t, 2) - 9396*pow(s, 2)*pow(t, 2) + 6048*pow(s, 3)*pow(t, 2) +
930*pow(t, 3) - 4572*s*pow(t, 3) + 4752*pow(s, 2)*pow(t, 3) - 720*pow(t, 4) +
1728*s*pow(t, 4) + 216*pow(t, 5)))/25.; 

  dphi(37,1) = (-648*(-1 + s + t)*pow(s,
2)*(-10 + 87*s + 184*t - 1237*s*t - 290*pow(s, 2) + 3003*t*pow(s, 2) +
465*pow(s, 3) - 3132*t*pow(s, 3) - 360*pow(s, 4) + 1188*t*pow(s, 4) + 108*pow(s,
5) - 947*pow(t, 2) + 4611*s*pow(t, 2) - 7236*pow(s, 2)*pow(t, 2) + 3672*pow(s,
3)*pow(t, 2) + 2073*pow(t, 3) - 6516*s*pow(t, 3) + 4968*pow(s, 2)*pow(t, 3) -
2052*pow(t, 4) + 3132*s*pow(t, 4) + 756*pow(t, 5)))/25.; 

  dphi(38,0) =
(81*s*t*(-1 + s + t)*(-20 + 361*s + 114*t - 1787*s*t - 1951*pow(s, 2) +
7434*t*pow(s, 2) + 4464*pow(s, 3) - 11412*t*pow(s, 3) - 4572*pow(s, 4) +
5832*t*pow(s, 4) + 1728*pow(s, 5) - 238*pow(t, 2) + 3186*s*pow(t, 2) -
9180*pow(s, 2)*pow(t, 2) + 7128*pow(s, 3)*pow(t, 2) + 216*pow(t, 3) -
2412*s*pow(t, 3) + 3672*pow(s, 2)*pow(t, 3) - 72*pow(t, 4) + 648*s*pow(t,
4)))/4.; 

  dphi(38,1) = (81*(-1 + 6*s)*(-1 + s + t)*pow(s, 2)*(10 - 57*s - 124*t +
523*s*t + 119*pow(s, 2) - 720*t*pow(s, 2) - 108*pow(s, 3) + 324*t*pow(s, 3) +
36*pow(s, 4) + 404*pow(t, 2) - 1116*s*pow(t, 2) + 756*pow(s, 2)*pow(t, 2) -
504*pow(t, 3) + 684*s*pow(t, 3) + 216*pow(t, 4)))/4.; 

  dphi(39,0) = (81*s*t*(-1 +
s + t)*(-10 + 203*s + 22*t - 405*s*t - 1289*pow(s, 2) + 2178*t*pow(s, 2) +
3456*pow(s, 3) - 4356*t*pow(s, 3) - 4068*pow(s, 4) + 2808*t*pow(s, 4) +
1728*pow(s, 5) - 12*pow(t, 2) + 198*s*pow(t, 2) - 864*pow(s, 2)*pow(t, 2) +
1080*pow(s, 3)*pow(t, 2)))/4.; 

  dphi(39,1) = (81*(-1 + 2*s)*(-1 + 3*s)*(-1 +
6*s)*(-1 + s + t)*pow(s, 2)*(5 - 11*s - 27*t + 30*s*t + 6*pow(s, 2) + 24*pow(t,
2)))/4.; 

  dphi(40,0) = (-648*s*t*(-1 + s + t)*(-4 + 83*s + 4*t - 75*s*t -
545*pow(s, 2) + 420*t*pow(s, 2) + 1530*pow(s, 3) - 900*t*pow(s, 3) - 1908*pow(s,
4) + 648*t*pow(s, 4) + 864*pow(s, 5)))/25.; 

  dphi(40,1) = (-648*(-1 + 2*s)*(-2 +
3*s)*(-1 + 3*s)*(-1 + 6*s)*(-1 + s + t)*(-1 + s + 3*t)*pow(s, 2))/25.;


  dphi(41,0) = (324*s*(-4 + 81*s + 4*t - 75*s*t - 520*pow(s, 2) + 420*t*pow(s, 2)
+ 1425*pow(s, 3) - 900*t*pow(s, 3) - 1728*pow(s, 4) + 648*t*pow(s, 4) +
756*pow(s, 5))*pow(t, 2)*sqrt(2))/25.; 

  dphi(41,1) = (324*(-1 + 2*s)*(-2 +
3*s)*(-1 + 3*s)*(-1 + 6*s)*t*(-2 + 2*s + 3*t)*pow(s, 2)*sqrt(2))/25.; 

  dphi(42,0)
= (81*s*(-1 + 6*t)*(2 - 36*s - 2*t + 33*s*t + 188*pow(s, 2) - 144*t*pow(s, 2) -
360*pow(s, 3) + 180*t*pow(s, 3) + 216*pow(s, 4))*pow(t, 2)*pow(sqrt(2), -1))/4.;


  dphi(42,1) = (81*(-1 + 2*s)*(-1 + 3*s)*(-1 + 6*s)*t*pow(s, 2)*(2 - 2*s - 21*t +
18*s*t + 24*pow(t, 2))*pow(sqrt(2), -1))/4.; 

  dphi(43,0) = (81*s*(-1 + 2*t)*(-1 +
3*t)*(-1 + 6*t)*(2 - 21*s - 2*t + 18*s*t + 24*pow(s, 2))*pow(t, 2)*pow(sqrt(2),
-1))/4.; 

  dphi(43,1) = (81*(-1 + 6*s)*t*pow(s, 2)*(2 - 2*s - 36*t + 33*s*t +
188*pow(t, 2) - 144*s*pow(t, 2) - 360*pow(t, 3) + 180*s*pow(t, 3) + 216*pow(t,
4))*pow(sqrt(2), -1))/4.; 

  dphi(44,0) = (324*s*(-1 + 2*t)*(-2 + 3*s + 2*t)*(-2 +
3*t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*sqrt(2))/25.; 

  dphi(44,1) = (324*t*pow(s,
2)*(-4 + 4*s + 81*t - 75*s*t - 520*pow(t, 2) + 420*s*pow(t, 2) + 1425*pow(t, 3)
- 900*s*pow(t, 3) - 1728*pow(t, 4) + 648*s*pow(t, 4) + 756*pow(t,
5))*sqrt(2))/25.; 

  dphi(45,0) = 5832*s*(-1 + s + t)*(-1 + 2*s + t)*(-1 + 2*t)*(-1
+ 3*t)*(-1 + 6*t)*pow(t, 2); 

  dphi(45,1) = 2916*t*(-1 + s + t)*pow(s, 2)*(2 - 2*s
- 37*t + 33*s*t + 199*pow(t, 2) - 144*s*pow(t, 2) - 396*pow(t, 3) + 180*s*pow(t,
3) + 252*pow(t, 4)); 

  dphi(46,0) = -2592*s*(-1 + s + t)*(-1 + 3*t)*(-1 +
6*t)*pow(t, 2)*(5 - 19*s - 11*t + 21*s*t + 15*pow(s, 2) + 6*pow(t, 2));


  dphi(46,1) = -1296*t*(-1 + s + t)*pow(s, 2)*(10 - 22*s - 173*t + 339*s*t +
12*pow(s, 2) - 162*t*pow(s, 2) + 831*pow(t, 2) - 1278*s*pow(t, 2) + 432*pow(s,
2)*pow(t, 2) - 1404*pow(t, 3) + 1188*s*pow(t, 3) + 756*pow(t, 4)); 

  dphi(47,0) =
1296*s*(-1 + s + t)*(-1 + 6*t)*pow(t, 2)*(-20 + 121*s + 74*t - 297*s*t -
207*pow(s, 2) + 252*t*pow(s, 2) + 108*pow(s, 3) - 90*pow(t, 2) + 180*s*pow(t, 2)
+ 36*pow(t, 3)); 

  dphi(47,1) = 1296*t*(-1 + s + t)*pow(s, 2)*(20 - 74*s - 301*t +
963*s*t + 90*pow(s, 2) - 990*t*pow(s, 2) - 36*pow(s, 3) + 324*t*pow(s, 3) +
1155*pow(t, 2) - 2574*s*pow(t, 2) + 1404*pow(s, 2)*pow(t, 2) - 1620*pow(t, 3) +
1836*s*pow(t, 3) + 756*pow(t, 4)); 

  dphi(48,0) = -2916*s*(-1 + s + t)*pow(t,
2)*(20 - 181*s - 114*t + 761*s*t + 523*pow(s, 2) - 1440*t*pow(s, 2) - 612*pow(s,
3) + 828*t*pow(s, 3) + 252*pow(s, 4) + 238*pow(t, 2) - 1044*s*pow(t, 2) +
972*pow(s, 2)*pow(t, 2) - 216*pow(t, 3) + 468*s*pow(t, 3) + 72*pow(t, 4));


  dphi(48,1) = -2916*t*(-1 + s + t)*pow(s, 2)*(20 - 114*s - 181*t + 761*s*t +
238*pow(s, 2) - 1044*t*pow(s, 2) - 216*pow(s, 3) + 468*t*pow(s, 3) + 72*pow(s,
4) + 523*pow(t, 2) - 1440*s*pow(t, 2) + 972*pow(s, 2)*pow(t, 2) - 612*pow(t, 3)
+ 828*s*pow(t, 3) + 252*pow(t, 4)); 

  dphi(49,0) = 1296*s*(-1 + s + t)*pow(t,
2)*(20 - 301*s - 74*t + 963*s*t + 1155*pow(s, 2) - 2574*t*pow(s, 2) -
1620*pow(s, 3) + 1836*t*pow(s, 3) + 756*pow(s, 4) + 90*pow(t, 2) - 990*s*pow(t,
2) + 1404*pow(s, 2)*pow(t, 2) - 36*pow(t, 3) + 324*s*pow(t, 3)); 

  dphi(49,1) =
1296*(-1 + 6*s)*t*(-1 + s + t)*pow(s, 2)*(-20 + 74*s + 121*t - 297*s*t -
90*pow(s, 2) + 180*t*pow(s, 2) + 36*pow(s, 3) - 207*pow(t, 2) + 252*s*pow(t, 2)
+ 108*pow(t, 3)); 

  dphi(50,0) = -1296*s*(-1 + s + t)*pow(t, 2)*(10 - 173*s - 22*t
+ 339*s*t + 831*pow(s, 2) - 1278*t*pow(s, 2) - 1404*pow(s, 3) + 1188*t*pow(s, 3)
+ 756*pow(s, 4) + 12*pow(t, 2) - 162*s*pow(t, 2) + 432*pow(s, 2)*pow(t, 2));


  dphi(50,1) = -2592*(-1 + 3*s)*(-1 + 6*s)*t*(-1 + s + t)*pow(s, 2)*(5 - 11*s -
19*t + 21*s*t + 6*pow(s, 2) + 15*pow(t, 2)); 

  dphi(51,0) = 2916*s*(-1 + s + t)*(2
- 37*s - 2*t + 33*s*t + 199*pow(s, 2) - 144*t*pow(s, 2) - 396*pow(s, 3) +
180*t*pow(s, 3) + 252*pow(s, 4))*pow(t, 2); 

  dphi(51,1) = 5832*(-1 + 2*s)*(-1 +
3*s)*(-1 + 6*s)*t*(-1 + s + t)*(-1 + s + 2*t)*pow(s, 2); 

  dphi(52,0) = 1296*s*(-1
+ s + t)*(-1 + 6*t)*(-2 + 31*s + 2*t - 27*s*t - 117*pow(s, 2) + 72*t*pow(s, 2) +
108*pow(s, 3))*pow(t, 2); 

  dphi(52,1) = 2592*(-1 + 3*s)*(-1 + 6*s)*t*(-1 + s +
t)*pow(s, 2)*(1 - s - 11*t + 9*s*t + 15*pow(t, 2)); 

  dphi(53,0) = 2592*s*(-1 + s
+ t)*(-1 + 3*t)*(-1 + 6*t)*(1 - 11*s - t + 9*s*t + 15*pow(s, 2))*pow(t, 2);


  dphi(53,1) = 1296*(-1 + 6*s)*t*(-1 + s + t)*pow(s, 2)*(-2 + 2*s + 31*t - 27*s*t
- 117*pow(t, 2) + 72*s*pow(t, 2) + 108*pow(t, 3)); 

  dphi(54,0) = -1458*s*(-1 + s
  + t)*(-1 + 2*s + t)*(-1 + 6*t)*(5 - 54*s - 6*t + 54*s*t + 54*pow(s, 2))*pow(t,
2); 

  dphi(54,1) = -1458*(-1 + 6*s)*t*(-1 + s + t)*(-1 + s + 2*t)*pow(s, 2)*(5 -
6*s - 54*t + 54*s*t + 54*pow(t, 2)); 
}

/// Precomputed d2basis polynomials for the 5th order boundary representation
template <> void  BernadouElementBasis<5>::d2full_basic_polynomials(const Vector<double>&
s_basic, DShape&d2phi) const 
{
  // For convenience 
  const double s=s_basic[0], t=s_basic[1];
  // Now fill in (automatically generated code)
  d2phi(0,0) = (3*(598360*s - 16261464*pow(s, 2) + 133811000*pow(s, 3) -
474346260*pow(s, 4) + 813567132*pow(s, 5) - 665594496*pow(s, 6) +
208225728*pow(s, 7) - 3415480*pow(t, 2) + 126762915*s*pow(t, 2) -
1077009696*pow(s, 2)*pow(t, 2) + 3315149910*pow(s, 3)*pow(t, 2) -
4117183020*pow(s, 4)*pow(t, 2) + 1756160028*pow(s, 5)*pow(t, 2) + 7389313*pow(t,
3) - 247491396*s*pow(t, 3) + 1848699990*pow(s, 2)*pow(t, 3) - 4323410640*pow(s,
3)*pow(t, 3) + 2960191980*pow(s, 4)*pow(t, 3) - 4195422*pow(t, 4) +
92584809*s*pow(t, 4) - 685821816*pow(s, 2)*pow(t, 4) + 1092981600*pow(s,
3)*pow(t, 4) + 1836081*pow(t, 5) + 42791328*s*pow(t, 5) - 67530672*pow(s,
2)*pow(t, 5) - 3266568*pow(t, 6) - 15545196*s*pow(t, 6) + 1652076*pow(t,
7)))/2000.; 

  d2phi(0,1) = (3*s*t*(-6830960 + 126762915*s + 22167939*t -
371237094*s*t - 718006464*pow(s, 2) + 1848699990*t*pow(s, 2) + 1657574955*pow(s,
3) - 3242557980*t*pow(s, 3) - 1646873208*pow(s, 4) + 1776115188*t*pow(s, 4) +
585386676*pow(s, 5) - 16781688*pow(t, 2) + 185169618*s*pow(t, 2) -
914429088*pow(s, 2)*pow(t, 2) + 1092981600*pow(s, 3)*pow(t, 2) + 9180405*pow(t,
3) + 106978320*s*pow(t, 3) - 112551120*pow(s, 2)*pow(t, 3) - 19599408*pow(t, 4)
- 46635588*s*pow(t, 4) + 11564532*pow(t, 5)))/2000.; 

  d2phi(0,2) = (3*pow(s,
  2)*(-3415480 + 42254305*s + 22167939*t - 247491396*s*t - 179501616*pow(s, 2) +
924349995*t*pow(s, 2) + 331514991*pow(s, 3) - 1297023192*t*pow(s, 3) -
274478868*pow(s, 4) + 592038396*t*pow(s, 4) + 83626668*pow(s, 5) -
25172532*pow(t, 2) + 185169618*s*pow(t, 2) - 685821816*pow(s, 2)*pow(t, 2) +
655788960*pow(s, 3)*pow(t, 2) + 18360810*pow(t, 3) + 142637760*s*pow(t, 3) -
112551120*pow(s, 2)*pow(t, 3) - 48998520*pow(t, 4) - 77725980*s*pow(t, 4) +
34693596*pow(t, 5)))/2000.; 

  d2phi(1,0) = (3*pow(t, 2)*(-3415480 + 22167939*s +
42254305*t - 247491396*s*t - 25172532*pow(s, 2) + 185169618*t*pow(s, 2) +
18360810*pow(s, 3) + 142637760*t*pow(s, 3) - 48998520*pow(s, 4) -
77725980*t*pow(s, 4) + 34693596*pow(s, 5) - 179501616*pow(t, 2) +
924349995*s*pow(t, 2) - 685821816*pow(s, 2)*pow(t, 2) - 112551120*pow(s,
3)*pow(t, 2) + 331514991*pow(t, 3) - 1297023192*s*pow(t, 3) + 655788960*pow(s,
2)*pow(t, 3) - 274478868*pow(t, 4) + 592038396*s*pow(t, 4) + 83626668*pow(t,
5)))/2000.; 

  d2phi(1,1) = (3*s*t*(-6830960 + 22167939*s + 126762915*t -
371237094*s*t - 16781688*pow(s, 2) + 185169618*t*pow(s, 2) + 9180405*pow(s, 3) +
106978320*t*pow(s, 3) - 19599408*pow(s, 4) - 46635588*t*pow(s, 4) +
11564532*pow(s, 5) - 718006464*pow(t, 2) + 1848699990*s*pow(t, 2) -
914429088*pow(s, 2)*pow(t, 2) - 112551120*pow(s, 3)*pow(t, 2) +
1657574955*pow(t, 3) - 3242557980*s*pow(t, 3) + 1092981600*pow(s, 2)*pow(t, 3) -
1646873208*pow(t, 4) + 1776115188*s*pow(t, 4) + 585386676*pow(t, 5)))/2000.;


  d2phi(1,2) = (3*(598360*t - 3415480*pow(s, 2) + 126762915*t*pow(s, 2) +
7389313*pow(s, 3) - 247491396*t*pow(s, 3) - 4195422*pow(s, 4) +
92584809*t*pow(s, 4) + 1836081*pow(s, 5) + 42791328*t*pow(s, 5) - 3266568*pow(s,
6) - 15545196*t*pow(s, 6) + 1652076*pow(s, 7) - 16261464*pow(t, 2) -
1077009696*pow(s, 2)*pow(t, 2) + 1848699990*pow(s, 3)*pow(t, 2) -
685821816*pow(s, 4)*pow(t, 2) - 67530672*pow(s, 5)*pow(t, 2) + 133811000*pow(t,
3) + 3315149910*pow(s, 2)*pow(t, 3) - 4323410640*pow(s, 3)*pow(t, 3) +
1092981600*pow(s, 4)*pow(t, 3) - 474346260*pow(t, 4) - 4117183020*pow(s,
2)*pow(t, 4) + 2960191980*pow(s, 3)*pow(t, 4) + 813567132*pow(t, 5) +
1756160028*pow(s, 2)*pow(t, 5) - 665594496*pow(t, 6) + 208225728*pow(t,
7)))/2000.; 

  d2phi(2,0) = (-3*(1986086*s - 31794072*pow(s, 2) + 174053920*pow(s,
3) - 448736220*pow(s, 4) + 596370222*pow(s, 5) - 395992800*pow(s, 6) +
104112864*pow(s, 7) - 5622406*pow(t, 2) + 81036633*s*pow(t, 2) -
205529076*pow(s, 2)*pow(t, 2) - 90608490*pow(s, 3)*pow(t, 2) + 618729030*pow(s,
4)*pow(t, 2) - 405776952*pow(s, 5)*pow(t, 2) + 27012211*pow(t, 3) -
435713364*s*pow(t, 3) + 1363972338*pow(s, 2)*pow(t, 3) - 1220397840*pow(s,
3)*pow(t, 3) + 181127340*pow(s, 4)*pow(t, 3) - 34254846*pow(t, 4) +
681986169*s*pow(t, 4) - 1789749432*pow(s, 2)*pow(t, 4) + 1068108120*pow(s,
3)*pow(t, 4) - 9060849*pow(t, 5) - 366119352*s*pow(t, 5) + 640864872*pow(s,
2)*pow(t, 5) + 41248602*pow(t, 6) + 36225468*s*pow(t, 6) - 19322712*pow(t,
7)))/1000.; 

  d2phi(2,1) = (3*s*t*(11244812 - 81036633*s - 81036633*t +
653570046*s*t + 137019384*pow(s, 2) - 1363972338*t*pow(s, 2) + 45304245*pow(s,
3) + 915298380*t*pow(s, 3) - 247491612*pow(s, 4) - 108676404*t*pow(s, 4) +
135258984*pow(s, 5) + 137019384*pow(t, 2) - 1363972338*s*pow(t, 2) +
2386332576*pow(s, 2)*pow(t, 2) - 1068108120*pow(s, 3)*pow(t, 2) +
45304245*pow(t, 3) + 915298380*s*pow(t, 3) - 1068108120*pow(s, 2)*pow(t, 3) -
247491612*pow(t, 4) - 108676404*s*pow(t, 4) + 135258984*pow(t, 5)))/1000.;


  d2phi(2,2) = (3*(-1986086*t + 5622406*pow(s, 2) - 81036633*t*pow(s, 2) -
27012211*pow(s, 3) + 435713364*t*pow(s, 3) + 34254846*pow(s, 4) -
681986169*t*pow(s, 4) + 9060849*pow(s, 5) + 366119352*t*pow(s, 5) -
41248602*pow(s, 6) - 36225468*t*pow(s, 6) + 19322712*pow(s, 7) + 31794072*pow(t,
2) + 205529076*pow(s, 2)*pow(t, 2) - 1363972338*pow(s, 3)*pow(t, 2) +
1789749432*pow(s, 4)*pow(t, 2) - 640864872*pow(s, 5)*pow(t, 2) -
174053920*pow(t, 3) + 90608490*pow(s, 2)*pow(t, 3) + 1220397840*pow(s, 3)*pow(t,
3) - 1068108120*pow(s, 4)*pow(t, 3) + 448736220*pow(t, 4) - 618729030*pow(s,
2)*pow(t, 4) - 181127340*pow(s, 3)*pow(t, 4) - 596370222*pow(t, 5) +
405776952*pow(s, 2)*pow(t, 5) + 395992800*pow(t, 6) - 104112864*pow(t,
7)))/1000.; 

  d2phi(3,0) = (-18840*s + 513816*pow(s, 2) - 4251000*pow(s, 3) +
15179940*pow(s, 4) - 26265708*pow(s, 5) + 21700224*pow(s, 6) - 6858432*pow(s, 7)
+ 103738*pow(t, 2) - 3854289*s*pow(t, 2) + 32818824*pow(s, 2)*pow(t, 2) -
101369790*pow(s, 3)*pow(t, 2) + 126520380*pow(s, 4)*pow(t, 2) - 54316332*pow(s,
5)*pow(t, 2) - 222565*pow(t, 3) + 7445574*s*pow(t, 3) - 55666710*pow(s,
2)*pow(t, 3) + 130358160*pow(s, 3)*pow(t, 3) - 89506620*pow(s, 4)*pow(t, 3) +
127188*pow(t, 4) - 2759481*s*pow(t, 4) + 20491704*pow(s, 2)*pow(t, 4) -
32594400*pow(s, 3)*pow(t, 4) - 59229*pow(t, 5) - 1288872*s*pow(t, 5) +
1986768*pow(s, 2)*pow(t, 5) + 102384*pow(t, 6) + 475308*s*pow(t, 6) -
51516*pow(t, 7))/200.; 

  d2phi(3,1) = -(s*t*(-207476 + 3854289*s + 667695*t -
11168361*s*t - 21879216*pow(s, 2) + 55666710*t*pow(s, 2) + 50684895*pow(s, 3) -
97768620*t*pow(s, 3) - 50608152*pow(s, 4) + 53703972*t*pow(s, 4) +
18105444*pow(s, 5) - 508752*pow(t, 2) + 5518962*s*pow(t, 2) - 27322272*pow(s,
2)*pow(t, 2) + 32594400*pow(s, 3)*pow(t, 2) + 296145*pow(t, 3) +
3222180*s*pow(t, 3) - 3311280*pow(s, 2)*pow(t, 3) - 614304*pow(t, 4) -
1425924*s*pow(t, 4) + 360612*pow(t, 5)))/200.; 

  d2phi(3,2) = -(pow(s, 2)*(-103738
+ 1284763*s + 667695*t - 7445574*s*t - 5469804*pow(s, 2) + 27833355*t*pow(s, 2)
+ 10136979*pow(s, 3) - 39107448*t*pow(s, 3) - 8434692*pow(s, 4) +
17901324*t*pow(s, 4) + 2586492*pow(s, 5) - 763128*pow(t, 2) + 5518962*s*pow(t,
2) - 20491704*pow(s, 2)*pow(t, 2) + 19556640*pow(s, 3)*pow(t, 2) + 592290*pow(t,
3) + 4296240*s*pow(t, 3) - 3311280*pow(s, 2)*pow(t, 3) - 1535760*pow(t, 4) -
2376540*s*pow(t, 4) + 1081836*pow(t, 5)))/200.; 

  d2phi(4,0) = (t*(-6680 +
293388*s + 9488*t - 345093*s*t - 3221616*pow(s, 2) + 2513700*t*pow(s, 2) +
14459400*pow(s, 3) - 6876270*t*pow(s, 3) - 30645000*pow(s, 4) + 9020160*t*pow(s,
4) + 30454704*pow(s, 5) - 4769604*t*pow(s, 5) - 11394432*pow(s, 6) -
55323*pow(t, 2) + 1915110*s*pow(t, 2) - 9840042*pow(s, 2)*pow(t, 2) +
11945880*pow(s, 3)*pow(t, 2) - 3018060*pow(s, 4)*pow(t, 2) + 122310*pow(t, 3) -
4806621*s*pow(t, 3) + 21387888*pow(s, 2)*pow(t, 3) - 17184960*pow(s, 3)*pow(t,
3) - 18927*pow(t, 4) + 3399084*s*pow(t, 4) - 10427616*pow(s, 2)*pow(t, 4) -
102384*pow(t, 5) - 475308*s*pow(t, 5) + 51516*pow(t, 6)))/200.; 

  d2phi(4,1) =
-(s*(6680 - 146694*s - 18976*t + 345093*s*t + 1073872*pow(s, 2) -
1675800*t*pow(s, 2) - 3614850*pow(s, 3) + 3438135*t*pow(s, 3) + 6129000*pow(s,
4) - 3608064*t*pow(s, 4) - 5075784*pow(s, 5) + 1589868*t*pow(s, 5) +
1627776*pow(s, 6) + 165969*pow(t, 2) - 2872665*s*pow(t, 2) + 9840042*pow(s,
2)*pow(t, 2) - 8959410*pow(s, 3)*pow(t, 2) + 1810836*pow(s, 4)*pow(t, 2) -
489240*pow(t, 3) + 9613242*s*pow(t, 3) - 28517184*pow(s, 2)*pow(t, 3) +
17184960*pow(s, 3)*pow(t, 3) + 94635*pow(t, 4) - 8497710*s*pow(t, 4) +
17379360*pow(s, 2)*pow(t, 4) + 614304*pow(t, 5) + 1425924*s*pow(t, 5) -
360612*pow(t, 6)))/200.; 

  d2phi(4,2) = -(pow(s, 2)*(-9488 + 115031*s + 165969*t -
1915110*s*t - 418950*pow(s, 2) + 4920021*t*pow(s, 2) + 687627*pow(s, 3) -
3583764*t*pow(s, 3) - 601344*pow(s, 4) + 603612*t*pow(s, 4) + 227124*pow(s, 5) -
733860*pow(t, 2) + 9613242*s*pow(t, 2) - 21387888*pow(s, 2)*pow(t, 2) +
10310976*pow(s, 3)*pow(t, 2) + 189270*pow(t, 3) - 11330280*s*pow(t, 3) +
17379360*pow(s, 2)*pow(t, 3) + 1535760*pow(t, 4) + 2376540*s*pow(t, 4) -
1081836*pow(t, 5)))/200.; 

  d2phi(5,0) = -(pow(t, 2)*(-9488 + 165969*s + 115031*t
- 1915110*s*t - 733860*pow(s, 2) + 9613242*t*pow(s, 2) + 189270*pow(s, 3) -
11330280*t*pow(s, 3) + 1535760*pow(s, 4) + 2376540*t*pow(s, 4) - 1081836*pow(s,
5) - 418950*pow(t, 2) + 4920021*s*pow(t, 2) - 21387888*pow(s, 2)*pow(t, 2) +
17379360*pow(s, 3)*pow(t, 2) + 687627*pow(t, 3) - 3583764*s*pow(t, 3) +
10310976*pow(s, 2)*pow(t, 3) - 601344*pow(t, 4) + 603612*s*pow(t, 4) +
227124*pow(t, 5)))/200.; 

  d2phi(5,1) = -(t*(6680 - 18976*s - 146694*t +
345093*s*t + 165969*pow(s, 2) - 2872665*t*pow(s, 2) - 489240*pow(s, 3) +
9613242*t*pow(s, 3) + 94635*pow(s, 4) - 8497710*t*pow(s, 4) + 614304*pow(s, 5) +
1425924*t*pow(s, 5) - 360612*pow(s, 6) + 1073872*pow(t, 2) - 1675800*s*pow(t, 2)
+ 9840042*pow(s, 2)*pow(t, 2) - 28517184*pow(s, 3)*pow(t, 2) + 17379360*pow(s,
4)*pow(t, 2) - 3614850*pow(t, 3) + 3438135*s*pow(t, 3) - 8959410*pow(s,
2)*pow(t, 3) + 17184960*pow(s, 3)*pow(t, 3) + 6129000*pow(t, 4) -
3608064*s*pow(t, 4) + 1810836*pow(s, 2)*pow(t, 4) - 5075784*pow(t, 5) +
1589868*s*pow(t, 5) + 1627776*pow(t, 6)))/200.; 

  d2phi(5,2) = (s*(-6680 + 9488*s
+ 293388*t - 345093*s*t - 55323*pow(s, 2) + 1915110*t*pow(s, 2) + 122310*pow(s,
3) - 4806621*t*pow(s, 3) - 18927*pow(s, 4) + 3399084*t*pow(s, 4) - 102384*pow(s,
5) - 475308*t*pow(s, 5) + 51516*pow(s, 6) - 3221616*pow(t, 2) + 2513700*s*pow(t,
2) - 9840042*pow(s, 2)*pow(t, 2) + 21387888*pow(s, 3)*pow(t, 2) -
10427616*pow(s, 4)*pow(t, 2) + 14459400*pow(t, 3) - 6876270*s*pow(t, 3) +
11945880*pow(s, 2)*pow(t, 3) - 17184960*pow(s, 3)*pow(t, 3) - 30645000*pow(t, 4)
+ 9020160*s*pow(t, 4) - 3018060*pow(s, 2)*pow(t, 4) + 30454704*pow(t, 5) -
4769604*s*pow(t, 5) - 11394432*pow(t, 6)))/200.; 

  d2phi(6,0) = -(pow(t,
2)*(-103738 + 667695*s + 1284763*t - 7445574*s*t - 763128*pow(s, 2) +
5518962*t*pow(s, 2) + 592290*pow(s, 3) + 4296240*t*pow(s, 3) - 1535760*pow(s, 4)
- 2376540*t*pow(s, 4) + 1081836*pow(s, 5) - 5469804*pow(t, 2) +
27833355*s*pow(t, 2) - 20491704*pow(s, 2)*pow(t, 2) - 3311280*pow(s, 3)*pow(t,
2) + 10136979*pow(t, 3) - 39107448*s*pow(t, 3) + 19556640*pow(s, 2)*pow(t, 3) -
8434692*pow(t, 4) + 17901324*s*pow(t, 4) + 2586492*pow(t, 5)))/200.; 

  d2phi(6,1)
= -(s*t*(-207476 + 667695*s + 3854289*t - 11168361*s*t - 508752*pow(s, 2) +
5518962*t*pow(s, 2) + 296145*pow(s, 3) + 3222180*t*pow(s, 3) - 614304*pow(s, 4)
- 1425924*t*pow(s, 4) + 360612*pow(s, 5) - 21879216*pow(t, 2) +
55666710*s*pow(t, 2) - 27322272*pow(s, 2)*pow(t, 2) - 3311280*pow(s, 3)*pow(t,
2) + 50684895*pow(t, 3) - 97768620*s*pow(t, 3) + 32594400*pow(s, 2)*pow(t, 3) -
50608152*pow(t, 4) + 53703972*s*pow(t, 4) + 18105444*pow(t, 5)))/200.;


  d2phi(6,2) = (-18840*t + 103738*pow(s, 2) - 3854289*t*pow(s, 2) - 222565*pow(s,
3) + 7445574*t*pow(s, 3) + 127188*pow(s, 4) - 2759481*t*pow(s, 4) - 59229*pow(s,
5) - 1288872*t*pow(s, 5) + 102384*pow(s, 6) + 475308*t*pow(s, 6) - 51516*pow(s,
7) + 513816*pow(t, 2) + 32818824*pow(s, 2)*pow(t, 2) - 55666710*pow(s, 3)*pow(t,
2) + 20491704*pow(s, 4)*pow(t, 2) + 1986768*pow(s, 5)*pow(t, 2) - 4251000*pow(t,
3) - 101369790*pow(s, 2)*pow(t, 3) + 130358160*pow(s, 3)*pow(t, 3) -
32594400*pow(s, 4)*pow(t, 3) + 15179940*pow(t, 4) + 126520380*pow(s, 2)*pow(t,
4) - 89506620*pow(s, 3)*pow(t, 4) - 26265708*pow(t, 5) - 54316332*pow(s,
2)*pow(t, 5) + 21700224*pow(t, 6) - 6858432*pow(t, 7))/200.; 

  d2phi(7,0) =
(-80334*s + 1170168*pow(s, 2) - 6114480*pow(s, 3) + 15345180*pow(s, 4) -
20045718*pow(s, 5) + 13154400*pow(s, 6) - 3429216*pow(s, 7) + 276743*pow(t, 2) -
2233623*s*pow(t, 2) + 1352538*pow(s, 2)*pow(t, 2) + 16665030*pow(s, 3)*pow(t, 2)
- 33653070*pow(s, 4)*pow(t, 2) + 17840088*pow(s, 5)*pow(t, 2) - 1963172*pow(t,
3) + 16216308*s*pow(t, 3) - 33488424*pow(s, 2)*pow(t, 3) + 13087980*pow(s,
3)*pow(t, 3) + 8937540*pow(s, 4)*pow(t, 3) + 5469291*pow(t, 4) -
36284409*s*pow(t, 4) + 60902604*pow(s, 2)*pow(t, 4) - 25654320*pow(s, 3)*pow(t,
4) - 7495902*pow(t, 5) + 33674454*s*pow(t, 5) - 30188376*pow(s, 2)*pow(t, 5) +
5059584*pow(t, 6) - 11302416*s*pow(t, 6) - 1346544*pow(t, 7))/100.; 

  d2phi(7,1) =
-(t*(30118 - 553486*s - 361179*t + 5889516*s*t + 2233623*pow(s, 2) -
24324462*t*pow(s, 2) - 901692*pow(s, 3) + 33488424*t*pow(s, 3) - 8332515*pow(s,
4) - 9815985*t*pow(s, 4) + 13461228*pow(s, 5) - 5362524*t*pow(s, 5) -
5946696*pow(s, 6) + 1680476*pow(t, 2) - 21877164*s*pow(t, 2) + 72568818*pow(s,
2)*pow(t, 2) - 81203472*pow(s, 3)*pow(t, 2) + 25654320*pow(s, 4)*pow(t, 2) -
3932775*pow(t, 3) + 37479510*s*pow(t, 3) - 84186135*pow(s, 2)*pow(t, 3) +
50313960*pow(s, 3)*pow(t, 3) + 4928796*pow(t, 4) - 30357504*s*pow(t, 4) +
33907248*pow(s, 2)*pow(t, 4) - 3159324*pow(t, 5) + 9425808*s*pow(t, 5) +
813888*pow(t, 6)))/100.; 

  d2phi(7,2) = (s*(-30118 + 276743*s + 722358*t -
5889516*s*t - 744541*pow(s, 2) + 16216308*t*pow(s, 2) + 225423*pow(s, 3) -
16744212*t*pow(s, 3) + 1666503*pow(s, 4) + 3926394*t*pow(s, 4) - 2243538*pow(s,
5) + 1787508*t*pow(s, 5) + 849528*pow(s, 6) - 5041428*pow(t, 2) +
32815746*s*pow(t, 2) - 72568818*pow(s, 2)*pow(t, 2) + 60902604*pow(s, 3)*pow(t,
2) - 15392592*pow(s, 4)*pow(t, 2) + 15731100*pow(t, 3) - 74959020*s*pow(t, 3) +
112248180*pow(s, 2)*pow(t, 3) - 50313960*pow(s, 3)*pow(t, 3) - 24643980*pow(t,
4) + 75893760*s*pow(t, 4) - 56512080*pow(s, 2)*pow(t, 4) + 18955944*pow(t, 5) -
28277424*s*pow(t, 5) - 5697216*pow(t, 6)))/100.; 

  d2phi(8,0) = (t*(-30118 +
722358*s + 276743*t - 5889516*s*t - 5041428*pow(s, 2) + 32815746*t*pow(s, 2) +
15731100*pow(s, 3) - 74959020*t*pow(s, 3) - 24643980*pow(s, 4) +
75893760*t*pow(s, 4) + 18955944*pow(s, 5) - 28277424*t*pow(s, 5) -
5697216*pow(s, 6) - 744541*pow(t, 2) + 16216308*s*pow(t, 2) - 72568818*pow(s,
2)*pow(t, 2) + 112248180*pow(s, 3)*pow(t, 2) - 56512080*pow(s, 4)*pow(t, 2) +
225423*pow(t, 3) - 16744212*s*pow(t, 3) + 60902604*pow(s, 2)*pow(t, 3) -
50313960*pow(s, 3)*pow(t, 3) + 1666503*pow(t, 4) + 3926394*s*pow(t, 4) -
15392592*pow(s, 2)*pow(t, 4) - 2243538*pow(t, 5) + 1787508*s*pow(t, 5) +
849528*pow(t, 6)))/100.; 

  d2phi(8,1) = -(s*(30118 - 361179*s - 553486*t +
5889516*s*t + 1680476*pow(s, 2) - 21877164*t*pow(s, 2) - 3932775*pow(s, 3) +
37479510*t*pow(s, 3) + 4928796*pow(s, 4) - 30357504*t*pow(s, 4) - 3159324*pow(s,
5) + 9425808*t*pow(s, 5) + 813888*pow(s, 6) + 2233623*pow(t, 2) -
24324462*s*pow(t, 2) + 72568818*pow(s, 2)*pow(t, 2) - 84186135*pow(s, 3)*pow(t,
2) + 33907248*pow(s, 4)*pow(t, 2) - 901692*pow(t, 3) + 33488424*s*pow(t, 3) -
81203472*pow(s, 2)*pow(t, 3) + 50313960*pow(s, 3)*pow(t, 3) - 8332515*pow(t, 4)
- 9815985*s*pow(t, 4) + 25654320*pow(s, 2)*pow(t, 4) + 13461228*pow(t, 5) -
5362524*s*pow(t, 5) - 5946696*pow(t, 6)))/100.; 

  d2phi(8,2) = (-80334*t +
276743*pow(s, 2) - 2233623*t*pow(s, 2) - 1963172*pow(s, 3) + 16216308*t*pow(s,
3) + 5469291*pow(s, 4) - 36284409*t*pow(s, 4) - 7495902*pow(s, 5) +
33674454*t*pow(s, 5) + 5059584*pow(s, 6) - 11302416*t*pow(s, 6) - 1346544*pow(s,
7) + 1170168*pow(t, 2) + 1352538*pow(s, 2)*pow(t, 2) - 33488424*pow(s, 3)*pow(t,
2) + 60902604*pow(s, 4)*pow(t, 2) - 30188376*pow(s, 5)*pow(t, 2) -
6114480*pow(t, 3) + 16665030*pow(s, 2)*pow(t, 3) + 13087980*pow(s, 3)*pow(t, 3)
- 25654320*pow(s, 4)*pow(t, 3) + 15345180*pow(t, 4) - 33653070*pow(s, 2)*pow(t,
4) + 8937540*pow(s, 3)*pow(t, 4) - 20045718*pow(t, 5) + 17840088*pow(s,
2)*pow(t, 5) + 13154400*pow(t, 6) - 3429216*pow(t, 7))/100.; 

  d2phi(9,0) = (120*s
- 3288*pow(s, 2) + 27400*pow(s, 3) - 98820*pow(s, 4) + 173124*pow(s, 5) -
145152*pow(s, 6) + 46656*pow(s, 7) - 628*pow(t, 2) + 23349*s*pow(t, 2) -
199152*pow(s, 2)*pow(t, 2) + 616770*pow(s, 3)*pow(t, 2) - 772740*pow(s,
4)*pow(t, 2) + 333396*pow(s, 5)*pow(t, 2) + 1339*pow(t, 3) - 44712*s*pow(t, 3) +
334530*pow(s, 2)*pow(t, 3) - 784080*pow(s, 3)*pow(t, 3) + 539460*pow(s,
4)*pow(t, 3) - 774*pow(t, 4) + 16443*s*pow(t, 4) - 122472*pow(s, 2)*pow(t, 4) +
194400*pow(s, 3)*pow(t, 4) + 387*pow(t, 5) + 7776*s*pow(t, 5) - 11664*pow(s,
2)*pow(t, 5) - 648*pow(t, 6) - 2916*s*pow(t, 6) + 324*pow(t, 7))/40.; 

  d2phi(9,1)
= (s*t*(-1256 + 23349*s + 4017*t - 67068*s*t - 132768*pow(s, 2) +
334530*t*pow(s, 2) + 308385*pow(s, 3) - 588060*t*pow(s, 3) - 309096*pow(s, 4) +
323676*t*pow(s, 4) + 111132*pow(s, 5) - 3096*pow(t, 2) + 32886*s*pow(t, 2) -
163296*pow(s, 2)*pow(t, 2) + 194400*pow(s, 3)*pow(t, 2) + 1935*pow(t, 3) +
19440*s*pow(t, 3) - 19440*pow(s, 2)*pow(t, 3) - 3888*pow(t, 4) - 8748*s*pow(t,
4) + 2268*pow(t, 5)))/40.; 

  d2phi(9,2) = (pow(s, 2)*(-628 + 7783*s + 4017*t -
44712*s*t - 33192*pow(s, 2) + 167265*t*pow(s, 2) + 61677*pow(s, 3) -
235224*t*pow(s, 3) - 51516*pow(s, 4) + 107892*t*pow(s, 4) + 15876*pow(s, 5) -
4644*pow(t, 2) + 32886*s*pow(t, 2) - 122472*pow(s, 2)*pow(t, 2) + 116640*pow(s,
3)*pow(t, 2) + 3870*pow(t, 3) + 25920*s*pow(t, 3) - 19440*pow(s, 2)*pow(t, 3) -
9720*pow(t, 4) - 14580*s*pow(t, 4) + 6804*pow(t, 5)))/40.; 

  d2phi(10,0) =
-(t*(-40 + 1764*s + 50*t - 1797*s*t - 19488*pow(s, 2) + 12420*t*pow(s, 2) +
88200*pow(s, 3) - 31410*t*pow(s, 3) - 189000*pow(s, 4) + 38880*t*pow(s, 4) +
190512*pow(s, 5) - 20412*t*pow(s, 5) - 72576*pow(s, 6) - 325*pow(t, 2) +
11340*s*pow(t, 2) - 57726*pow(s, 2)*pow(t, 2) + 68040*pow(s, 3)*pow(t, 2) -
14580*pow(s, 4)*pow(t, 2) + 720*pow(t, 3) - 28863*s*pow(t, 3) + 128304*pow(s,
2)*pow(t, 3) - 103680*pow(s, 3)*pow(t, 3) - 81*pow(t, 4) + 20412*s*pow(t, 4) -
62208*pow(s, 2)*pow(t, 4) - 648*pow(t, 5) - 2916*s*pow(t, 5) + 324*pow(t,
6)))/20.; 

  d2phi(10,1) = (s*(40 - 882*s - 100*t + 1797*s*t + 6496*pow(s, 2) -
8280*t*pow(s, 2) - 22050*pow(s, 3) + 15705*t*pow(s, 3) + 37800*pow(s, 4) -
15552*t*pow(s, 4) - 31752*pow(s, 5) + 6804*t*pow(s, 5) + 10368*pow(s, 6) +
975*pow(t, 2) - 17010*s*pow(t, 2) + 57726*pow(s, 2)*pow(t, 2) - 51030*pow(s,
3)*pow(t, 2) + 8748*pow(s, 4)*pow(t, 2) - 2880*pow(t, 3) + 57726*s*pow(t, 3) -
171072*pow(s, 2)*pow(t, 3) + 103680*pow(s, 3)*pow(t, 3) + 405*pow(t, 4) -
51030*s*pow(t, 4) + 103680*pow(s, 2)*pow(t, 4) + 3888*pow(t, 5) + 8748*s*pow(t,
5) - 2268*pow(t, 6)))/20.; 

  d2phi(10,2) = (pow(s, 2)*(-50 + 599*s + 975*t -
11340*s*t - 2070*pow(s, 2) + 28863*t*pow(s, 2) + 3141*pow(s, 3) - 20412*t*pow(s,
3) - 2592*pow(s, 4) + 2916*t*pow(s, 4) + 972*pow(s, 5) - 4320*pow(t, 2) +
57726*s*pow(t, 2) - 128304*pow(s, 2)*pow(t, 2) + 62208*pow(s, 3)*pow(t, 2) +
810*pow(t, 3) - 68040*s*pow(t, 3) + 103680*pow(s, 2)*pow(t, 3) + 9720*pow(t, 4)
+ 14580*s*pow(t, 4) - 6804*pow(t, 5)))/20.; 

  d2phi(11,0) = (pow(t, 2)*(10 - 153*s
- 325*t + 11340*s*t - 3780*pow(s, 2) - 57726*t*pow(s, 2) + 29790*pow(s, 3) +
68040*t*pow(s, 3) - 58320*pow(s, 4) - 14580*t*pow(s, 4) + 34020*pow(s, 5) +
720*pow(t, 2) - 28863*s*pow(t, 2) + 128304*pow(s, 2)*pow(t, 2) - 103680*pow(s,
3)*pow(t, 2) - 81*pow(t, 3) + 20412*s*pow(t, 3) - 62208*pow(s, 2)*pow(t, 3) -
648*pow(t, 4) - 2916*s*pow(t, 4) + 324*pow(t, 5)))/40.; 

  d2phi(11,1) = (s*t*(20 -
153*s - 975*t + 17010*s*t - 2520*pow(s, 2) - 57726*t*pow(s, 2) + 14895*pow(s, 3)
+ 51030*t*pow(s, 3) - 23328*pow(s, 4) - 8748*t*pow(s, 4) + 11340*pow(s, 5) +
2880*pow(t, 2) - 57726*s*pow(t, 2) + 171072*pow(s, 2)*pow(t, 2) - 103680*pow(s,
3)*pow(t, 2) - 405*pow(t, 3) + 51030*s*pow(t, 3) - 103680*pow(s, 2)*pow(t, 3) -
3888*pow(t, 4) - 8748*s*pow(t, 4) + 2268*pow(t, 5)))/40.; 

  d2phi(11,2) = (pow(s,
2)*(10 - 51*s - 975*t + 11340*s*t - 630*pow(s, 2) - 28863*t*pow(s, 2) +
2979*pow(s, 3) + 20412*t*pow(s, 3) - 3888*pow(s, 4) - 2916*t*pow(s, 4) +
1620*pow(s, 5) + 4320*pow(t, 2) - 57726*s*pow(t, 2) + 128304*pow(s, 2)*pow(t, 2)
- 62208*pow(s, 3)*pow(t, 2) - 810*pow(t, 3) + 68040*s*pow(t, 3) - 103680*pow(s,
2)*pow(t, 3) - 9720*pow(t, 4) - 14580*s*pow(t, 4) + 6804*pow(t, 5)))/40.;


  d2phi(12,0) = (pow(t, 2)*(10 - 975*s - 51*t + 11340*s*t + 4320*pow(s, 2) -
57726*t*pow(s, 2) - 810*pow(s, 3) + 68040*t*pow(s, 3) - 9720*pow(s, 4) -
14580*t*pow(s, 4) + 6804*pow(s, 5) - 630*pow(t, 2) - 28863*s*pow(t, 2) +
128304*pow(s, 2)*pow(t, 2) - 103680*pow(s, 3)*pow(t, 2) + 2979*pow(t, 3) +
20412*s*pow(t, 3) - 62208*pow(s, 2)*pow(t, 3) - 3888*pow(t, 4) - 2916*s*pow(t,
4) + 1620*pow(t, 5)))/40.; 

  d2phi(12,1) = (s*t*(20 - 975*s - 153*t + 17010*s*t +
2880*pow(s, 2) - 57726*t*pow(s, 2) - 405*pow(s, 3) + 51030*t*pow(s, 3) -
3888*pow(s, 4) - 8748*t*pow(s, 4) + 2268*pow(s, 5) - 2520*pow(t, 2) -
57726*s*pow(t, 2) + 171072*pow(s, 2)*pow(t, 2) - 103680*pow(s, 3)*pow(t, 2) +
14895*pow(t, 3) + 51030*s*pow(t, 3) - 103680*pow(s, 2)*pow(t, 3) - 23328*pow(t,
4) - 8748*s*pow(t, 4) + 11340*pow(t, 5)))/40.; 

  d2phi(12,2) = (pow(s, 2)*(10 -
325*s - 153*t + 11340*s*t + 720*pow(s, 2) - 28863*t*pow(s, 2) - 81*pow(s, 3) +
20412*t*pow(s, 3) - 648*pow(s, 4) - 2916*t*pow(s, 4) + 324*pow(s, 5) -
3780*pow(t, 2) - 57726*s*pow(t, 2) + 128304*pow(s, 2)*pow(t, 2) - 62208*pow(s,
3)*pow(t, 2) + 29790*pow(t, 3) + 68040*s*pow(t, 3) - 103680*pow(s, 2)*pow(t, 3)
- 58320*pow(t, 4) - 14580*s*pow(t, 4) + 34020*pow(t, 5)))/40.; 

  d2phi(13,0) =
  (pow(t, 2)*(-50 + 975*s + 599*t - 11340*s*t - 4320*pow(s, 2) + 57726*t*pow(s,
2) + 810*pow(s, 3) - 68040*t*pow(s, 3) + 9720*pow(s, 4) + 14580*t*pow(s, 4) -
6804*pow(s, 5) - 2070*pow(t, 2) + 28863*s*pow(t, 2) - 128304*pow(s, 2)*pow(t, 2)
+ 103680*pow(s, 3)*pow(t, 2) + 3141*pow(t, 3) - 20412*s*pow(t, 3) + 62208*pow(s,
2)*pow(t, 3) - 2592*pow(t, 4) + 2916*s*pow(t, 4) + 972*pow(t, 5)))/20.;


  d2phi(13,1) = (t*(40 - 100*s - 882*t + 1797*s*t + 975*pow(s, 2) - 17010*t*pow(s,
2) - 2880*pow(s, 3) + 57726*t*pow(s, 3) + 405*pow(s, 4) - 51030*t*pow(s, 4) +
3888*pow(s, 5) + 8748*t*pow(s, 5) - 2268*pow(s, 6) + 6496*pow(t, 2) -
8280*s*pow(t, 2) + 57726*pow(s, 2)*pow(t, 2) - 171072*pow(s, 3)*pow(t, 2) +
103680*pow(s, 4)*pow(t, 2) - 22050*pow(t, 3) + 15705*s*pow(t, 3) - 51030*pow(s,
2)*pow(t, 3) + 103680*pow(s, 3)*pow(t, 3) + 37800*pow(t, 4) - 15552*s*pow(t, 4)
+ 8748*pow(s, 2)*pow(t, 4) - 31752*pow(t, 5) + 6804*s*pow(t, 5) + 10368*pow(t,
6)))/20.; 

  d2phi(13,2) = -(s*(-40 + 50*s + 1764*t - 1797*s*t - 325*pow(s, 2) +
11340*t*pow(s, 2) + 720*pow(s, 3) - 28863*t*pow(s, 3) - 81*pow(s, 4) +
20412*t*pow(s, 4) - 648*pow(s, 5) - 2916*t*pow(s, 5) + 324*pow(s, 6) -
19488*pow(t, 2) + 12420*s*pow(t, 2) - 57726*pow(s, 2)*pow(t, 2) + 128304*pow(s,
3)*pow(t, 2) - 62208*pow(s, 4)*pow(t, 2) + 88200*pow(t, 3) - 31410*s*pow(t, 3) +
68040*pow(s, 2)*pow(t, 3) - 103680*pow(s, 3)*pow(t, 3) - 189000*pow(t, 4) +
38880*s*pow(t, 4) - 14580*pow(s, 2)*pow(t, 4) + 190512*pow(t, 5) -
20412*s*pow(t, 5) - 72576*pow(t, 6)))/20.; 

  d2phi(14,0) = (pow(t, 2)*(-628 +
4017*s + 7783*t - 44712*s*t - 4644*pow(s, 2) + 32886*t*pow(s, 2) + 3870*pow(s,
3) + 25920*t*pow(s, 3) - 9720*pow(s, 4) - 14580*t*pow(s, 4) + 6804*pow(s, 5) -
33192*pow(t, 2) + 167265*s*pow(t, 2) - 122472*pow(s, 2)*pow(t, 2) - 19440*pow(s,
3)*pow(t, 2) + 61677*pow(t, 3) - 235224*s*pow(t, 3) + 116640*pow(s, 2)*pow(t, 3)
- 51516*pow(t, 4) + 107892*s*pow(t, 4) + 15876*pow(t, 5)))/40.; 

  d2phi(14,1) =
  (s*t*(-1256 + 4017*s + 23349*t - 67068*s*t - 3096*pow(s, 2) + 32886*t*pow(s,
2) + 1935*pow(s, 3) + 19440*t*pow(s, 3) - 3888*pow(s, 4) - 8748*t*pow(s, 4) +
2268*pow(s, 5) - 132768*pow(t, 2) + 334530*s*pow(t, 2) - 163296*pow(s, 2)*pow(t,
2) - 19440*pow(s, 3)*pow(t, 2) + 308385*pow(t, 3) - 588060*s*pow(t, 3) +
194400*pow(s, 2)*pow(t, 3) - 309096*pow(t, 4) + 323676*s*pow(t, 4) +
111132*pow(t, 5)))/40.; 

  d2phi(14,2) = (120*t - 628*pow(s, 2) + 23349*t*pow(s, 2)
+ 1339*pow(s, 3) - 44712*t*pow(s, 3) - 774*pow(s, 4) + 16443*t*pow(s, 4) +
387*pow(s, 5) + 7776*t*pow(s, 5) - 648*pow(s, 6) - 2916*t*pow(s, 6) + 324*pow(s,
7) - 3288*pow(t, 2) - 199152*pow(s, 2)*pow(t, 2) + 334530*pow(s, 3)*pow(t, 2) -
122472*pow(s, 4)*pow(t, 2) - 11664*pow(s, 5)*pow(t, 2) + 27400*pow(t, 3) +
616770*pow(s, 2)*pow(t, 3) - 784080*pow(s, 3)*pow(t, 3) + 194400*pow(s,
4)*pow(t, 3) - 98820*pow(t, 4) - 772740*pow(s, 2)*pow(t, 4) + 539460*pow(s,
3)*pow(t, 4) + 173124*pow(t, 5) + 333396*pow(s, 2)*pow(t, 5) - 145152*pow(t, 6)
+ 46656*pow(t, 7))/40.; 

  d2phi(15,0) = (20 - 882*s + 9864*pow(s, 2) -
46640*pow(s, 3) + 111240*pow(s, 4) - 140994*pow(s, 5) + 90720*pow(s, 6) -
23328*pow(s, 7) - 339*pow(t, 2) + 12573*s*pow(t, 2) - 105786*pow(s, 2)*pow(t, 2)
+ 324090*pow(s, 3)*pow(t, 2) - 405810*pow(s, 4)*pow(t, 2) + 176904*pow(s,
5)*pow(t, 2) + 832*pow(t, 3) - 28026*s*pow(t, 3) + 196128*pow(s, 2)*pow(t, 3) -
426060*pow(s, 3)*pow(t, 3) + 277020*pow(s, 4)*pow(t, 3) - 747*pow(t, 4) +
22653*s*pow(t, 4) - 125388*pow(s, 2)*pow(t, 4) + 149040*pow(s, 3)*pow(t, 4) +
234*pow(t, 5) - 6318*s*pow(t, 5) + 25272*pow(s, 2)*pow(t, 5))/20.; 

  d2phi(15,1) =
(3*s*t*(-226 + 4191*s + 832*t - 14013*s*t - 23508*pow(s, 2) + 65376*t*pow(s, 2)
+ 54015*pow(s, 3) - 106515*t*pow(s, 3) - 54108*pow(s, 4) + 55404*t*pow(s, 4) +
19656*pow(s, 5) - 996*pow(t, 2) + 15102*s*pow(t, 2) - 55728*pow(s, 2)*pow(t, 2)
+ 49680*pow(s, 3)*pow(t, 2) + 390*pow(t, 3) - 5265*s*pow(t, 3) + 14040*pow(s,
2)*pow(t, 3)))/20.; 

  d2phi(15,2) = (3*(-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*(-113 +
380*s + 832*t - 1854*s*t - 423*pow(s, 2) + 1026*t*pow(s, 2) + 156*pow(s, 3) -
1494*pow(t, 2) + 1656*s*pow(t, 2) + 780*pow(t, 3)))/20.; 

  d2phi(16,0) = -(t*(157
- 2877*s - 1918*t + 27153*s*t + 18102*pow(s, 2) - 128520*t*pow(s, 2) -
53550*pow(s, 3) + 270900*t*pow(s, 3) + 81270*pow(s, 4) - 262440*t*pow(s, 4) -
61236*pow(s, 5) + 95256*t*pow(s, 5) + 18144*pow(s, 6) + 9051*pow(t, 2) -
96390*s*pow(t, 2) + 325080*pow(s, 2)*pow(t, 2) - 437400*pow(s, 3)*pow(t, 2) +
204120*pow(s, 4)*pow(t, 2) - 21420*pow(t, 3) + 162540*s*pow(t, 3) -
349920*pow(s, 2)*pow(t, 3) + 226800*pow(s, 3)*pow(t, 3) + 27090*pow(t, 4) -
131220*s*pow(t, 4) + 136080*pow(s, 2)*pow(t, 4) - 17496*pow(t, 5) +
40824*s*pow(t, 5) + 4536*pow(t, 6)))/5.; 

  d2phi(16,1) = (10 - 314*s - 314*t +
7672*s*t + 2877*pow(s, 2) - 54306*t*pow(s, 2) - 12068*pow(s, 3) +
171360*t*pow(s, 3) + 26775*pow(s, 4) - 270900*t*pow(s, 4) - 32508*pow(s, 5) +
209952*t*pow(s, 5) + 20412*pow(s, 6) - 63504*t*pow(s, 6) - 5184*pow(s, 7) +
2877*pow(t, 2) - 54306*s*pow(t, 2) + 289170*pow(s, 2)*pow(t, 2) - 650160*pow(s,
3)*pow(t, 2) + 656100*pow(s, 4)*pow(t, 2) - 244944*pow(s, 5)*pow(t, 2) -
12068*pow(t, 3) + 171360*s*pow(t, 3) - 650160*pow(s, 2)*pow(t, 3) +
933120*pow(s, 3)*pow(t, 3) - 453600*pow(s, 4)*pow(t, 3) + 26775*pow(t, 4) -
270900*s*pow(t, 4) + 656100*pow(s, 2)*pow(t, 4) - 453600*pow(s, 3)*pow(t, 4) -
32508*pow(t, 5) + 209952*s*pow(t, 5) - 244944*pow(s, 2)*pow(t, 5) + 20412*pow(t,
6) - 63504*s*pow(t, 6) - 5184*pow(t, 7))/10.; 

  d2phi(16,2) = -(s*(157 - 1918*s -
2877*t + 27153*s*t + 9051*pow(s, 2) - 96390*t*pow(s, 2) - 21420*pow(s, 3) +
162540*t*pow(s, 3) + 27090*pow(s, 4) - 131220*t*pow(s, 4) - 17496*pow(s, 5) +
40824*t*pow(s, 5) + 4536*pow(s, 6) + 18102*pow(t, 2) - 128520*s*pow(t, 2) +
325080*pow(s, 2)*pow(t, 2) - 349920*pow(s, 3)*pow(t, 2) + 136080*pow(s,
4)*pow(t, 2) - 53550*pow(t, 3) + 270900*s*pow(t, 3) - 437400*pow(s, 2)*pow(t, 3)
+ 226800*pow(s, 3)*pow(t, 3) + 81270*pow(t, 4) - 262440*s*pow(t, 4) +
204120*pow(s, 2)*pow(t, 4) - 61236*pow(t, 5) + 95256*s*pow(t, 5) + 18144*pow(t,
6)))/5.; 

  d2phi(17,0) = (3*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(-113 + 832*s + 380*t
- 1854*s*t - 1494*pow(s, 2) + 1656*t*pow(s, 2) + 780*pow(s, 3) - 423*pow(t, 2) +
1026*s*pow(t, 2) + 156*pow(t, 3)))/20.; 

  d2phi(17,1) = (3*s*t*(-226 + 832*s +
4191*t - 14013*s*t - 996*pow(s, 2) + 15102*t*pow(s, 2) + 390*pow(s, 3) -
5265*t*pow(s, 3) - 23508*pow(t, 2) + 65376*s*pow(t, 2) - 55728*pow(s, 2)*pow(t,
2) + 14040*pow(s, 3)*pow(t, 2) + 54015*pow(t, 3) - 106515*s*pow(t, 3) +
49680*pow(s, 2)*pow(t, 3) - 54108*pow(t, 4) + 55404*s*pow(t, 4) + 19656*pow(t,
5)))/20.; 

  d2phi(17,2) = (20 - 882*t - 339*pow(s, 2) + 12573*t*pow(s, 2) +
832*pow(s, 3) - 28026*t*pow(s, 3) - 747*pow(s, 4) + 22653*t*pow(s, 4) +
234*pow(s, 5) - 6318*t*pow(s, 5) + 9864*pow(t, 2) - 105786*pow(s, 2)*pow(t, 2) +
196128*pow(s, 3)*pow(t, 2) - 125388*pow(s, 4)*pow(t, 2) + 25272*pow(s, 5)*pow(t,
2) - 46640*pow(t, 3) + 324090*pow(s, 2)*pow(t, 3) - 426060*pow(s, 3)*pow(t, 3) +
149040*pow(s, 4)*pow(t, 3) + 111240*pow(t, 4) - 405810*pow(s, 2)*pow(t, 4) +
277020*pow(s, 3)*pow(t, 4) - 140994*pow(t, 5) + 176904*pow(s, 2)*pow(t, 5) +
90720*pow(t, 6) - 23328*pow(t, 7))/20.; 

  d2phi(18,0) = -32*(-1 + 3*t)*(-1 +
6*t)*pow(t, 2)*(-47 + 246*s + 164*t - 567*s*t - 378*pow(s, 2) + 432*t*pow(s, 2)
+ 180*pow(s, 3) - 189*pow(t, 2) + 324*s*pow(t, 2) + 72*pow(t, 3)); 

  d2phi(18,1) =
-16*t*(20 - 188*s - 411*t + 3522*s*t + 492*pow(s, 2) - 8343*t*pow(s, 2) -
504*pow(s, 3) + 7668*t*pow(s, 3) + 180*pow(s, 4) - 2430*t*pow(s, 4) +
2740*pow(t, 2) - 20088*s*pow(t, 2) + 39420*pow(s, 2)*pow(t, 2) - 28512*pow(s,
3)*pow(t, 2) + 6480*pow(s, 4)*pow(t, 2) - 8235*pow(t, 3) + 47250*s*pow(t, 3) -
65610*pow(s, 2)*pow(t, 3) + 25920*pow(s, 3)*pow(t, 3) + 12366*pow(t, 4) -
48600*s*pow(t, 4) + 34992*pow(s, 2)*pow(t, 4) - 9072*pow(t, 5) + 18144*s*pow(t,
5) + 2592*pow(t, 6)); 

  d2phi(18,2) = -32*s*(10 - 47*s - 411*t + 1761*s*t +
82*pow(s, 2) - 2781*t*pow(s, 2) - 63*pow(s, 3) + 1917*t*pow(s, 3) + 18*pow(s, 4)
- 486*t*pow(s, 4) + 4110*pow(t, 2) - 15066*s*pow(t, 2) + 19710*pow(s, 2)*pow(t,
2) - 10692*pow(s, 3)*pow(t, 2) + 1944*pow(s, 4)*pow(t, 2) - 16470*pow(t, 3) +
47250*s*pow(t, 3) - 43740*pow(s, 2)*pow(t, 3) + 12960*pow(s, 3)*pow(t, 3) +
30915*pow(t, 4) - 60750*s*pow(t, 4) + 29160*pow(s, 2)*pow(t, 4) - 27216*pow(t,
5) + 27216*s*pow(t, 5) + 9072*pow(t, 6)); 

  d2phi(19,0) = -32*t*(10 - 411*s - 47*t
+ 1761*s*t + 4110*pow(s, 2) - 15066*t*pow(s, 2) - 16470*pow(s, 3) +
47250*t*pow(s, 3) + 30915*pow(s, 4) - 60750*t*pow(s, 4) - 27216*pow(s, 5) +
27216*t*pow(s, 5) + 9072*pow(s, 6) + 82*pow(t, 2) - 2781*s*pow(t, 2) +
19710*pow(s, 2)*pow(t, 2) - 43740*pow(s, 3)*pow(t, 2) + 29160*pow(s, 4)*pow(t,
2) - 63*pow(t, 3) + 1917*s*pow(t, 3) - 10692*pow(s, 2)*pow(t, 3) + 12960*pow(s,
3)*pow(t, 3) + 18*pow(t, 4) - 486*s*pow(t, 4) + 1944*pow(s, 2)*pow(t, 4));


  d2phi(19,1) = -16*s*(20 - 411*s - 188*t + 3522*s*t + 2740*pow(s, 2) -
20088*t*pow(s, 2) - 8235*pow(s, 3) + 47250*t*pow(s, 3) + 12366*pow(s, 4) -
48600*t*pow(s, 4) - 9072*pow(s, 5) + 18144*t*pow(s, 5) + 2592*pow(s, 6) +
492*pow(t, 2) - 8343*s*pow(t, 2) + 39420*pow(s, 2)*pow(t, 2) - 65610*pow(s,
3)*pow(t, 2) + 34992*pow(s, 4)*pow(t, 2) - 504*pow(t, 3) + 7668*s*pow(t, 3) -
28512*pow(s, 2)*pow(t, 3) + 25920*pow(s, 3)*pow(t, 3) + 180*pow(t, 4) -
2430*s*pow(t, 4) + 6480*pow(s, 2)*pow(t, 4)); 

  d2phi(19,2) = -32*(-1 + 3*s)*(-1 +
6*s)*pow(s, 2)*(-47 + 164*s + 246*t - 567*s*t - 189*pow(s, 2) + 324*t*pow(s, 2)
+ 72*pow(s, 3) - 378*pow(t, 2) + 432*s*pow(t, 2) + 180*pow(t, 3)); 

  d2phi(20,0) =
16*(-1 + 3*t)*(-1 + 6*t)*(-1 + 30*s + t - 27*s*t - 162*pow(s, 2) + 108*t*pow(s,
2) + 180*pow(s, 3))*pow(t, 2)*sqrt(2); 

  d2phi(20,1) = 8*s*t*(-4 + 60*s + 60*t -
891*s*t - 216*pow(s, 2) + 3132*t*pow(s, 2) + 180*pow(s, 3) - 2430*t*pow(s, 3) -
216*pow(t, 2) + 3132*s*pow(t, 2) - 10368*pow(s, 2)*pow(t, 2) + 6480*pow(s,
3)*pow(t, 2) + 180*pow(t, 3) - 2430*s*pow(t, 3) + 6480*pow(s, 2)*pow(t,
3))*sqrt(2); 

  d2phi(20,2) = 16*(-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*(-1 + s + 30*t -
27*s*t - 162*pow(t, 2) + 108*s*pow(t, 2) + 180*pow(t, 3))*sqrt(2); 

  d2phi(21,0) =
(-11664*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(-93 + 792*s + 290*t - 1714*s*t -
1494*pow(s, 2) + 1656*t*pow(s, 2) + 780*pow(s, 3) - 293*pow(t, 2) + 906*s*pow(t,
2) + 96*pow(t, 3)))/125.; 

  d2phi(21,1) = (-11664*s*t*(-186 + 792*s + 3381*t -
13263*s*t - 996*pow(s, 2) + 15102*t*pow(s, 2) + 390*pow(s, 3) - 5265*t*pow(s, 3)
- 18308*pow(t, 2) + 61176*s*pow(t, 2) - 55728*pow(s, 2)*pow(t, 2) + 14040*pow(s,
3)*pow(t, 2) + 39765*pow(t, 3) - 97515*s*pow(t, 3) + 49680*pow(s, 2)*pow(t, 3) -
36828*pow(t, 4) + 48924*s*pow(t, 4) + 12096*pow(t, 5)))/125.; 

  d2phi(21,2) =
(-11664*(16*t - 93*pow(s, 2) + 3381*t*pow(s, 2) + 264*pow(s, 3) - 8842*t*pow(s,
3) - 249*pow(s, 4) + 7551*t*pow(s, 4) + 78*pow(s, 5) - 2106*t*pow(s, 5) -
432*pow(t, 2) - 27462*pow(s, 2)*pow(t, 2) + 61176*pow(s, 3)*pow(t, 2) -
41796*pow(s, 4)*pow(t, 2) + 8424*pow(s, 5)*pow(t, 2) + 3520*pow(t, 3) +
79530*pow(s, 2)*pow(t, 3) - 130020*pow(s, 3)*pow(t, 3) + 49680*pow(s, 4)*pow(t,
3) - 12320*pow(t, 4) - 92070*pow(s, 2)*pow(t, 4) + 81540*pow(s, 3)*pow(t, 4) +
20832*pow(t, 5) + 36288*pow(s, 2)*pow(t, 5) - 16800*pow(t, 6) + 5184*pow(t,
7)))/125.; 

  d2phi(22,0) = (729*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(-61 + 600*s +
168*t - 1150*s*t - 1302*pow(s, 2) + 1272*t*pow(s, 2) + 780*pow(s, 3) -
151*pow(t, 2) + 546*s*pow(t, 2) + 44*pow(t, 3)))/8.; 

  d2phi(22,1) =
(729*s*t*(-122 + 600*s + 2151*t - 9825*s*t - 868*pow(s, 2) + 12990*t*pow(s, 2) +
390*pow(s, 3) - 5265*t*pow(s, 3) - 11044*pow(t, 2) + 43392*s*pow(t, 2) -
46512*pow(s, 2)*pow(t, 2) + 14040*pow(s, 3)*pow(t, 2) + 22135*pow(t, 3) -
64035*s*pow(t, 3) + 38160*pow(s, 2)*pow(t, 3) - 18684*pow(t, 4) + 29484*s*pow(t,
4) + 5544*pow(t, 5)))/8.; 

  d2phi(22,2) = (729*(10*t - 61*pow(s, 2) +
2151*t*pow(s, 2) + 200*pow(s, 3) - 6550*t*pow(s, 3) - 217*pow(s, 4) +
6495*t*pow(s, 4) + 78*pow(s, 5) - 2106*t*pow(s, 5) - 264*pow(t, 2) -
16566*pow(s, 2)*pow(t, 2) + 43392*pow(s, 3)*pow(t, 2) - 34884*pow(s, 4)*pow(t,
2) + 8424*pow(s, 5)*pow(t, 2) + 2080*pow(t, 3) + 44270*pow(s, 2)*pow(t, 3) -
85380*pow(s, 3)*pow(t, 3) + 38160*pow(s, 4)*pow(t, 3) - 6980*pow(t, 4) -
46710*pow(s, 2)*pow(t, 4) + 49140*pow(s, 3)*pow(t, 4) + 11298*pow(t, 5) +
16632*pow(s, 2)*pow(t, 5) - 8736*pow(t, 6) + 2592*pow(t, 7)))/8.; 

  d2phi(23,0) =
(-81*(-1 + 6*t)*pow(t, 2)*(1003 - 14712*s - 3635*t + 43134*s*t + 52890*pow(s, 2)
- 106758*t*pow(s, 2) - 69300*pow(s, 3) + 71100*t*pow(s, 3) + 30240*pow(s, 4) +
4603*pow(t, 2) - 40932*s*pow(t, 2) + 52920*pow(s, 2)*pow(t, 2) - 2295*pow(t, 3)
+ 12474*s*pow(t, 3) + 324*pow(t, 4)))/8.; 

  d2phi(23,1) = (-81*s*t*(-2006 +
14712*s + 28959*t - 197109*s*t - 35260*pow(s, 2) + 424098*t*pow(s, 2) +
34650*pow(s, 3) - 365175*t*pow(s, 3) - 12096*pow(s, 4) + 108864*t*pow(s, 4) -
105652*pow(t, 2) + 599472*s*pow(t, 2) - 924624*pow(s, 2)*pow(t, 2) +
426600*pow(s, 3)*pow(t, 2) + 149565*pow(t, 3) - 645165*s*pow(t, 3) +
529200*pow(s, 2)*pow(t, 3) - 84564*pow(t, 4) + 224532*s*pow(t, 4) + 13608*pow(t,
5)))/8.; 

  d2phi(23,2) = (-81*(180*t - 1003*pow(s, 2) + 28959*t*pow(s, 2) +
4904*pow(s, 3) - 131406*t*pow(s, 3) - 8815*pow(s, 4) + 212049*t*pow(s, 4) +
6930*pow(s, 5) - 146070*t*pow(s, 5) - 2016*pow(s, 6) + 36288*t*pow(s, 6) -
4212*pow(t, 2) - 158478*pow(s, 2)*pow(t, 2) + 599472*pow(s, 3)*pow(t, 2) -
693468*pow(s, 4)*pow(t, 2) + 255960*pow(s, 5)*pow(t, 2) + 28260*pow(t, 3) +
299130*pow(s, 2)*pow(t, 3) - 860220*pow(s, 3)*pow(t, 3) + 529200*pow(s,
4)*pow(t, 3) - 82710*pow(t, 4) - 211410*pow(s, 2)*pow(t, 4) + 374220*pow(s,
3)*pow(t, 4) + 119826*pow(t, 5) + 40824*pow(s, 2)*pow(t, 5) - 84672*pow(t, 6) +
23328*pow(t, 7)))/8.; 

  d2phi(24,0) = (-648*pow(t, 2)*(5929 - 111993*s - 32876*t +
553884*s*t + 553578*pow(s, 2) - 2146824*t*pow(s, 2) - 1123410*pow(s, 3) +
2966040*t*pow(s, 3) + 1013580*pow(s, 4) - 1351080*t*pow(s, 4) - 337932*pow(s, 5)
+ 58369*pow(t, 2) - 910791*s*pow(t, 2) + 2527632*pow(s, 2)*pow(t, 2) -
1804680*pow(s, 3)*pow(t, 2) - 31746*pow(t, 3) + 584172*s*pow(t, 3) -
914976*pow(s, 2)*pow(t, 3) - 10044*pow(t, 4) - 114696*s*pow(t, 4) + 10368*pow(t,
5)))/125.; 

  d2phi(24,1) = (648*s*t*(-11858 + 111993*s + 98628*t - 830826*s*t -
369052*pow(s, 2) + 2146824*t*pow(s, 2) + 561705*pow(s, 3) - 2224530*t*pow(s, 3)
- 405432*pow(s, 4) + 810648*t*pow(s, 4) + 112644*pow(s, 5) - 233476*pow(t, 2) +
1821582*s*pow(t, 2) - 3370176*pow(s, 2)*pow(t, 2) + 1804680*pow(s, 3)*pow(t, 2)
+ 158730*pow(t, 3) - 1460430*s*pow(t, 3) + 1524960*pow(s, 2)*pow(t, 3) +
60264*pow(t, 4) + 344088*s*pow(t, 4) - 72576*pow(t, 5)))/125.; 

  d2phi(24,2) =
(648*(1440*t - 5929*pow(s, 2) + 98628*t*pow(s, 2) + 37331*pow(s, 3) -
553884*t*pow(s, 3) - 92263*pow(s, 4) + 1073412*t*pow(s, 4) + 112341*pow(s, 5) -
889812*t*pow(s, 5) - 67572*pow(s, 6) + 270216*t*pow(s, 6) + 16092*pow(s, 7) -
25056*pow(t, 2) - 350214*pow(s, 2)*pow(t, 2) + 1821582*pow(s, 3)*pow(t, 2) -
2527632*pow(s, 4)*pow(t, 2) + 1082808*pow(s, 5)*pow(t, 2) + 144000*pow(t, 3) +
317460*pow(s, 2)*pow(t, 3) - 1947240*pow(s, 3)*pow(t, 3) + 1524960*pow(s,
4)*pow(t, 3) - 383040*pow(t, 4) + 150660*pow(s, 2)*pow(t, 4) + 573480*pow(s,
3)*pow(t, 4) + 520128*pow(t, 5) - 217728*pow(s, 2)*pow(t, 5) - 350784*pow(t, 6)
+ 93312*pow(t, 7)))/125.; 

  d2phi(25,0) = (648*(1440*s - 25056*pow(s, 2) +
144000*pow(s, 3) - 383040*pow(s, 4) + 520128*pow(s, 5) - 350784*pow(s, 6) +
93312*pow(s, 7) - 5929*pow(t, 2) + 98628*s*pow(t, 2) - 350214*pow(s, 2)*pow(t,
2) + 317460*pow(s, 3)*pow(t, 2) + 150660*pow(s, 4)*pow(t, 2) - 217728*pow(s,
5)*pow(t, 2) + 37331*pow(t, 3) - 553884*s*pow(t, 3) + 1821582*pow(s, 2)*pow(t,
3) - 1947240*pow(s, 3)*pow(t, 3) + 573480*pow(s, 4)*pow(t, 3) - 92263*pow(t, 4)
+ 1073412*s*pow(t, 4) - 2527632*pow(s, 2)*pow(t, 4) + 1524960*pow(s, 3)*pow(t,
4) + 112341*pow(t, 5) - 889812*s*pow(t, 5) + 1082808*pow(s, 2)*pow(t, 5) -
67572*pow(t, 6) + 270216*s*pow(t, 6) + 16092*pow(t, 7)))/125.; 

  d2phi(25,1) =
(-648*s*t*(11858 - 98628*s - 111993*t + 830826*s*t + 233476*pow(s, 2) -
1821582*t*pow(s, 2) - 158730*pow(s, 3) + 1460430*t*pow(s, 3) - 60264*pow(s, 4) -
344088*t*pow(s, 4) + 72576*pow(s, 5) + 369052*pow(t, 2) - 2146824*s*pow(t, 2) +
3370176*pow(s, 2)*pow(t, 2) - 1524960*pow(s, 3)*pow(t, 2) - 561705*pow(t, 3) +
2224530*s*pow(t, 3) - 1804680*pow(s, 2)*pow(t, 3) + 405432*pow(t, 4) -
810648*s*pow(t, 4) - 112644*pow(t, 5)))/125.; 

  d2phi(25,2) = (-648*pow(s,
2)*(5929 - 32876*s - 111993*t + 553884*s*t + 58369*pow(s, 2) - 910791*t*pow(s,
2) - 31746*pow(s, 3) + 584172*t*pow(s, 3) - 10044*pow(s, 4) - 114696*t*pow(s, 4)
+ 10368*pow(s, 5) + 553578*pow(t, 2) - 2146824*s*pow(t, 2) + 2527632*pow(s,
2)*pow(t, 2) - 914976*pow(s, 3)*pow(t, 2) - 1123410*pow(t, 3) + 2966040*s*pow(t,
3) - 1804680*pow(s, 2)*pow(t, 3) + 1013580*pow(t, 4) - 1351080*s*pow(t, 4) -
337932*pow(t, 5)))/125.; 

  d2phi(26,0) = (-81*(180*s - 4212*pow(s, 2) +
28260*pow(s, 3) - 82710*pow(s, 4) + 119826*pow(s, 5) - 84672*pow(s, 6) +
23328*pow(s, 7) - 1003*pow(t, 2) + 28959*s*pow(t, 2) - 158478*pow(s, 2)*pow(t,
2) + 299130*pow(s, 3)*pow(t, 2) - 211410*pow(s, 4)*pow(t, 2) + 40824*pow(s,
5)*pow(t, 2) + 4904*pow(t, 3) - 131406*s*pow(t, 3) + 599472*pow(s, 2)*pow(t, 3)
- 860220*pow(s, 3)*pow(t, 3) + 374220*pow(s, 4)*pow(t, 3) - 8815*pow(t, 4) +
212049*s*pow(t, 4) - 693468*pow(s, 2)*pow(t, 4) + 529200*pow(s, 3)*pow(t, 4) +
6930*pow(t, 5) - 146070*s*pow(t, 5) + 255960*pow(s, 2)*pow(t, 5) - 2016*pow(t,
6) + 36288*s*pow(t, 6)))/8.; 

  d2phi(26,1) = (-81*s*t*(-2006 + 28959*s + 14712*t -
197109*s*t - 105652*pow(s, 2) + 599472*t*pow(s, 2) + 149565*pow(s, 3) -
645165*t*pow(s, 3) - 84564*pow(s, 4) + 224532*t*pow(s, 4) + 13608*pow(s, 5) -
35260*pow(t, 2) + 424098*s*pow(t, 2) - 924624*pow(s, 2)*pow(t, 2) +
529200*pow(s, 3)*pow(t, 2) + 34650*pow(t, 3) - 365175*s*pow(t, 3) +
426600*pow(s, 2)*pow(t, 3) - 12096*pow(t, 4) + 108864*s*pow(t, 4)))/8.;


  d2phi(26,2) = (-81*(-1 + 6*s)*pow(s, 2)*(1003 - 3635*s - 14712*t + 43134*s*t +
4603*pow(s, 2) - 40932*t*pow(s, 2) - 2295*pow(s, 3) + 12474*t*pow(s, 3) +
324*pow(s, 4) + 52890*pow(t, 2) - 106758*s*pow(t, 2) + 52920*pow(s, 2)*pow(t, 2)
- 69300*pow(t, 3) + 71100*s*pow(t, 3) + 30240*pow(t, 4)))/8.; 

  d2phi(27,0) =
  (729*(10*s - 264*pow(s, 2) + 2080*pow(s, 3) - 6980*pow(s, 4) + 11298*pow(s, 5)
- 8736*pow(s, 6) + 2592*pow(s, 7) - 61*pow(t, 2) + 2151*s*pow(t, 2) -
16566*pow(s, 2)*pow(t, 2) + 44270*pow(s, 3)*pow(t, 2) - 46710*pow(s, 4)*pow(t,
2) + 16632*pow(s, 5)*pow(t, 2) + 200*pow(t, 3) - 6550*s*pow(t, 3) + 43392*pow(s,
2)*pow(t, 3) - 85380*pow(s, 3)*pow(t, 3) + 49140*pow(s, 4)*pow(t, 3) -
217*pow(t, 4) + 6495*s*pow(t, 4) - 34884*pow(s, 2)*pow(t, 4) + 38160*pow(s,
3)*pow(t, 4) + 78*pow(t, 5) - 2106*s*pow(t, 5) + 8424*pow(s, 2)*pow(t, 5)))/8.;


  d2phi(27,1) = (729*s*t*(-122 + 2151*s + 600*t - 9825*s*t - 11044*pow(s, 2) +
43392*t*pow(s, 2) + 22135*pow(s, 3) - 64035*t*pow(s, 3) - 18684*pow(s, 4) +
29484*t*pow(s, 4) + 5544*pow(s, 5) - 868*pow(t, 2) + 12990*s*pow(t, 2) -
46512*pow(s, 2)*pow(t, 2) + 38160*pow(s, 3)*pow(t, 2) + 390*pow(t, 3) -
5265*s*pow(t, 3) + 14040*pow(s, 2)*pow(t, 3)))/8.; 

  d2phi(27,2) = (729*(-1 +
3*s)*(-1 + 6*s)*pow(s, 2)*(-61 + 168*s + 600*t - 1150*s*t - 151*pow(s, 2) +
546*t*pow(s, 2) + 44*pow(s, 3) - 1302*pow(t, 2) + 1272*s*pow(t, 2) + 780*pow(t,
3)))/8.; 

  d2phi(28,0) = (-11664*(16*s - 432*pow(s, 2) + 3520*pow(s, 3) -
12320*pow(s, 4) + 20832*pow(s, 5) - 16800*pow(s, 6) + 5184*pow(s, 7) - 93*pow(t,
2) + 3381*s*pow(t, 2) - 27462*pow(s, 2)*pow(t, 2) + 79530*pow(s, 3)*pow(t, 2) -
92070*pow(s, 4)*pow(t, 2) + 36288*pow(s, 5)*pow(t, 2) + 264*pow(t, 3) -
8842*s*pow(t, 3) + 61176*pow(s, 2)*pow(t, 3) - 130020*pow(s, 3)*pow(t, 3) +
81540*pow(s, 4)*pow(t, 3) - 249*pow(t, 4) + 7551*s*pow(t, 4) - 41796*pow(s,
2)*pow(t, 4) + 49680*pow(s, 3)*pow(t, 4) + 78*pow(t, 5) - 2106*s*pow(t, 5) +
8424*pow(s, 2)*pow(t, 5)))/125.; 

  d2phi(28,1) = (-11664*s*t*(-186 + 3381*s +
792*t - 13263*s*t - 18308*pow(s, 2) + 61176*t*pow(s, 2) + 39765*pow(s, 3) -
97515*t*pow(s, 3) - 36828*pow(s, 4) + 48924*t*pow(s, 4) + 12096*pow(s, 5) -
996*pow(t, 2) + 15102*s*pow(t, 2) - 55728*pow(s, 2)*pow(t, 2) + 49680*pow(s,
3)*pow(t, 2) + 390*pow(t, 3) - 5265*s*pow(t, 3) + 14040*pow(s, 2)*pow(t,
3)))/125.; 

  d2phi(28,2) = (-11664*(-1 + 3*s)*(-1 + 6*s)*pow(s, 2)*(-93 + 290*s +
792*t - 1714*s*t - 293*pow(s, 2) + 906*t*pow(s, 2) + 96*pow(s, 3) - 1494*pow(t,
2) + 1656*s*pow(t, 2) + 780*pow(t, 3)))/125.; 

  d2phi(29,0) = (-324*pow(t,
2)*(-518 + 22449*s - 1313*t + 42000*s*t - 244200*pow(s, 2) - 131346*t*pow(s, 2)
+ 970590*pow(s, 3) - 200160*t*pow(s, 3) - 1512000*pow(s, 4) + 445500*t*pow(s, 4)
+ 786996*pow(s, 5) + 4540*pow(t, 2) - 170253*s*pow(t, 2) + 756864*pow(s,
2)*pow(t, 2) - 589680*pow(s, 3)*pow(t, 2) - 1161*pow(t, 3) + 121392*s*pow(t, 3)
- 379728*pow(s, 2)*pow(t, 3) - 3168*pow(t, 4) - 15876*s*pow(t, 4) + 1620*pow(t,
5)))/125.; 

  d2phi(29,1) = (-324*s*t*(-1036 + 22449*s - 3939*t + 63000*s*t -
162800*pow(s, 2) - 131346*t*pow(s, 2) + 485295*pow(s, 3) - 150120*t*pow(s, 3) -
604800*pow(s, 4) + 267300*t*pow(s, 4) + 262332*pow(s, 5) + 18160*pow(t, 2) -
340506*s*pow(t, 2) + 1009152*pow(s, 2)*pow(t, 2) - 589680*pow(s, 3)*pow(t, 2) -
5805*pow(t, 3) + 303480*s*pow(t, 3) - 632880*pow(s, 2)*pow(t, 3) - 19008*pow(t,
4) - 47628*s*pow(t, 4) + 11340*pow(t, 5)))/125.; 

  d2phi(29,2) = (-324*pow(s,
2)*(-518 + 7483*s - 3939*t + 42000*s*t - 40700*pow(s, 2) - 65673*t*pow(s, 2) +
97059*pow(s, 3) - 60048*t*pow(s, 3) - 100800*pow(s, 4) + 89100*t*pow(s, 4) +
37476*pow(s, 5) + 27240*pow(t, 2) - 340506*s*pow(t, 2) + 756864*pow(s, 2)*pow(t,
2) - 353808*pow(s, 3)*pow(t, 2) - 11610*pow(t, 3) + 404640*s*pow(t, 3) -
632880*pow(s, 2)*pow(t, 3) - 47520*pow(t, 4) - 79380*s*pow(t, 4) + 34020*pow(t,
5)))/125.; 

  d2phi(30,0) = (81*pow(t, 2)*(274 - 9501*s - 2021*t + 68064*s*t +
67416*pow(s, 2) - 432522*t*pow(s, 2) - 158670*pow(s, 3) + 820440*t*pow(s, 3) +
137160*pow(s, 4) - 458460*t*pow(s, 4) - 34020*pow(s, 5) + 3382*pow(t, 2) -
113301*s*pow(t, 2) + 566352*pow(s, 2)*pow(t, 2) - 587520*pow(s, 3)*pow(t, 2) -
1167*pow(t, 3) + 60804*s*pow(t, 3) - 202176*pow(s, 2)*pow(t, 3) - 1008*pow(t, 4)
- 6156*s*pow(t, 4) + 540*pow(t, 5)))/16.; 

  d2phi(30,1) = (-81*s*t*(-548 + 9501*s
  + 6063*t - 102096*s*t - 44944*pow(s, 2) + 432522*t*pow(s, 2) + 79335*pow(s, 3)
- 615330*t*pow(s, 3) - 54864*pow(s, 4) + 275076*t*pow(s, 4) + 11340*pow(s, 5) -
13528*pow(t, 2) + 226602*s*pow(t, 2) - 755136*pow(s, 2)*pow(t, 2) +
587520*pow(s, 3)*pow(t, 2) + 5835*pow(t, 3) - 152010*s*pow(t, 3) + 336960*pow(s,
2)*pow(t, 3) + 6048*pow(t, 4) + 18468*s*pow(t, 4) - 3780*pow(t, 5)))/16.;


  d2phi(30,2) = (-81*pow(s, 2)*(-274 + 3167*s + 6063*t - 68064*s*t - 11236*pow(s,
2) + 216261*t*pow(s, 2) + 15867*pow(s, 3) - 246132*t*pow(s, 3) - 9144*pow(s, 4)
+ 91692*t*pow(s, 4) + 1620*pow(s, 5) - 20292*pow(t, 2) + 226602*s*pow(t, 2) -
566352*pow(s, 2)*pow(t, 2) + 352512*pow(s, 3)*pow(t, 2) + 11670*pow(t, 3) -
202680*s*pow(t, 3) + 336960*pow(s, 2)*pow(t, 3) + 15120*pow(t, 4) +
30780*s*pow(t, 4) - 11340*pow(t, 5)))/16.; 

  d2phi(31,0) = (-81*pow(t, 2)*(-274 +
6063*s + 3167*t - 68064*s*t - 20292*pow(s, 2) + 226602*t*pow(s, 2) +
11670*pow(s, 3) - 202680*t*pow(s, 3) + 15120*pow(s, 4) + 30780*t*pow(s, 4) -
11340*pow(s, 5) - 11236*pow(t, 2) + 216261*s*pow(t, 2) - 566352*pow(s, 2)*pow(t,
2) + 336960*pow(s, 3)*pow(t, 2) + 15867*pow(t, 3) - 246132*s*pow(t, 3) +
352512*pow(s, 2)*pow(t, 3) - 9144*pow(t, 4) + 91692*s*pow(t, 4) + 1620*pow(t,
5)))/16.; 

  d2phi(31,1) = (81*s*t*(548 - 6063*s - 9501*t + 102096*s*t +
13528*pow(s, 2) - 226602*t*pow(s, 2) - 5835*pow(s, 3) + 152010*t*pow(s, 3) -
6048*pow(s, 4) - 18468*t*pow(s, 4) + 3780*pow(s, 5) + 44944*pow(t, 2) -
432522*s*pow(t, 2) + 755136*pow(s, 2)*pow(t, 2) - 336960*pow(s, 3)*pow(t, 2) -
79335*pow(t, 3) + 615330*s*pow(t, 3) - 587520*pow(s, 2)*pow(t, 3) + 54864*pow(t,
4) - 275076*s*pow(t, 4) - 11340*pow(t, 5)))/16.; 

  d2phi(31,2) = (81*pow(s,
2)*(274 - 2021*s - 9501*t + 68064*s*t + 3382*pow(s, 2) - 113301*t*pow(s, 2) -
1167*pow(s, 3) + 60804*t*pow(s, 3) - 1008*pow(s, 4) - 6156*t*pow(s, 4) +
540*pow(s, 5) + 67416*pow(t, 2) - 432522*s*pow(t, 2) + 566352*pow(s, 2)*pow(t,
2) - 202176*pow(s, 3)*pow(t, 2) - 158670*pow(t, 3) + 820440*s*pow(t, 3) -
587520*pow(s, 2)*pow(t, 3) + 137160*pow(t, 4) - 458460*s*pow(t, 4) -
34020*pow(t, 5)))/16.; 

  d2phi(32,0) = (-324*pow(t, 2)*(-518 - 3939*s + 7483*t +
42000*s*t + 27240*pow(s, 2) - 340506*t*pow(s, 2) - 11610*pow(s, 3) +
404640*t*pow(s, 3) - 47520*pow(s, 4) - 79380*t*pow(s, 4) + 34020*pow(s, 5) -
40700*pow(t, 2) - 65673*s*pow(t, 2) + 756864*pow(s, 2)*pow(t, 2) - 632880*pow(s,
3)*pow(t, 2) + 97059*pow(t, 3) - 60048*s*pow(t, 3) - 353808*pow(s, 2)*pow(t, 3)
- 100800*pow(t, 4) + 89100*s*pow(t, 4) + 37476*pow(t, 5)))/125.; 

  d2phi(32,1) =
  (-324*s*t*(-1036 - 3939*s + 22449*t + 63000*s*t + 18160*pow(s, 2) -
340506*t*pow(s, 2) - 5805*pow(s, 3) + 303480*t*pow(s, 3) - 19008*pow(s, 4) -
47628*t*pow(s, 4) + 11340*pow(s, 5) - 162800*pow(t, 2) - 131346*s*pow(t, 2) +
1009152*pow(s, 2)*pow(t, 2) - 632880*pow(s, 3)*pow(t, 2) + 485295*pow(t, 3) -
150120*s*pow(t, 3) - 589680*pow(s, 2)*pow(t, 3) - 604800*pow(t, 4) +
267300*s*pow(t, 4) + 262332*pow(t, 5)))/125.; 

  d2phi(32,2) = (-324*pow(s,
2)*(-518 - 1313*s + 22449*t + 42000*s*t + 4540*pow(s, 2) - 170253*t*pow(s, 2) -
1161*pow(s, 3) + 121392*t*pow(s, 3) - 3168*pow(s, 4) - 15876*t*pow(s, 4) +
1620*pow(s, 5) - 244200*pow(t, 2) - 131346*s*pow(t, 2) + 756864*pow(s, 2)*pow(t,
2) - 379728*pow(s, 3)*pow(t, 2) + 970590*pow(t, 3) - 200160*s*pow(t, 3) -
589680*pow(s, 2)*pow(t, 3) - 1512000*pow(t, 4) + 445500*s*pow(t, 4) +
786996*pow(t, 5)))/125.; 

  d2phi(33,0) = (-1296*(-1 + 2*t)*(-2 + 3*s + 2*t)*(-2 +
3*t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2))/25.; 

  d2phi(33,1) = (-648*t*(4 - 16*s -
87*t + 324*s*t + 12*pow(s, 2) - 225*t*pow(s, 2) + 628*pow(t, 2) - 2080*s*pow(t,
2) + 1260*pow(s, 2)*pow(t, 2) - 2075*pow(t, 3) + 5700*s*pow(t, 3) - 2700*pow(s,
2)*pow(t, 3) + 3438*pow(t, 4) - 6912*s*pow(t, 4) + 1944*pow(s, 2)*pow(t, 4) -
2772*pow(t, 5) + 3024*s*pow(t, 5) + 864*pow(t, 6)))/25.; 

  d2phi(33,2) =
(-1296*s*(2 - 4*s - 87*t + 162*s*t + 2*pow(s, 2) - 75*t*pow(s, 2) + 942*pow(t,
2) - 1560*s*pow(t, 2) + 630*pow(s, 2)*pow(t, 2) - 4150*pow(t, 3) + 5700*s*pow(t,
3) - 1800*pow(s, 2)*pow(t, 3) + 8595*pow(t, 4) - 8640*s*pow(t, 4) + 1620*pow(s,
2)*pow(t, 4) - 8316*pow(t, 5) + 4536*s*pow(t, 5) + 3024*pow(t, 6)))/25.;


  d2phi(34,0) = (81*(-1 + 2*t)*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(16 - 51*s - 34*t +
54*s*t + 36*pow(s, 2) + 18*pow(t, 2)))/2.; 

  d2phi(34,1) = (81*t*(10 - 64*s -
213*t + 1260*s*t + 102*pow(s, 2) - 1845*t*pow(s, 2) - 48*pow(s, 3) +
792*t*pow(s, 3) + 1492*pow(t, 2) - 7744*s*pow(t, 2) + 9720*pow(s, 2)*pow(t, 2) -
3456*pow(s, 3)*pow(t, 2) - 4745*pow(t, 3) + 19980*s*pow(t, 3) - 18900*pow(s,
2)*pow(t, 3) + 4320*pow(s, 3)*pow(t, 3) + 7524*pow(t, 4) - 22464*s*pow(t, 4) +
11664*pow(s, 2)*pow(t, 4) - 5796*pow(t, 5) + 9072*s*pow(t, 5) + 1728*pow(t,
6)))/4.; 

  d2phi(34,2) = (81*s*(5 - 16*s - 213*t + 630*s*t + 17*pow(s, 2) -
615*t*pow(s, 2) - 6*pow(s, 3) + 198*t*pow(s, 3) + 2238*pow(t, 2) - 5808*s*pow(t,
2) + 4860*pow(s, 2)*pow(t, 2) - 1296*pow(s, 3)*pow(t, 2) - 9490*pow(t, 3) +
19980*s*pow(t, 3) - 12600*pow(s, 2)*pow(t, 3) + 2160*pow(s, 3)*pow(t, 3) +
18810*pow(t, 4) - 28080*s*pow(t, 4) + 9720*pow(s, 2)*pow(t, 4) - 17388*pow(t, 5)
+ 13608*s*pow(t, 5) + 6048*pow(t, 6)))/2.; 

  d2phi(35,0) = (81*(-1 + 6*t)*pow(t,
2)*(67 - 528*s - 352*t + 2043*s*t + 1362*pow(s, 2) - 3456*t*pow(s, 2) -
1440*pow(s, 3) + 1800*t*pow(s, 3) + 540*pow(s, 4) + 681*pow(t, 2) -
2592*s*pow(t, 2) + 2160*pow(s, 2)*pow(t, 2) - 576*pow(t, 3) + 1080*s*pow(t, 3) +
180*pow(t, 4)))/2.; 

  d2phi(35,1) = (81*t*(20 - 268*s - 381*t + 4524*s*t +
1056*pow(s, 2) - 15633*t*pow(s, 2) - 1816*pow(s, 3) + 23256*t*pow(s, 3) +
1440*pow(s, 4) - 15660*t*pow(s, 4) - 432*pow(s, 5) + 3888*t*pow(s, 5) +
2312*pow(t, 2) - 22344*s*pow(t, 2) + 59400*pow(s, 2)*pow(t, 2) - 61056*pow(s,
3)*pow(t, 2) + 21600*pow(s, 4)*pow(t, 2) - 6415*pow(t, 3) + 46620*s*pow(t, 3) -
83160*pow(s, 2)*pow(t, 3) + 43200*pow(s, 3)*pow(t, 3) + 9036*pow(t, 4) -
43632*s*pow(t, 4) + 38880*pow(s, 2)*pow(t, 4) - 6300*pow(t, 5) + 15120*s*pow(t,
5) + 1728*pow(t, 6)))/4.; 

  d2phi(35,2) = (81*s*(10 - 67*s - 381*t + 2262*s*t +
176*pow(s, 2) - 5211*t*pow(s, 2) - 227*pow(s, 3) + 5814*t*pow(s, 3) + 144*pow(s,
4) - 3132*t*pow(s, 4) - 36*pow(s, 5) + 648*t*pow(s, 5) + 3468*pow(t, 2) -
16758*s*pow(t, 2) + 29700*pow(s, 2)*pow(t, 2) - 22896*pow(s, 3)*pow(t, 2) +
6480*pow(s, 4)*pow(t, 2) - 12830*pow(t, 3) + 46620*s*pow(t, 3) - 55440*pow(s,
2)*pow(t, 3) + 21600*pow(s, 3)*pow(t, 3) + 22590*pow(t, 4) - 54540*s*pow(t, 4) +
32400*pow(s, 2)*pow(t, 4) - 18900*pow(t, 5) + 22680*s*pow(t, 5) + 6048*pow(t,
6)))/2.; 

  d2phi(36,0) = (-1296*pow(t, 2)*(-97 + 1131*s + 754*t - 6795*s*t -
4530*pow(s, 2) + 19800*t*pow(s, 2) + 8250*pow(s, 3) - 23400*t*pow(s, 3) -
7020*pow(s, 4) + 9720*t*pow(s, 4) + 2268*pow(s, 5) - 2265*pow(t, 2) +
14850*s*pow(t, 2) - 28080*pow(s, 2)*pow(t, 2) + 16200*pow(s, 3)*pow(t, 2) +
3300*pow(t, 3) - 14040*s*pow(t, 3) + 12960*pow(s, 2)*pow(t, 3) - 2340*pow(t, 4)
+ 4860*s*pow(t, 4) + 648*pow(t, 5)))/25.; 

  d2phi(36,1) = (-648*t*(20 - 388*s -
291*t + 4524*s*t + 2262*pow(s, 2) - 20385*t*pow(s, 2) - 6040*pow(s, 3) +
39600*t*pow(s, 3) + 8250*pow(s, 4) - 35100*t*pow(s, 4) - 5616*pow(s, 5) +
11664*t*pow(s, 5) + 1512*pow(s, 6) + 1508*pow(t, 2) - 18120*s*pow(t, 2) +
59400*pow(s, 2)*pow(t, 2) - 74880*pow(s, 3)*pow(t, 2) + 32400*pow(s, 4)*pow(t,
2) - 3775*pow(t, 3) + 33000*s*pow(t, 3) - 70200*pow(s, 2)*pow(t, 3) +
43200*pow(s, 3)*pow(t, 3) + 4950*pow(t, 4) - 28080*s*pow(t, 4) + 29160*pow(s,
2)*pow(t, 4) - 3276*pow(t, 5) + 9072*s*pow(t, 5) + 864*pow(t, 6)))/25.;


  d2phi(36,2) = (-1296*s*(10 - 97*s - 291*t + 2262*s*t + 377*pow(s, 2) -
6795*t*pow(s, 2) - 755*pow(s, 3) + 9900*t*pow(s, 3) + 825*pow(s, 4) -
7020*t*pow(s, 4) - 468*pow(s, 5) + 1944*t*pow(s, 5) + 108*pow(s, 6) +
2262*pow(t, 2) - 13590*s*pow(t, 2) + 29700*pow(s, 2)*pow(t, 2) - 28080*pow(s,
3)*pow(t, 2) + 9720*pow(s, 4)*pow(t, 2) - 7550*pow(t, 3) + 33000*s*pow(t, 3) -
46800*pow(s, 2)*pow(t, 3) + 21600*pow(s, 3)*pow(t, 3) + 12375*pow(t, 4) -
35100*s*pow(t, 4) + 24300*pow(s, 2)*pow(t, 4) - 9828*pow(t, 5) + 13608*s*pow(t,
5) + 3024*pow(t, 6)))/25.; 

  d2phi(37,0) = (-1296*t*(10 - 291*s - 97*t + 2262*s*t
+ 2262*pow(s, 2) - 13590*t*pow(s, 2) - 7550*pow(s, 3) + 33000*t*pow(s, 3) +
12375*pow(s, 4) - 35100*t*pow(s, 4) - 9828*pow(s, 5) + 13608*t*pow(s, 5) +
3024*pow(s, 6) + 377*pow(t, 2) - 6795*s*pow(t, 2) + 29700*pow(s, 2)*pow(t, 2) -
46800*pow(s, 3)*pow(t, 2) + 24300*pow(s, 4)*pow(t, 2) - 755*pow(t, 3) +
9900*s*pow(t, 3) - 28080*pow(s, 2)*pow(t, 3) + 21600*pow(s, 3)*pow(t, 3) +
825*pow(t, 4) - 7020*s*pow(t, 4) + 9720*pow(s, 2)*pow(t, 4) - 468*pow(t, 5) +
1944*s*pow(t, 5) + 108*pow(t, 6)))/25.; 

  d2phi(37,1) = (-648*s*(20 - 291*s -
388*t + 4524*s*t + 1508*pow(s, 2) - 18120*t*pow(s, 2) - 3775*pow(s, 3) +
33000*t*pow(s, 3) + 4950*pow(s, 4) - 28080*t*pow(s, 4) - 3276*pow(s, 5) +
9072*t*pow(s, 5) + 864*pow(s, 6) + 2262*pow(t, 2) - 20385*s*pow(t, 2) +
59400*pow(s, 2)*pow(t, 2) - 70200*pow(s, 3)*pow(t, 2) + 29160*pow(s, 4)*pow(t,
2) - 6040*pow(t, 3) + 39600*s*pow(t, 3) - 74880*pow(s, 2)*pow(t, 3) +
43200*pow(s, 3)*pow(t, 3) + 8250*pow(t, 4) - 35100*s*pow(t, 4) + 32400*pow(s,
2)*pow(t, 4) - 5616*pow(t, 5) + 11664*s*pow(t, 5) + 1512*pow(t, 6)))/25.;


  d2phi(37,2) = (-1296*pow(s, 2)*(-97 + 754*s + 1131*t - 6795*s*t - 2265*pow(s, 2)
+ 14850*t*pow(s, 2) + 3300*pow(s, 3) - 14040*t*pow(s, 3) - 2340*pow(s, 4) +
4860*t*pow(s, 4) + 648*pow(s, 5) - 4530*pow(t, 2) + 19800*s*pow(t, 2) -
28080*pow(s, 2)*pow(t, 2) + 12960*pow(s, 3)*pow(t, 2) + 8250*pow(t, 3) -
23400*s*pow(t, 3) + 16200*pow(s, 2)*pow(t, 3) - 7020*pow(t, 4) + 9720*s*pow(t,
4) + 2268*pow(t, 5)))/25.; 

  d2phi(38,0) = (81*t*(10 - 381*s - 67*t + 2262*s*t +
3468*pow(s, 2) - 16758*t*pow(s, 2) - 12830*pow(s, 3) + 46620*t*pow(s, 3) +
22590*pow(s, 4) - 54540*t*pow(s, 4) - 18900*pow(s, 5) + 22680*t*pow(s, 5) +
6048*pow(s, 6) + 176*pow(t, 2) - 5211*s*pow(t, 2) + 29700*pow(s, 2)*pow(t, 2) -
55440*pow(s, 3)*pow(t, 2) + 32400*pow(s, 4)*pow(t, 2) - 227*pow(t, 3) +
5814*s*pow(t, 3) - 22896*pow(s, 2)*pow(t, 3) + 21600*pow(s, 3)*pow(t, 3) +
144*pow(t, 4) - 3132*s*pow(t, 4) + 6480*pow(s, 2)*pow(t, 4) - 36*pow(t, 5) +
648*s*pow(t, 5)))/2.; 

  d2phi(38,1) = (81*s*(20 - 381*s - 268*t + 4524*s*t +
2312*pow(s, 2) - 22344*t*pow(s, 2) - 6415*pow(s, 3) + 46620*t*pow(s, 3) +
9036*pow(s, 4) - 43632*t*pow(s, 4) - 6300*pow(s, 5) + 15120*t*pow(s, 5) +
1728*pow(s, 6) + 1056*pow(t, 2) - 15633*s*pow(t, 2) + 59400*pow(s, 2)*pow(t, 2)
- 83160*pow(s, 3)*pow(t, 2) + 38880*pow(s, 4)*pow(t, 2) - 1816*pow(t, 3) +
23256*s*pow(t, 3) - 61056*pow(s, 2)*pow(t, 3) + 43200*pow(s, 3)*pow(t, 3) +
1440*pow(t, 4) - 15660*s*pow(t, 4) + 21600*pow(s, 2)*pow(t, 4) - 432*pow(t, 5) +
3888*s*pow(t, 5)))/4.; 

  d2phi(38,2) = (81*(-1 + 6*s)*pow(s, 2)*(67 - 352*s -
528*t + 2043*s*t + 681*pow(s, 2) - 2592*t*pow(s, 2) - 576*pow(s, 3) +
1080*t*pow(s, 3) + 180*pow(s, 4) + 1362*pow(t, 2) - 3456*s*pow(t, 2) +
2160*pow(s, 2)*pow(t, 2) - 1440*pow(t, 3) + 1800*s*pow(t, 3) + 540*pow(t,
4)))/2.; 

  d2phi(39,0) = (81*t*(5 - 213*s - 16*t + 630*s*t + 2238*pow(s, 2) -
5808*t*pow(s, 2) - 9490*pow(s, 3) + 19980*t*pow(s, 3) + 18810*pow(s, 4) -
28080*t*pow(s, 4) - 17388*pow(s, 5) + 13608*t*pow(s, 5) + 6048*pow(s, 6) +
17*pow(t, 2) - 615*s*pow(t, 2) + 4860*pow(s, 2)*pow(t, 2) - 12600*pow(s,
3)*pow(t, 2) + 9720*pow(s, 4)*pow(t, 2) - 6*pow(t, 3) + 198*s*pow(t, 3) -
1296*pow(s, 2)*pow(t, 3) + 2160*pow(s, 3)*pow(t, 3)))/2.; 

  d2phi(39,1) =
(81*s*(10 - 213*s - 64*t + 1260*s*t + 1492*pow(s, 2) - 7744*t*pow(s, 2) -
4745*pow(s, 3) + 19980*t*pow(s, 3) + 7524*pow(s, 4) - 22464*t*pow(s, 4) -
5796*pow(s, 5) + 9072*t*pow(s, 5) + 1728*pow(s, 6) + 102*pow(t, 2) -
1845*s*pow(t, 2) + 9720*pow(s, 2)*pow(t, 2) - 18900*pow(s, 3)*pow(t, 2) +
11664*pow(s, 4)*pow(t, 2) - 48*pow(t, 3) + 792*s*pow(t, 3) - 3456*pow(s,
2)*pow(t, 3) + 4320*pow(s, 3)*pow(t, 3)))/4.; 

  d2phi(39,2) = (81*(-1 + 2*s)*(-1 +
3*s)*(-1 + 6*s)*pow(s, 2)*(16 - 34*s - 51*t + 54*s*t + 18*pow(s, 2) + 36*pow(t,
2)))/2.; 

  d2phi(40,0) = (-1296*t*(2 - 87*s - 4*t + 162*s*t + 942*pow(s, 2) -
1560*t*pow(s, 2) - 4150*pow(s, 3) + 5700*t*pow(s, 3) + 8595*pow(s, 4) -
8640*t*pow(s, 4) - 8316*pow(s, 5) + 4536*t*pow(s, 5) + 3024*pow(s, 6) + 2*pow(t,
2) - 75*s*pow(t, 2) + 630*pow(s, 2)*pow(t, 2) - 1800*pow(s, 3)*pow(t, 2) +
1620*pow(s, 4)*pow(t, 2)))/25.; 

  d2phi(40,1) = (-648*s*(4 - 87*s - 16*t + 324*s*t
+ 628*pow(s, 2) - 2080*t*pow(s, 2) - 2075*pow(s, 3) + 5700*t*pow(s, 3) +
3438*pow(s, 4) - 6912*t*pow(s, 4) - 2772*pow(s, 5) + 3024*t*pow(s, 5) +
864*pow(s, 6) + 12*pow(t, 2) - 225*s*pow(t, 2) + 1260*pow(s, 2)*pow(t, 2) -
2700*pow(s, 3)*pow(t, 2) + 1944*pow(s, 4)*pow(t, 2)))/25.; 

  d2phi(40,2) =
(-1296*(-1 + 2*s)*(-2 + 3*s)*(-1 + 3*s)*(-1 + 6*s)*(-2 + 2*s + 3*t)*pow(s,
2))/25.; 

  d2phi(41,0) = (648*(-2 + 81*s + 2*t - 75*s*t - 780*pow(s, 2) +
630*t*pow(s, 2) + 2850*pow(s, 3) - 1800*t*pow(s, 3) - 4320*pow(s, 4) +
1620*t*pow(s, 4) + 2268*pow(s, 5))*pow(t, 2)*sqrt(2))/25.; 

  d2phi(41,1) =
(324*s*t*(-8 + 162*s + 12*t - 225*s*t - 1040*pow(s, 2) + 1260*t*pow(s, 2) +
2850*pow(s, 3) - 2700*t*pow(s, 3) - 3456*pow(s, 4) + 1944*t*pow(s, 4) +
1512*pow(s, 5))*sqrt(2))/25.; 

  d2phi(41,2) = (648*(-1 + 2*s)*(-2 + 3*s)*(-1 +
3*s)*(-1 + 6*s)*(-1 + s + 3*t)*pow(s, 2)*sqrt(2))/25.; 

  d2phi(42,0) = (81*(-1 +
6*t)*(1 - 36*s - t + 33*s*t + 282*pow(s, 2) - 216*t*pow(s, 2) - 720*pow(s, 3) +
360*t*pow(s, 3) + 540*pow(s, 4))*pow(t, 2)*pow(sqrt(2), -1))/2.; 

  d2phi(42,1) =
(81*s*t*(-4 + 72*s + 42*t - 747*s*t - 376*pow(s, 2) + 3816*t*pow(s, 2) +
720*pow(s, 3) - 7020*t*pow(s, 3) - 432*pow(s, 4) + 3888*t*pow(s, 4) - 48*pow(t,
2) + 792*s*pow(t, 2) - 3456*pow(s, 2)*pow(t, 2) + 4320*pow(s, 3)*pow(t,
2))*pow(sqrt(2), -1))/4.; 

  d2phi(42,2) = (81*(-1 + 2*s)*(-1 + 3*s)*(-1 +
6*s)*pow(s, 2)*(1 - s - 21*t + 18*s*t + 36*pow(t, 2))*pow(sqrt(2), -1))/2.;


  d2phi(43,0) = (81*(-1 + 2*t)*(-1 + 3*t)*(-1 + 6*t)*(1 - 21*s - t + 18*s*t +
36*pow(s, 2))*pow(t, 2)*pow(sqrt(2), -1))/2.; 

  d2phi(43,1) = (81*s*t*(-4 + 42*s +
72*t - 747*s*t - 48*pow(s, 2) + 792*t*pow(s, 2) - 376*pow(t, 2) + 3816*s*pow(t,
2) - 3456*pow(s, 2)*pow(t, 2) + 720*pow(t, 3) - 7020*s*pow(t, 3) + 4320*pow(s,
2)*pow(t, 3) - 432*pow(t, 4) + 3888*s*pow(t, 4))*pow(sqrt(2), -1))/4.;


  d2phi(43,2) = (81*(-1 + 6*s)*pow(s, 2)*(1 - s - 36*t + 33*s*t + 282*pow(t, 2) -
216*s*pow(t, 2) - 720*pow(t, 3) + 360*s*pow(t, 3) + 540*pow(t, 4))*pow(sqrt(2),
-1))/2.; 

  d2phi(44,0) = (648*(-1 + 3*s + t)*(-1 + 2*t)*(-2 + 3*t)*(-1 + 3*t)*(-1
+ 6*t)*pow(t, 2)*sqrt(2))/25; 

  d2phi(44,1) = (324*s*t*(-8 + 12*s + 162*t -
225*s*t - 1040*pow(t, 2) + 1260*s*pow(t, 2) + 2850*pow(t, 3) - 2700*s*pow(t, 3)
- 3456*pow(t, 4) + 1944*s*pow(t, 4) + 1512*pow(t, 5))*sqrt(2))/25; 

  d2phi(44,2) =
(648*pow(s, 2)*(-2 + 2*s + 81*t - 75*s*t - 780*pow(t, 2) + 630*s*pow(t, 2) +
2850*pow(t, 3) - 1800*s*pow(t, 3) - 4320*pow(t, 4) + 1620*s*pow(t, 4) +
2268*pow(t, 5))*sqrt(2))/25; 

  d2phi(45,0) = 5832*(-1 + 2*t)*(-1 + 3*t)*(-1 +
6*t)*pow(t, 2)*(1 - 6*s - 2*t + 6*s*t + 6*pow(s, 2) + pow(t, 2)); 

  d2phi(45,1) =
5832*s*t*(-2 + 6*s + 39*t - 108*s*t - 4*pow(s, 2) + 66*t*pow(s, 2) - 236*pow(t,
2) + 564*s*pow(t, 2) - 288*pow(s, 2)*pow(t, 2) + 595*pow(t, 3) - 1080*s*pow(t,
3) + 360*pow(s, 2)*pow(t, 3) - 648*pow(t, 4) + 648*s*pow(t, 4) + 252*pow(t, 5));


  d2phi(45,2) = 5832*pow(s, 2)*(-1 + 2*s + 39*t - 72*s*t - pow(s, 2) + 33*t*pow(s,
2) - 354*pow(t, 2) + 564*s*pow(t, 2) - 216*pow(s, 2)*pow(t, 2) + 1190*pow(t, 3)
- 1440*s*pow(t, 3) + 360*pow(s, 2)*pow(t, 3) - 1620*pow(t, 4) + 1080*s*pow(t, 4)
+ 756*pow(t, 5)); 

  d2phi(46,0) = -2592*(-1 + 3*t)*(-1 + 6*t)*pow(t, 2)*(-5 + 48*s
+ 16*t - 102*s*t - 102*pow(s, 2) + 108*t*pow(s, 2) + 60*pow(s, 3) - 17*pow(t, 2)
+ 54*s*pow(t, 2) + 6*pow(t, 3)); 

  d2phi(46,1) = -2592*s*t*(-10 + 48*s + 183*t -
801*s*t - 68*pow(s, 2) + 1026*t*pow(s, 2) + 30*pow(s, 3) - 405*t*pow(s, 3) -
1004*pow(t, 2) + 3672*s*pow(t, 2) - 3744*pow(s, 2)*pow(t, 2) + 1080*pow(s,
3)*pow(t, 2) + 2235*pow(t, 3) - 5805*s*pow(t, 3) + 3240*pow(s, 2)*pow(t, 3) -
2160*pow(t, 4) + 2916*s*pow(t, 4) + 756*pow(t, 5)); 

  d2phi(46,2) = -2592*pow(s,
2)*(-5 + 16*s + 183*t - 534*s*t - 17*pow(s, 2) + 513*t*pow(s, 2) + 6*pow(s, 3) -
162*t*pow(s, 3) - 1506*pow(t, 2) + 3672*s*pow(t, 2) - 2808*pow(s, 2)*pow(t, 2) +
648*pow(s, 3)*pow(t, 2) + 4470*pow(t, 3) - 7740*s*pow(t, 3) + 3240*pow(s,
2)*pow(t, 3) - 5400*pow(t, 4) + 4860*s*pow(t, 4) + 2268*pow(t, 5)); 

  d2phi(47,0)
= 2592*(-1 + 6*t)*pow(t, 2)*(10 - 141*s - 47*t + 492*s*t + 492*pow(s, 2) -
1134*t*pow(s, 2) - 630*pow(s, 3) + 720*t*pow(s, 3) + 270*pow(s, 4) + 82*pow(t,
2) - 567*s*pow(t, 2) + 648*pow(s, 2)*pow(t, 2) - 63*pow(t, 3) + 216*s*pow(t, 3)
+ 18*pow(t, 4)); 

  d2phi(47,1) = 2592*s*t*(-20 + 141*s + 321*t - 2007*s*t -
328*pow(s, 2) + 4086*t*pow(s, 2) + 315*pow(s, 3) - 3375*t*pow(s, 3) - 108*pow(s,
4) + 972*t*pow(s, 4) - 1456*pow(t, 2) + 7038*s*pow(t, 2) - 9936*pow(s, 2)*pow(t,
2) + 4320*pow(s, 3)*pow(t, 2) + 2775*pow(t, 3) - 9045*s*pow(t, 3) + 6480*pow(s,
2)*pow(t, 3) - 2376*pow(t, 4) + 3888*s*pow(t, 4) + 756*pow(t, 5)); 

  d2phi(47,2) =
2592*pow(s, 2)*(-10 + 47*s + 321*t - 1338*s*t - 82*pow(s, 2) + 2043*t*pow(s, 2)
+ 63*pow(s, 3) - 1350*t*pow(s, 3) - 18*pow(s, 4) + 324*t*pow(s, 4) - 2184*pow(t,
2) + 7038*s*pow(t, 2) - 7452*pow(s, 2)*pow(t, 2) + 2592*pow(s, 3)*pow(t, 2) +
5550*pow(t, 3) - 12060*s*pow(t, 3) + 6480*pow(s, 2)*pow(t, 3) - 5940*pow(t, 4) +
6480*s*pow(t, 4) + 2268*pow(t, 5)); 

  d2phi(48,0) = -5832*pow(t, 2)*(-10 + 201*s +
67*t - 1056*s*t - 1056*pow(s, 2) + 4086*t*pow(s, 2) + 2270*pow(s, 3) -
5760*t*pow(s, 3) - 2160*pow(s, 4) + 2700*t*pow(s, 4) + 756*pow(s, 5) -
176*pow(t, 2) + 2043*s*pow(t, 2) - 5184*pow(s, 2)*pow(t, 2) + 3600*pow(s,
3)*pow(t, 2) + 227*pow(t, 3) - 1728*s*pow(t, 3) + 2160*pow(s, 2)*pow(t, 3) -
144*pow(t, 4) + 540*s*pow(t, 4) + 36*pow(t, 5)); 

  d2phi(48,1) = -5832*s*t*(-20 +
201*s + 201*t - 1584*s*t - 704*pow(s, 2) + 4086*t*pow(s, 2) + 1135*pow(s, 3) -
4320*t*pow(s, 3) - 864*pow(s, 4) + 1620*t*pow(s, 4) + 252*pow(s, 5) - 704*pow(t,
2) + 4086*s*pow(t, 2) - 6912*pow(s, 2)*pow(t, 2) + 3600*pow(s, 3)*pow(t, 2) +
1135*pow(t, 3) - 4320*s*pow(t, 3) + 3600*pow(s, 2)*pow(t, 3) - 864*pow(t, 4) +
1620*s*pow(t, 4) + 252*pow(t, 5)); 

  d2phi(48,2) = -5832*pow(s, 2)*(-10 + 67*s +
201*t - 1056*s*t - 176*pow(s, 2) + 2043*t*pow(s, 2) + 227*pow(s, 3) -
1728*t*pow(s, 3) - 144*pow(s, 4) + 540*t*pow(s, 4) + 36*pow(s, 5) - 1056*pow(t,
2) + 4086*s*pow(t, 2) - 5184*pow(s, 2)*pow(t, 2) + 2160*pow(s, 3)*pow(t, 2) +
2270*pow(t, 3) - 5760*s*pow(t, 3) + 3600*pow(s, 2)*pow(t, 3) - 2160*pow(t, 4) +
2700*s*pow(t, 4) + 756*pow(t, 5)); 

  d2phi(49,0) = 2592*pow(t, 2)*(-10 + 321*s +
47*t - 1338*s*t - 2184*pow(s, 2) + 7038*t*pow(s, 2) + 5550*pow(s, 3) -
12060*t*pow(s, 3) - 5940*pow(s, 4) + 6480*t*pow(s, 4) + 2268*pow(s, 5) -
82*pow(t, 2) + 2043*s*pow(t, 2) - 7452*pow(s, 2)*pow(t, 2) + 6480*pow(s,
3)*pow(t, 2) + 63*pow(t, 3) - 1350*s*pow(t, 3) + 2592*pow(s, 2)*pow(t, 3) -
18*pow(t, 4) + 324*s*pow(t, 4)); 

  d2phi(49,1) = 2592*s*t*(-20 + 321*s + 141*t -
2007*s*t - 1456*pow(s, 2) + 7038*t*pow(s, 2) + 2775*pow(s, 3) - 9045*t*pow(s, 3)
- 2376*pow(s, 4) + 3888*t*pow(s, 4) + 756*pow(s, 5) - 328*pow(t, 2) +
4086*s*pow(t, 2) - 9936*pow(s, 2)*pow(t, 2) + 6480*pow(s, 3)*pow(t, 2) +
315*pow(t, 3) - 3375*s*pow(t, 3) + 4320*pow(s, 2)*pow(t, 3) - 108*pow(t, 4) +
972*s*pow(t, 4)); 

  d2phi(49,2) = 2592*(-1 + 6*s)*pow(s, 2)*(10 - 47*s - 141*t +
492*s*t + 82*pow(s, 2) - 567*t*pow(s, 2) - 63*pow(s, 3) + 216*t*pow(s, 3) +
18*pow(s, 4) + 492*pow(t, 2) - 1134*s*pow(t, 2) + 648*pow(s, 2)*pow(t, 2) -
630*pow(t, 3) + 720*s*pow(t, 3) + 270*pow(t, 4)); 

  d2phi(50,0) = -2592*pow(t,
2)*(-5 + 183*s + 16*t - 534*s*t - 1506*pow(s, 2) + 3672*t*pow(s, 2) +
4470*pow(s, 3) - 7740*t*pow(s, 3) - 5400*pow(s, 4) + 4860*t*pow(s, 4) +
2268*pow(s, 5) - 17*pow(t, 2) + 513*s*pow(t, 2) - 2808*pow(s, 2)*pow(t, 2) +
3240*pow(s, 3)*pow(t, 2) + 6*pow(t, 3) - 162*s*pow(t, 3) + 648*pow(s, 2)*pow(t,
3)); 

  d2phi(50,1) = -2592*s*t*(-10 + 183*s + 48*t - 801*s*t - 1004*pow(s, 2) +
3672*t*pow(s, 2) + 2235*pow(s, 3) - 5805*t*pow(s, 3) - 2160*pow(s, 4) +
2916*t*pow(s, 4) + 756*pow(s, 5) - 68*pow(t, 2) + 1026*s*pow(t, 2) - 3744*pow(s,
2)*pow(t, 2) + 3240*pow(s, 3)*pow(t, 2) + 30*pow(t, 3) - 405*s*pow(t, 3) +
1080*pow(s, 2)*pow(t, 3)); 

  d2phi(50,2) = -2592*(-1 + 3*s)*(-1 + 6*s)*pow(s,
2)*(-5 + 16*s + 48*t - 102*s*t - 17*pow(s, 2) + 54*t*pow(s, 2) + 6*pow(s, 3) -
102*pow(t, 2) + 108*s*pow(t, 2) + 60*pow(t, 3)); 

  d2phi(51,0) = 5832*pow(t,
2)*(-1 + 39*s + 2*t - 72*s*t - 354*pow(s, 2) + 564*t*pow(s, 2) + 1190*pow(s, 3)
- 1440*t*pow(s, 3) - 1620*pow(s, 4) + 1080*t*pow(s, 4) + 756*pow(s, 5) - pow(t,
2) + 33*s*pow(t, 2) - 216*pow(s, 2)*pow(t, 2) + 360*pow(s, 3)*pow(t, 2));


  d2phi(51,1) = 5832*s*t*(-2 + 39*s + 6*t - 108*s*t - 236*pow(s, 2) + 564*t*pow(s,
2) + 595*pow(s, 3) - 1080*t*pow(s, 3) - 648*pow(s, 4) + 648*t*pow(s, 4) +
252*pow(s, 5) - 4*pow(t, 2) + 66*s*pow(t, 2) - 288*pow(s, 2)*pow(t, 2) +
360*pow(s, 3)*pow(t, 2)); 

  d2phi(51,2) = 5832*(-1 + 2*s)*(-1 + 3*s)*(-1 +
6*s)*pow(s, 2)*(1 - 2*s - 6*t + 6*s*t + pow(s, 2) + 6*pow(t, 2)); 

  d2phi(52,0) =
2592*(-1 + 6*t)*pow(t, 2)*(1 - 33*s - 2*t + 60*s*t + 222*pow(s, 2) -
324*t*pow(s, 2) - 450*pow(s, 3) + 360*t*pow(s, 3) + 270*pow(s, 4) + pow(t, 2) -
27*s*pow(t, 2) + 108*pow(s, 2)*pow(t, 2)); 

  d2phi(52,1) = 2592*s*t*(-2 + 33*s +
24*t - 387*s*t - 148*pow(s, 2) + 1656*t*pow(s, 2) + 225*pow(s, 3) -
2295*t*pow(s, 3) - 108*pow(s, 4) + 972*t*pow(s, 4) - 52*pow(t, 2) + 774*s*pow(t,
2) - 2736*pow(s, 2)*pow(t, 2) + 2160*pow(s, 3)*pow(t, 2) + 30*pow(t, 3) -
405*s*pow(t, 3) + 1080*pow(s, 2)*pow(t, 3)); 

  d2phi(52,2) = 2592*(-1 + 3*s)*(-1 +
6*s)*pow(s, 2)*(-1 + 2*s + 24*t - 42*s*t - pow(s, 2) + 18*t*pow(s, 2) -
78*pow(t, 2) + 72*s*pow(t, 2) + 60*pow(t, 3)); 

  d2phi(53,0) = 2592*(-1 + 3*t)*(-1
+ 6*t)*pow(t, 2)*(-1 + 24*s + 2*t - 42*s*t - 78*pow(s, 2) + 72*t*pow(s, 2) +
60*pow(s, 3) - pow(t, 2) + 18*s*pow(t, 2)); 

  d2phi(53,1) = 2592*s*t*(-2 + 24*s +
33*t - 387*s*t - 52*pow(s, 2) + 774*t*pow(s, 2) + 30*pow(s, 3) - 405*t*pow(s, 3)
- 148*pow(t, 2) + 1656*s*pow(t, 2) - 2736*pow(s, 2)*pow(t, 2) + 1080*pow(s,
3)*pow(t, 2) + 225*pow(t, 3) - 2295*s*pow(t, 3) + 2160*pow(s, 2)*pow(t, 3) -
108*pow(t, 4) + 972*s*pow(t, 4)); 

  d2phi(53,2) = 2592*(-1 + 6*s)*pow(s, 2)*(1 -
2*s - 33*t + 60*s*t + pow(s, 2) - 27*t*pow(s, 2) + 222*pow(t, 2) - 324*s*pow(t,
2) + 108*pow(s, 2)*pow(t, 2) - 450*pow(t, 3) + 360*s*pow(t, 3) + 270*pow(t, 4));


  d2phi(54,0) = -1458*(-1 + 6*t)*pow(t, 2)*(5 - 138*s - 16*t + 390*s*t +
678*pow(s, 2) - 1332*t*pow(s, 2) - 1080*pow(s, 3) + 1080*t*pow(s, 3) +
540*pow(s, 4) + 17*pow(t, 2) - 360*s*pow(t, 2) + 648*pow(s, 2)*pow(t, 2) -
6*pow(t, 3) + 108*s*pow(t, 3)); 

  d2phi(54,1) = -1458*s*t*(-10 + 138*s + 138*t -
1827*s*t - 452*pow(s, 2) + 5400*t*pow(s, 2) + 540*pow(s, 3) - 5670*t*pow(s, 3) -
216*pow(s, 4) + 1944*t*pow(s, 4) - 452*pow(t, 2) + 5400*s*pow(t, 2) -
11520*pow(s, 2)*pow(t, 2) + 6480*pow(s, 3)*pow(t, 2) + 540*pow(t, 3) -
5670*s*pow(t, 3) + 6480*pow(s, 2)*pow(t, 3) - 216*pow(t, 4) + 1944*s*pow(t, 4));


  d2phi(54,2) = -1458*(-1 + 6*s)*pow(s, 2)*(5 - 16*s - 138*t + 390*s*t + 17*pow(s,
2) - 360*t*pow(s, 2) - 6*pow(s, 3) + 108*t*pow(s, 3) + 678*pow(t, 2) -
1332*s*pow(t, 2) + 648*pow(s, 2)*pow(t, 2) - 1080*pow(t, 3) + 1080*s*pow(t, 3) +
540*pow(t, 4)); }

} }
