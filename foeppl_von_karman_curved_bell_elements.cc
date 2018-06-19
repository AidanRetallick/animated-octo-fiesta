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
