// Non--inline functions for BellBiharmonic elements
#include "foeppl_von_karman.h"

namespace oomph
{
namespace MyC1CurvedElements
{
namespace TLagrangeShape
{
// Define templates
template<unsigned DIM, unsigned NNODE> 
void shape_u(const Vector<double> &s, Shape &psi);

template<unsigned DIM, unsigned NNODE>
void dshape_u_local(const Vector<double> &s, Shape &psi, DShape &dpsids);

template<unsigned DIM, unsigned NNODE>
void local_coordinate_of_node(const unsigned& j, Vector<double>& s);
//=======================================================================
/// Shape function for specific TElement<2,2>
//=======================================================================
 template<>
 void shape_u<2,2>(const Vector<double> &s, Shape &psi)
   {
    psi[0] = s[0];
    psi[1] = s[1];
    psi[2] = 1.0-s[0]-s[1];
   }
//=======================================================================
/// Derivatives of shape functions for specific TElement<2,2>
//=======================================================================
 template<>
 void dshape_u_local<2,2>(const Vector<double> &s,
                    Shape &psi, DShape &dpsids)
   {
    shape_u<2,2>(s, psi);
    
    // Derivatives
    dpsids(0,0) = 1.0;
    dpsids(0,1) = 0.0;
    dpsids(1,0) = 0.0;
    dpsids(1,1) = 1.0;
    dpsids(2,0) = -1.0;
    dpsids(2,1) = -1.0;
   }

  template<>
  void local_coordinate_of_node<2,2>(const unsigned& j,
                                Vector<double>& s) 
   {
    s.resize(2);

    switch (j)
     {
     case 0:
      s[0]=1.0;
      s[1]=0.0;
      break;
      
     case 1:
      s[0]=0.0;
      s[1]=1.0;
      break;
      
     case 2:
      s[0]=0.0;
      s[1]=0.0;
      break;
      
     default:
      std::ostringstream error_message;
      error_message << "Element only has three nodes; called with node number " 
                    << j << std::endl;
      
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
   }
  

//=======================================================================
/// Shape function for specific TElement<2,3>
//=======================================================================
 template<>
 void shape_u<2,3>(const Vector<double> &s, Shape &psi)
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
 void dshape_u_local<2,3>(const Vector<double> &s,
                    Shape &psi, DShape &dpsids)
   {
 //ALH: Don't know why object qualifier is needed
 shape_u<2,3>(s, psi);

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
/// Return local coordinates of node j
//=======================================================================
  template<>
  void local_coordinate_of_node<2,3>(const unsigned& j,
                                Vector<double>& s) 
   {
    s.resize(2);

    switch (j)
     {
     case 0:
      s[0]=1.0;
      s[1]=0.0;
      break;
      
     case 1:
      s[0]=0.0;
      s[1]=1.0;
      break;
      
     case 2:
      s[0]=0.0;
      s[1]=0.0;
      break;
      
     case 3:
      s[0]=0.5;
      s[1]=0.5;
      break;
      
     case 4:
      s[0]=0.0;
      s[1]=0.5;
      break;
      
     case 5:
      s[0]=0.5;
      s[1]=0.0;
      break;
      
     default:
      std::ostringstream error_message;
      error_message << "Element only has six nodes; called with node number " 
                    << j << std::endl;
      
      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
   }

//=======================================================================
/// Shape function for specific TElement<2,4>
//=======================================================================
template<>
void shape_u<2,4>(const Vector<double> &s, Shape &psi)
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
void dshape_u_local<2,4>(const Vector<double> &s,
                   Shape &psi, DShape &dpsids)
{
  
 //ALH: Don't know why object qualifier is needed
 shape_u<2,4>(s, psi);
 
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
  template<>
  void local_coordinate_of_node<2,4>(const unsigned& j,
                                Vector<double>& s)
   {
    s.resize(2);

    switch (j)
    {
    case 0:
     s[0]=1.0;
     s[1]=0.0;
     break;

    case 1:
     s[0]=0.0;
     s[1]=1.0;
     break;

    case 2:
     s[0]=0.0;
     s[1]=0.0;
     break;

    case 3:
     s[0]=2.0/3.0;
     s[1]=1.0/3.0;
     break;

    case 4:
     s[0]=1.0/3.0;
     s[1]=2.0/3.0;
     break;

    case 5:
     s[0]=0.0;
     s[1]=2.0/3.0;
     break;

    case 6:
     s[0]=0.0;
     s[1]=1.0/3.0;
     break;

    case 8:
     s[0]=2.0/3.0;
     s[1]=0.0;
     break;

    case 7:
     s[0]=1.0/3.0;
     s[1]=0.0;
     break;

    case 9:
     s[0]=1.0/3.0;
     s[1]=1.0/3.0;
     break;

    default:
     std::ostringstream error_message;
     error_message << "Element only has ten nodes; called with node number " 
                   << j << std::endl;
     
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   }
}
}
//======================================================================
/// Set the data for the number of Variables at each node
//======================================================================
// template<unsigned DIM, unsigned NNODE_1D>
// const unsigned FoepplVonKarmanC1CurvedBellElement<DIM,NNODE_1D>::Initial_Nvalue[1] = 8;

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

//  template<unsigned BOUNDARY_ORDER>
// const unsigned FoepplVonKarmanC1CurvedBellElement<2,4,BOUNDARY_ORDER>::Initial_Nvalue[9] = {8,8,8,2,2,2,2,2,2,};

//=======================================================================
/// Shape function for specific TElement<DIM,NNODE,BOUNDARY_ORDER>
//=======================================================================
 template<unsigned DIM, unsigned NNODE, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE,BOUNDARY_ORDER>::shape_u(const Vector<double> &s, Shape &psi) const
   {
    // Call the template function 
    MyC1CurvedElements::TLagrangeShape::
     shape_u<DIM,NNODE>(s,psi);
   }
//=======================================================================
/// Derivatives of shape functions for specific TElement<2,2,BOUNDARY_ORDER>
//=======================================================================
 template<unsigned DIM, unsigned NNODE, unsigned BOUNDARY_ORDER>
 void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE,BOUNDARY_ORDER>::dshape_u_local(const Vector<double> &s,
                    Shape &psi, DShape &dpsids) const
   {
    // Call the template function 
    MyC1CurvedElements::TLagrangeShape::
     dshape_u_local<DIM,NNODE>(s,psi,dpsids);
  }

////=======================================================================
///// Return local coordinates of node j
////=======================================================================
// template<unsigned DIM, unsigned NNODE, unsigned BOUNDARY_ORDER>
// void FoepplVonKarmanC1CurvedBellElement<DIM,NNODE,BOUNDARY_ORDER>::
// local_coordinate_of_node(const unsigned& j, Vector<double>& s) const;

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
