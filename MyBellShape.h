#ifndef OOMPH_MY_BELL_ELEMENTS
#define OOMPH_MY_BELL_ELEMENTS

//oomph-lib headers
#include "generic/Vector.h"
#include "generic/shape.h"

namespace oomph {

namespace MyShape {

//===========================================================================
/// Two by two specialisation of function to calculate inverse of a matrix
//===========================================================================
inline double invert_two_by_two(const DenseMatrix<double> &jacobian,
                    DenseMatrix<double> &inverse_jacobian)
 {
  //Calculate the determinant of the matrix
  const double det = jacobian(0,0)*jacobian(1,1) - jacobian(0,1)*jacobian(1,0);

//Report if Matrix is singular or negative
#ifdef PARANOID
 if( fabs(det)<1e-12)
   {
    std::stringstream error_stream;
    error_stream << "The matrix is singular to machine precision : det(M) = "
         <<det<<".\n";
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }
#endif

  //Calculate the inverse of the 2x2 matrix
  inverse_jacobian(0,0) =  jacobian(1,1)/det;
  inverse_jacobian(0,1) = -jacobian(0,1)/det;
  inverse_jacobian(1,0) = -jacobian(1,0)/det;
  inverse_jacobian(1,1) =  jacobian(0,0)/det;

 return det;
 }

//=============================================================================
/// Three-by three specialisation of function to calculate inverse of a matrix
//=============================================================================
inline void invert_three_by_three(const DenseMatrix<double> &jacobian,
                    DenseMatrix<double> &inverse_jacobian)
 {
  //Calculate the determinant of the matrix
  const double det = jacobian(0,0)*jacobian(1,1)*jacobian(2,2)
   + jacobian(0,1)*jacobian(1,2)*jacobian(2,0)
   + jacobian(0,2)*jacobian(1,0)*jacobian(2,1)
   - jacobian(0,0)*jacobian(1,2)*jacobian(2,1)
   - jacobian(0,1)*jacobian(1,0)*jacobian(2,2)
   - jacobian(0,2)*jacobian(1,1)*jacobian(2,0);

  //Report if Matrix is singular or negative
#ifdef PARANOID
 if( fabs(det)<1e-12)
   {
    std::stringstream error_stream;
    error_stream << "The matrix is singular to machine precision : det(M) = "
      <<det<<".\n";
     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
   }
#endif

  //Calculate the inverse of the 3x3 matrix
  inverse_jacobian(0,0) =  (jacobian(1,1)*jacobian(2,2)
                            - jacobian(1,2)*jacobian(2,1))/det;
  inverse_jacobian(0,1) = -(jacobian(0,1)*jacobian(2,2)
                            - jacobian(0,2)*jacobian(2,1))/det;
  inverse_jacobian(0,2) =  (jacobian(0,1)*jacobian(1,2)
                            - jacobian(0,2)*jacobian(1,1))/det;
  inverse_jacobian(1,0) = -(jacobian(1,0)*jacobian(2,2)
                            - jacobian(1,2)*jacobian(2,0))/det;
  inverse_jacobian(1,1) =  (jacobian(0,0)*jacobian(2,2)
                            - jacobian(0,2)*jacobian(2,0))/det;
  inverse_jacobian(1,2) = -(jacobian(0,0)*jacobian(1,2)
                            - jacobian(0,2)*jacobian(1,0))/det;
  inverse_jacobian(2,0) =  (jacobian(1,0)*jacobian(2,1)
                            - jacobian(1,1)*jacobian(2,0))/det;
  inverse_jacobian(2,1) = -(jacobian(0,0)*jacobian(2,1)
                            - jacobian(0,1)*jacobian(2,0))/det;
  inverse_jacobian(2,2) =  (jacobian(0,0)*jacobian(1,1)
                            - jacobian(0,1)*jacobian(1,0))/det;
 }

// Get the (twice) area of the triangle from the vertices
double get_twice_triangle_area(const Vector<double>& v0, const Vector<double>&
  v1, const Vector<double>& v2);

// Get outer normal of side between vertices v0 and v1, assumes
// counter-clockwise triangle vertices.
Vector<double> get_outer_normal(const Vector<double>& v0, const
 Vector<double>& v1);

/// Basis on a reference element. This follows exactly the notation of M. Okabe
//  in Comput. Methods Appl. Mech. 117 (1994) 411-421
void d2_basis(const Vector<double>& s,const Vector<Vector<double> >& v,
  Shape& psi, DShape& dpsi, DShape& d2psi);


/// Basis on a reference element. This follows exactly the notation of M. Okabe
//  in Comput. Methods Appl. Mech. 117 (1994) 411-421
double d2_basis_eulerian(const Vector<double>& s,const Vector<Vector<double> >& v,
  Shape& psi, DShape& dpsi, DShape& d2psi);

}
}

#endif
