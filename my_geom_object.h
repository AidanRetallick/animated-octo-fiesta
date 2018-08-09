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
#ifndef OOMPH_MY_GEOM_OBJECTS_HEADER
#define OOMPH_MY_GEOM_OBJECTS_HEADER
#include "generic/geom_objects.h"
namespace oomph
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
// My Geometric object
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// \short Specialisation of GeomObject that has a scalar parametric zeta and 
/// a 2D vector r coordinate r = (r_1(zeta),r_2(zeta)). 
/// This is suitable for use in meshing/construction of planar objects.
class CurvilineGeomObject : GeomObject

{
public:

 /// \short Constructor: Pass dimension of geometric object (# of Eulerian
 /// coords = # of Lagrangian coords; no time history available/needed)
 CurvilineGeomObject() : GeomObject(1,2)
  { }

 /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
 /// and pointer to time-stepper which is used to handle the
 /// position at previous timesteps and allows the evaluation
 /// of veloc/acceleration etc. in cases where the GeomData
 /// varies with time.
 CurvilineGeomObject(TimeStepper* time_stepper_pt) : GeomObject(1, 2, time_stepper_pt)
   { }

 /// Broken copy constructor
 CurvilineGeomObject(const CurvilineGeomObject& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineGeomObject");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineGeomObject&) 
  {
   BrokenCopy::broken_assign("CurvilineGeomObject");
  }

 /// (Empty) destructor
 virtual ~CurvilineGeomObject(){}

 /// Parametrised position on object at current time: r(zeta).
 virtual void position(const Vector<double>& zeta, Vector<double>& r) 
  const 
  {
   throw OomphLibError(
    "You must specify position() with a vector zeta for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// \short Position Vector w.r.t. to zeta: 
 virtual void position(const double& zeta, 
                        Vector<double> &r) const = 0;

 /// \short Derivative of position Vector w.r.t. to zeta: 
 virtual void dposition(const double& zeta, 
                        Vector<double> &drdzeta) const
  {
   throw OomphLibError(
    "You must specify dposition() for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }


 /// \short 2nd derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const double& zeta, 
                         Vector<double> &drdzeta) const
  {
   throw OomphLibError(
    "You must specify d2position() for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Get s from x for part 0 of the boundary (inverse mapping - for convenience)
 virtual double get_zeta(const Vector<double>& x)
 {
   throw OomphLibError(
    "You must specify get_zeta() for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
 }

};


/// \short Specialisation of CurvilineGeomObject for half a circle.
class CurvilineCircleTop : public CurvilineGeomObject
{
public:
 /// \short Constructor: Pass dimension of geometric object (# of Eulerian
 /// coords = # of Lagrangian coords; no time history available/needed)
 CurvilineCircleTop() : CurvilineGeomObject()
  { }

 /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
 /// and pointer to time-stepper which is used to handle the
 /// position at previous timesteps and allows the evaluation
 /// of veloc/acceleration etc. in cases where the GeomData
 /// varies with time.
 CurvilineCircleTop(TimeStepper* time_stepper_pt) : CurvilineGeomObject(time_stepper_pt)
   { }

 /// Broken copy constructor
 CurvilineCircleTop(const CurvilineCircleTop& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineCircleTop");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineCircleTop&) 
  {
   BrokenCopy::broken_assign("CurvilineCircleTop");
  }

 /// (Empty) destructor
 virtual ~CurvilineCircleTop(){}

 /// Parametrised position on object at current time: r(zeta).
 virtual void position(const Vector<double>& zeta, Vector<double>& r) 
  const 
  {
   throw OomphLibError(
    "You must specify position() with a vector zeta for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// \short Position Vector w.r.t. to zeta: 
 virtual void position(const double& zeta, 
                        Vector<double> &r) const
  { r[0] =-std::sin(zeta);  r[1] = std::cos(zeta);}

 /// \short Derivative of position Vector w.r.t. to zeta: 
 virtual void dposition(const double& zeta, 
                        Vector<double> &drdzeta) const
  { drdzeta[0] =-std::cos(zeta);  drdzeta[1] =-std::sin(zeta);}


 /// \short 2nd derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const double& zeta, 
                         Vector<double> &drdzeta) const
  { drdzeta[0] = std::sin(zeta);  drdzeta[1] =-std::cos(zeta);}

 /// Get s from x for part 0 of the boundary (inverse mapping - for convenience)
 double get_zeta(const Vector<double>& x)
 {
 // The arc length (parametric parameter) for the upper semi circular arc
  return atan2(-x[0],x[1]);
 }
  

};

/// \short Specialisation of CurvilineGeomObject for half a circle.
class CurvilineCircleBottom : public CurvilineGeomObject
{
public:
 /// \short Constructor: Pass dimension of geometric object (# of Eulerian
 /// coords = # of Lagrangian coords; no time history available/needed)
 CurvilineCircleBottom() : CurvilineGeomObject()
  { }

 /// \short Constructor: pass # of Eulerian and Lagrangian coordinates
 /// and pointer to time-stepper which is used to handle the
 /// position at previous timesteps and allows the evaluation
 /// of veloc/acceleration etc. in cases where the GeomData
 /// varies with time.
 CurvilineCircleBottom(TimeStepper* time_stepper_pt) : CurvilineGeomObject(time_stepper_pt)
   { }

 /// Broken copy constructor
 CurvilineCircleBottom(const CurvilineCircleBottom& dummy) 
  { 
   BrokenCopy::broken_copy("CurvilineCircleBottom");
  } 
 
 /// Broken assignment operator
 void operator=(const CurvilineCircleBottom&) 
  {
   BrokenCopy::broken_assign("CurvilineCircleBottom");
  }

 /// (Empty) destructor
 virtual ~CurvilineCircleBottom(){}

 /// Parametrised position on object at current time: r(zeta).
 virtual void position(const Vector<double>& zeta, Vector<double>& r) 
  const 
  {
   throw OomphLibError(
    "You must specify position() with a vector zeta for your own object! \n",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// \short Position Vector w.r.t. to zeta: 
 virtual void position(const double& zeta, 
                        Vector<double> &r) const
  { r[0] = std::sin(zeta);  r[1] =-std::cos(zeta);}

 /// \short Derivative of position Vector w.r.t. to zeta: 
 virtual void dposition(const double& zeta, 
                        Vector<double> &drdzeta) const
  { drdzeta[0] = std::cos(zeta);  drdzeta[1] = std::sin(zeta);}


 /// \short 2nd derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const double& zeta, 
                         Vector<double> &drdzeta) const
  { drdzeta[0] =-std::sin(zeta);  drdzeta[1] = std::cos(zeta);}

 /// Get s from x for part 0 of the boundary (inverse mapping - for convenience)
 double get_zeta(const Vector<double>& x)
 {
 // The arc length (parametric parameter) for the upper semi circular arc
  return atan2(x[0],-x[1]);
 }
  

};
}

#endif
