#include "SimTKmath.h"
#include "Wrapping.h"

using namespace SimTK;

//==============================================================================
//            SURFACE
//==============================================================================

Surface::Surface(const MobilizedBody& mobod, const Transform& X_BS, const ContactGeometry& geometry)
{

}

const ContactGeometry& getGeometry() const;
const Transform& calcSurfaceToGroundTransform(const State& state) const;
