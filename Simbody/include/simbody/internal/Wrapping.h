#ifndef SimTK_SIMBODY_WRAPPING_PATH_SUBSYSTEM_H_
#define SimTK_SIMBODY_WRAPPING_PATH_SUBSYSTEM_H_

#include "SimTKcommon/internal/State.h"
#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/common.h"
#include "simmath/internal/ContactGeometry.h"
#include <functional>
#include <stdexcept>

namespace SimTK
{

SimTK_DEFINE_UNIQUE_INDEX_TYPE(CurveSegmentIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(CableSpanIndex);

class MultibodySystem;
class CableSubsystem;
class CurveSegment;
class CableSpan;

//==============================================================================
//                                CURVE SEGMENT
//==============================================================================
// A curved segment on a surface that is part of a `CableSpan`.
class SimTK_SIMBODY_EXPORT CurveSegment final
{
private:
    CurveSegment()                                   = default;

public:
    CurveSegment(const CurveSegment&)                = default;
    CurveSegment& operator=(const CurveSegment&)     = default;
    CurveSegment(CurveSegment&&) noexcept            = default;
    CurveSegment& operator=(CurveSegment&&) noexcept = default;
    ~CurveSegment()                                  = default;

    /** Construct a CurveSegment representing a segment of a CableSpan that
    wraps over a ContactGeometry.
    @param cable The cable this segment belongs to.
    @param mobod The body that the contact geometry is rigidly attached to.
    @param X_BS Transform specifying the location and orientation of the
    contact geometry's origin frame with respect to the mobilized body.
    @param geometry The contact geometry over which this segment wraps.
    @param initialContactPointHint A guess of the contact point of the cable
    span and the contact geometry to compute the initial cable path. This point
    is defined in the local contact geometry's frame. The point will be used as
    a starting point when computing the initial cable path. As such, it does
    not have to lie on the contact geometry's surface, nor does it have to
    belong to a valid cable path.*/
    CurveSegment(
        CableSpan cable,
        const MobilizedBody& mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 initialContactPointHint);

    /** A helper class, representing the wrapping status of this segment in
      relation to the contact geometry.

      Status::Ok indicates that the cable is in contact with the surface.
      Status::Lifted indicates that the cable is not in contact with the surface.
      Status::Disabled indicates that the surface obstacle is "disabled",
      preventing any interaction with the
      cable. */
    enum class Status
    {
        Ok,
        Lifted,
        Disabled,
    };

//------------------------------------------------------------------------------
//                Parameter interface
//------------------------------------------------------------------------------
    const CableSpan& getCable() const;

    const ContactGeometry& getContactGeometry() const;

    const Mobod& getMobilizedBody() const;

    const Transform& getXformSurfaceToBody() const;
    void setXformSurfaceToBody(Transform X_BS);

//------------------------------------------------------------------------------
//                State dependent getters.
//------------------------------------------------------------------------------
    Real getSegmentLength(const State& state) const;

    const Transform& getFrenetFrameStart(const State& state) const;

    const Transform& getFrenetFrameEnd(const State& state) const;

    Status getStatus(const State& state) const;

    int getNumberOfIntegratorStepsTaken(const State& state);

    Real getInitialIntegratorStepSize(const State& state);

    // TODO useful?
    /* void setDisabled(const State& state) const; */
    /* void setEnabled(const State& state) const; */

//------------------------------------------------------------------------------
//                State dependent computations.
//------------------------------------------------------------------------------
    void calcUnitForce(const State& state, SpatialVec& unitForce_G) const;

    // Compute the curve points in ground frame.
    //
    // Use `nPoints` to resample over the curve length at equal intervals using Hermite interpolation.
    // If `nPoints=0` no interpolation is applied, and the points from the numerical integration are used directly.
    //
    // Special cases that will override the `nPoints` argument:
    // - If the curve length is zero, a single point is written.
    // - If the curve is not active (disabled or lifted), no points are written.
    //
    // Note:
    // If `nPoints=1`, and the curve length is not zero, an exception is thrown.
    //
    // Returns the number of points written.
    int calcPoints(const State& state, std::vector<Vec3>& points_G, int nPoints = 0) const;

    bool isActive(const State& state) const
    {
        return getStatus(state) == Status::Ok;
    }

//------------------------------------------------------------------------------
    // TODO public? private?
    class Impl;

private:
    friend CableSpan;
    const Impl& getImpl() const
    {
        return *m_Impl;
    }
    Impl& updImpl()
    {
        return *m_Impl;
    }

    std::shared_ptr<Impl> m_Impl = nullptr;
};

//==============================================================================
//                          CABLE SPAN
//==============================================================================
// Representation of a cable spanning from one body to another.
// The cable can adopt obstacles that must be wrapped over.
//
// NOTE: The interaction with the obstacles is **ordered**. The cable can wrap
// over the obstacles in the order that they are stored. This means that the
// cable can not wrap over the first obstacle twice, even though it might
// spatially intersect it twice. This greatly simplifies the implementation,
// while covering many use-cases.
class SimTK_SIMBODY_EXPORT CableSpan final
{
public:

//------------------------------------------------------------------------------

    // Helper struct representing the cable segments that do not lie on a surface.
    struct LineSegment final
    {
        LineSegment() = default;

        LineSegment(Vec3 a, Vec3 b) : l((b - a).norm()), d((b - a) / l)
        {}

        Real l = NaN;
        UnitVec3 d{NaN, NaN, NaN};
    };

//------------------------------------------------------------------------------
    CableSpan()                                = default;
    ~CableSpan()                               = default;
    CableSpan(const CableSpan&)                = default;
    CableSpan& operator=(const CableSpan&)     = default;
    CableSpan(CableSpan&&) noexcept            = default;
    CableSpan& operator=(CableSpan&&) noexcept = default;

//------------------------------------------------------------------------------
//                Parameter interface
//------------------------------------------------------------------------------
    CableSpan(
        CableSubsystem& subsystem,
        const MobilizedBody& originBody,
        const Vec3& defaultOriginPoint,
        const MobilizedBody& terminationBody,
        const Vec3& defaultTerminationPoint);

    void adoptWrappingObstacle(
        const MobilizedBody& mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 contactPointHint = {1., 0., 0.});

    int getNumCurveSegments() const;

    const CurveSegment& getCurveSegment(CurveSegmentIndex ix) const;

//------------------------------------------------------------------------------
//                State dependent interface
//------------------------------------------------------------------------------
    Real getLength(const State& state) const;

    Real getLengthDot(const State& state) const;

    size_t countActive(const State& state) const;

//------------------------------------------------------------------------------
//                State dependent calculations
//------------------------------------------------------------------------------
    Real calcCablePower(const State& state, Real tension) const;

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    // Calls `CurveSegment::calcPoints` on each active curve segment.
    int calcPoints(const State& state, std::vector<Vec3>& points_G, int nPointsPerCurveSegment = 0) const;

//------------------------------------------------------------------------------
//                TODO TEMPORARY FOR UNIT TESTING
//------------------------------------------------------------------------------
    // TODO this is useful for unit testing, but we really really dont want to expose it...
    void calcPathErrorJacobian(
        const State& state,
        Vector& pathError,
        Matrix& pathErrorJacobian) const;

    // TODO this is useful for unit testing, but we really really dont want to expose it...
    void applyCorrection(const State& state, const Vector& correction) const;

//------------------------------------------------------------------------------
    class Impl; // TODO  private? public?

private:
    const Impl& getImpl() const
    {
        return *m_Impl;
    }

    Impl& updImpl()
    {
        return *m_Impl;
    }

    std::shared_ptr<Impl> m_Impl = nullptr;

    friend CurveSegment;
    friend CurveSegment::Impl;
    friend CableSubsystem;
};

//==============================================================================
//                                SUBSYSTEM
//==============================================================================

// TODO Alternative names: CableTrackerSubSystem, WrappingPathSubsystem.
class SimTK_SIMBODY_EXPORT CableSubsystem : public Subsystem
{
public:
    CableSubsystem();
    explicit CableSubsystem(MultibodySystem&);

    int getNumPaths() const;
    const CableSpan& getPath(CableSpanIndex idx) const;
    CableSpan& updPath(CableSpanIndex idx);

    /* private: */
    SimTK_PIMPL_DOWNCAST(CableSubsystem, Subsystem);
    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;
};

} // namespace SimTK

#endif
