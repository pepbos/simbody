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
// The total cable path/span consists of LineSegments and CurveSegments
class SimTK_SIMBODY_EXPORT CurveSegment
{
public:
    CurveSegment()                                   = default;
    CurveSegment(const CurveSegment&)                = default;
    CurveSegment& operator=(const CurveSegment&)     = default;
    CurveSegment(CurveSegment&&) noexcept            = default;
    CurveSegment& operator=(CurveSegment&&) noexcept = default;
    ~CurveSegment()                                  = default;

    CurveSegment(
        CableSpan cable,
        const MobilizedBody& mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 xHint);

    enum class Status
    {
        Ok,
        Liftoff,
        Disabled,
    };

    // TODO remove? keep? make private?
    const CableSpan& getCable() const;

    Real getSegmentLength(const State& s) const;
    const Transform& getFrenetFrameStart(const State& s) const {
        throw std::runtime_error("NOTYETIMPLEMENTED");
    }
    const Transform& getFrenetFrameEnd(const State& s) const {
        throw std::runtime_error("NOTYETIMPLEMENTED");
    }

    Status getStatus(const State& state) const;

    const ContactGeometry& getContactGeometry() const;
    const Transform& getSurfaceToGroundTransform(const State& s) const;

    Real calcNormalCurvature(Vec3 point, Vec3 tangent) const {
        throw std::runtime_error("NOTYETIMPLEMENTED");
        return NaN;
    }

    Real calcGeodesicTorsion(Vec3 point, Vec3 tangent) const {
        throw std::runtime_error("NOTYETIMPLEMENTED");
        return NaN;
    }

    UnitVec3 calcSurfaceNormal(Vec3 point) const {
        throw std::runtime_error("NOTYETIMPLEMENTED");
        return {NaN, NaN, NaN};
    }

    Real calcSurfaceValue(Vec3 point) const {
        throw std::runtime_error("NOTYETIMPLEMENTED");
        return NaN;
    }

    // TODO seems useful:
    /* void setDisabled(const State& state) const; */
    /* void setEnabled(const State& state) const; */

    /* const MobilizedBody& getMobilizedBody(const State& state) const; */

    /* int calcSegmentPoints(const State& state, std::vector<Vec3>& points); */
    /* int calcSegmentFrenetFrames( */
    /*     const State& state, */
    /*     std::vector<ContactGeometry::FrenetFrame>& frames); */

    class Impl;

private:
    friend CableSpan;
    const Impl& getImpl() const {return *m_Impl;}
    Impl& updImpl() {return *m_Impl;}

    std::shared_ptr<Impl> m_Impl = nullptr;
};

//==============================================================================
//                          CABLE SPAN
//==============================================================================
class SimTK_SIMBODY_EXPORT CableSpan
{
public:
    struct LineSegment
    {
        LineSegment() = default;

        LineSegment(Vec3 a, Vec3 b): l((b-a).norm()), d((b-a)/l) {}

        Real l = NaN;
        UnitVec3 d{NaN, NaN, NaN};
    };

public:
    CableSpan(
        CableSubsystem& subsystem,
        const MobilizedBody& originBody,
        const Vec3& defaultOriginPoint,
        const MobilizedBody& terminationBody,
        const Vec3& defaultTerminationPoint);

    CurveSegmentIndex adoptSegment(const CurveSegment& segment);

    void adoptWrappingObstacle(
        const MobilizedBody& mobod,
        Transform X_BS,
        const ContactGeometry& geometry,
        Vec3 contactPointHint = {1., 0., 0.});

    int getNumCurveSegments() const;

    const CurveSegment& getCurveSegment(CurveSegmentIndex ix) const;

    Real getLength(const State& state) const;
    Real getLengthDot(const State& state) const;
    Real calcCablePower(const State& state, Real tension) const {return NaN;}

    void applyBodyForces(
        const State& state,
        Real tension,
        Vector_<SpatialVec>& bodyForcesInG) const;

    /* int calcPathPoints(const State& state, std::vector<Vec3>& points); */
    /* int calcPathFrenetFrames( */
    /*     const State& state, */
    /*     std::vector<ContactGeometry::FrenetFrame>& frames); */

    class Impl;

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

    friend CurveSegment::Impl;
    friend CableSubsystem;
};

//==============================================================================
//                                SUBSYSTEM
//==============================================================================

// Rename to CableTrackerSubSystem? WrappingPathSubsystem?
class SimTK_SIMBODY_EXPORT CableSubsystem : public Subsystem
{
public:
    CableSubsystem();
    explicit CableSubsystem(MultibodySystem&);

    int getNumPaths() const;
    const CableSpan& getPath(CableSpanIndex idx) const;
    CableSpan& updPath(CableSpanIndex idx);

    size_t writePathPoints(std::vector<Vec3>& points) const;
    size_t writePathFrames(std::vector<Transform>& frenetFrames) const;

    /* private: */
    SimTK_PIMPL_DOWNCAST(CableSubsystem, Subsystem);
    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;
};

} // namespace SimTK

#endif
