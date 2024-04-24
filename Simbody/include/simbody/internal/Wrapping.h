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

SimTK_DEFINE_UNIQUE_INDEX_TYPE(WrappingPathIndex);
SimTK_DEFINE_UNIQUE_INDEX_TYPE(WrapObstacleIndex);

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
    CurveSegment()                                 = default;
    CurveSegment(const CurveSegment&)              = delete;
    CurveSegment& operator = (const CurveSegment&) = delete;
    CurveSegment(CurveSegment&&) noexcept                  = default;
    CurveSegment& operator = (CurveSegment&&) noexcept     = default;
    ~CurveSegment()                                = default;

    CurveSegment(const MobilizedBody& mobod, Transform X_BS, const ContactGeometry& geometry, Vec3 xHint);
    /* bool isActive(const State& state) const; */

    class Impl;

    enum class Status
    {
        Ok,
        Liftoff,
        Disabled,
    };

private:
    explicit CurveSegment(std::unique_ptr<Impl> impl);

    friend CableSpan;
    const Impl& getImpl() const;
    Impl& updImpl();

    std::unique_ptr<Impl> impl;
};

//==============================================================================
//                          PATH - or SPAN
//==============================================================================
// Cable, wire, rope, cord
class SimTK_SIMBODY_EXPORT CableSpan
{
public:
        struct LineSegment
        {
            UnitVec3 d {NaN, NaN, NaN};
            Real l = NaN;
        };

public:
    CableSpan(
        CableSubsystem& subsystem,
        const MobilizedBody& originBody,
        const Vec3& defaultOriginPoint,
        const MobilizedBody& terminationBody,
        const Vec3& defaultTerminationPoint);

    std::vector<CurveSegment>& updObstacles();
    const std::vector<CurveSegment>& getObstacles();

    Real getLength(const State& state) const;
    Real getLengthDot(const State& state) const;

    class Impl;
private:
    explicit CableSpan(std::unique_ptr<Impl> impl);

    const Impl& getImpl() const { return *impl; }
    Impl& updImpl() { return *impl; }

    std::shared_ptr<Impl> impl = nullptr;

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
    const CableSpan& getPath(WrappingPathIndex idx) const;
    CableSpan& updPath(WrappingPathIndex idx);

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
