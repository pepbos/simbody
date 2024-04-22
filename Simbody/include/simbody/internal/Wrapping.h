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
class WrappingPathSubsystem;
class WrapObstacle;
class WrappingPath;

//==============================================================================
//                      SURFACE
//==============================================================================
// Represents the local surface wrapping problem.
// Caches last computed geodesic as a warmstart.
// Not exposed outside of simbody.
// Not shared amongst different paths or obstacles.
class SimTK_SIMBODY_EXPORT Surface
{
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;

public:
    Surface()                              = default;
    ~Surface() = default;
    Surface(Surface&&) noexcept            = default;
    Surface& operator=(Surface&&) noexcept = default;
    Surface(const Surface&)                = delete;
    Surface& operator=(const Surface&)     = delete;

    Surface(WrappingPathSubsystem subsystem, const ContactGeometry& geometry, Vec3 xHint);

    // TODO move to impl to hide?
    struct LocalGeodesic
    {
        FrenetFrame KP {};
        FrenetFrame KQ {};

        Real length = NaN;

        Variation dKP {};
        Variation dKQ {};
    };

    // TODO move to impl to hide?
    struct WrappingStatus
    {
        Vec3 pointOnLine {NaN, NaN, NaN};

        bool liftoff = false;
        bool disabled = false;
    };

    const LocalGeodesic& calcGeodesic(State& s, Vec3 x, Vec3 t, Real l);
    const LocalGeodesic& applyGeodesicCorrection(const State& s, const Correction& c);
    const LocalGeodesic& getGeodesic(const State& s);
    Vec3 getPointOnLineNearSurface(const State& s);

    class Impl;
private:
    explicit Surface(std::unique_ptr<Impl> impl);
    const Impl& getImpl() const;
    Impl& updImpl();

    friend Impl;
    friend WrapObstacle;

    std::unique_ptr<Impl> impl = nullptr;
};

//==============================================================================
//                                OBSTACLE
//==============================================================================
// Although cheap to copy, we cannot hand them out because they have a cache entry associated with them.
// The surface they hold a pointer to can be reused in the model.
class SimTK_SIMBODY_EXPORT WrapObstacle
{
public:
    WrapObstacle()                                 = default;
    WrapObstacle(const WrapObstacle&)              = delete;
    WrapObstacle& operator = (const WrapObstacle&) = delete;
    WrapObstacle(WrapObstacle&&) noexcept                  = default;
    WrapObstacle& operator = (WrapObstacle&&) noexcept     = default;
    ~WrapObstacle()                                = default;

    WrapObstacle(const MobilizedBody& mobod, Transform X_BS, const ContactGeometry& geometry, Vec3 xHint);
    bool isActive(const State& state) const;

    class Impl;
private:
    explicit WrapObstacle(std::unique_ptr<Impl> impl);

    friend WrappingPath;
    const Impl& getImpl() const;
    Impl& updImpl();

    std::unique_ptr<Impl> impl;
};

//==============================================================================
//                                PATH
//==============================================================================
class SimTK_SIMBODY_EXPORT WrappingPath
{
public:
        struct LineSegment
        {
            UnitVec3 d {NaN, NaN, NaN};
            Real l = NaN;
        };

public:
    WrappingPath(
        WrappingPathSubsystem& subsystem,
        const MobilizedBody& originBody,
        const Vec3& defaultOriginPoint,
        const MobilizedBody& terminationBody,
        const Vec3& defaultTerminationPoint);

    std::vector<WrapObstacle>& updObstacles();
    const std::vector<WrapObstacle>& getObstacles();

    Real getLength(const State& state) const;

    class Impl;
private:
    explicit WrappingPath(std::unique_ptr<Impl> impl);

    friend WrappingPathSubsystem;
    const Impl& getImpl() const { return *impl; }
    Impl& updImpl() { return *impl; }

    std::shared_ptr<Impl> impl = nullptr;
};

//==============================================================================
//                                SUBSYSTEM
//==============================================================================

class SimTK_SIMBODY_EXPORT WrappingPathSubsystem : public Subsystem
{
public:
    WrappingPathSubsystem();
    explicit WrappingPathSubsystem(MultibodySystem&);

    int getNumPaths() const;
    const WrappingPath& getPath(WrappingPathIndex idx) const;
    WrappingPath& updPath(WrappingPathIndex idx);

    size_t writePathPoints(std::vector<Vec3>& points) const;
    size_t writePathFrames(std::vector<Transform>& frenetFrames) const;

/* private: */
    SimTK_PIMPL_DOWNCAST(WrappingPathSubsystem, Subsystem);
    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;
};

} // namespace SimTK

#endif
