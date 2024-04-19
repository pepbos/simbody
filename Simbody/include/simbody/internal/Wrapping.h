#ifndef SimTK_SIMBODY_WRAPPING_PATH_SUBSYSTEM_H_
#define SimTK_SIMBODY_WRAPPING_PATH_SUBSYSTEM_H_

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
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

//==============================================================================
//                      SURFACE
//==============================================================================
class SurfaceImpl {

private:
    SurfaceImpl()                              = default;
public:
    ~SurfaceImpl() = default;
    SurfaceImpl(SurfaceImpl&&) noexcept            = default;
    SurfaceImpl& operator=(SurfaceImpl&&) noexcept = default;
    SurfaceImpl(const SurfaceImpl&)                = default;
    SurfaceImpl& operator=(const SurfaceImpl&)     = default;

    SurfaceImpl(MobilizedBody mobod, Transform X_BS, ContactGeometry geometry) :
        m_Mobod(std::move(mobod)),
        m_Offset(std::move(X_BS)),
        m_Geometry(geometry) {}

    const ContactGeometry& getGeometry() const {return m_Geometry;}
    Transform calcSurfaceToGroundTransform(const State& state) const
    { throw std::runtime_error("check transform order");
        return m_Mobod.getBodyTransform(state).compose(m_Offset);}

private:
    MobilizedBody m_Mobod;
    Transform m_Offset;
    ContactGeometry m_Geometry;
};

// Binds {ContactGeometry, Body, OffsetFrame}
// Cheap to copy, reuseable in the model.
class SimTK_SIMBODY_EXPORT Surface
{
private:
    Surface()                              = default;
public:
    ~Surface() = default;
    Surface(Surface&&) noexcept            = default;
    Surface& operator=(Surface&&) noexcept = default;
    Surface(const Surface&)                = default;
    Surface& operator=(const Surface&)     = default;

    Surface(const MobilizedBody& mobod, const Transform& X_BS, const ContactGeometry& geometry)
        : impl(std::make_shared<SurfaceImpl>(mobod, X_BS, geometry)) {}

    const ContactGeometry& getGeometry() const {return impl->getGeometry();}
    Transform calcSurfaceToGroundTransform(const State& state) const {return impl->calcSurfaceToGroundTransform(state);}

private:
    std::shared_ptr<SurfaceImpl> impl = nullptr;
};

//==============================================================================
//                                OBSTACLE
//==============================================================================
// Although cheap to copy, we cannot hand them out because they have a cache entry associated with them.
// The surface they hold a pointer to can be reused in the model.
class SimTK_SIMBODY_EXPORT WrapObstacle
{
public:
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;
    using BoundaryFrameVariation  = WrapObstacle::Correction;

private:
    WrapObstacle()                                      = default;

public:
    WrapObstacle(const WrapObstacle& source)            = default;
    WrapObstacle& operator=(const WrapObstacle& source) = default;
    ~WrapObstacle()                                     = default;

    explicit WrapObstacle(Surface surface);

    // Solution previously commputed, all in local surface frame coordinates.
    struct WarmStartInfo
    {
        FrenetFrame KP {};
        FrenetFrame KQ {};

        Real length = NaN;

        BoundaryFrameVariation dKP {};
        BoundaryFrameVariation dKQ {};

        std::vector<Vec3> points {};
    };

    // Ground frame solution.
    struct PosInfo
    {
        FrenetFrame KP {};
        FrenetFrame KQ {};

        Real length = NaN;

        BoundaryFrameVariation dKP {};
        BoundaryFrameVariation dKQ {};
    };

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    void realizeInstance(const State& state) const;
    void realizePosition(const State& state) const;
    void realizeVelocity(const State& state) const;
    void realizeAcceleration(const State& state) const;
    void invalidateTopology();

    const PosInfo& getPosInfo(const State& state) const;
    PosInfo& updPosInfo(const State& state) const;

    size_t writeGeodesicPoints(const State& state, std::vector<Vec3> points) const;
    void applyCorrection(const State& state) const;

private:
    const WarmStartInfo& getWarmStartInfo(const State& state) const;
    WarmStartInfo& updWarmStartInfo(const State& state) const;
    void calcPosInfo(PosInfo& posInfo) const;

    // Required for accessing the discrete variable?
    std::shared_ptr<WrappingPathSubsystem> subsystem = nullptr;

    std::vector<WrapObstacle> obstacles {};

    // TOPOLOGY CACHE (set during realizeTopology())
    DiscreteVariableIndex       warmStartInfoIx;
    DiscreteVariableIndex       posInfoIx;
};

//==============================================================================
//                                PATH
//==============================================================================
class SimTK_SIMBODY_EXPORT WrappingPath
{
public:
    WrappingPath(
        WrappingPathSubsystem& subsystem,
        const MobilizedBody& originBody,
        const Vec3& defaultOriginPoint,
        const MobilizedBody& terminationBody,
        const Vec3& defaultTerminationPoint);

    int getNumObstacles() const;
    const WrapObstacle& getObstacle(WrapObstacleIndex obstacleIx) const;
    WrapObstacleIndex adoptObstacle(WrapObstacle);

    Real getLength(const State& state) const;
    Real getLengthDot(const State& state) const;
    void applyBodyForces( const State& state, Real tension, Vector_<SpatialVec>& bodyForcesInG) const;
    Real calcPower(const State& state, Real tension) const;

    WrapObstacle& setDecorativeGeometry(const DecorativeGeometry& viz);
    WrapObstacle& setNearPoint(const Vec3& point);
    WrapObstacle& setContactPointHints(
        const Vec3& startHint,
        const Vec3& endHint);

    class Impl;

private:
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

    int getNumCablePaths() const;
    const WrappingPath& getPath(WrappingPathIndex idx) const;
    WrappingPath& updPath(WrappingPathIndex idx);

    /** @cond **/ // Hide from Doxygen.
    SimTK_PIMPL_DOWNCAST(WrappingPathSubsystem, Subsystem);
    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;
    /** @endcond **/
};

//==============================================================================
//                         PATH :: IMPL
//==============================================================================

class WrappingPath::Impl {
public:

    struct LineSegment
    {
        UnitVec3 d {NaN, NaN, NaN};
        Real l = NaN;
    };

    struct PosInfo
    {
        Vec3 xO {NaN, NaN, NaN};
        Vec3 xI {NaN, NaN, NaN};

        Real l = NaN;

        std::vector<LineSegment> lines;

        size_t loopIter = 0;
    };

    struct VizInfo
    {
        std::vector<std::pair<Vec3, UnitVec3>> points {};
    };

    struct VelInfo
    {
        Real lDot = NaN;
    };

    int getNumObstacles() const {return obstacles.size();}
    const WrapObstacle& getObstacle(WrapObstacleIndex ix) const;
    WrapObstacleIndex adoptObstacle(WrapObstacle& obstacle);

    void calcPath(State& state, bool preventLiftOff = false) const;
    void calcInitPath(State& state, std::function<Vec3(WrapObstacleIndex)> pointHints);

    Real getPathError(const State& state) const;
    Real getLength(const State& state) const;
    Real getLengthDot(const State& state) const;
    void applyBodyForces(const State& state, Real tension, Vector_<SpatialVec>& bodyForcesInG) const;
    Real calcCablePower(const State& state, Real tension) const;

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    void realizeInstance(const State& state) const;
    void realizePosition(const State& state) const;
    void realizeVelocity(const State& state) const;
    void realizeAcceleration(const State& state) const;
    /* void calcEventTriggerInfo (const State&, Array_<EventTriggerInfo>&) const; */
    /* void handleEvents (State&, Event::Cause, const Array_<EventId>& eventIds, const HandleEventsOptions& options, HandleEventsResults& results) const; */

    const PosInfo& getPosInfo(const State& state) const;
    const VelInfo& getVelInfo(const State& state) const;
    const VizInfo& getVizInfo(const State& state) const;

    PosInfo& updPosEntry(const State& state) const;
    VelInfo& updVelEntry(const State& state) const;
    VizInfo& updVizInfo(const State& state) const;

    void calcPosInfo(PosInfo& posInfo) const;
    void calcVelInfo(const PosInfo& posInfo, VelInfo& velInfo) const;
    void calcVizInfo(const PosInfo& posInfo, VizInfo& vizInfo) const;

    // Be sure to call this whenever you make a topology-level change to
    // the cable definition, like adding an obstacle or modifying one in
    // a significant way.
    void invalidateTopology()
    {   if (subsystem) subsystem->invalidateSubsystemTopologyCache(); }

private:
friend class WrappingPath;
    std::shared_ptr<WrappingPathSubsystem> subsystem = nullptr;
    WrappingPathIndex          index {};

    std::vector<WrapObstacle> obstacles {};

    // TOPOLOGY CACHE (set during realizeTopology())
    DiscreteVariableIndex       posInfoIx;
    DiscreteVariableIndex       velInfoIx;
    DiscreteVariableIndex       vizInfoIx;
};

} // namespace SimTK

#endif
