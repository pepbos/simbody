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
//                                SUBSYSTEM
//==============================================================================
class WrappingPath;

class SimTK_SIMBODY_EXPORT WrappingPathSubsystem : public Subsystem
{
public:
    WrappingPathSubsystem();
    explicit WrappingPathSubsystem(MultibodySystem&);

    int getNumPaths() const;
    const WrappingPath& getPath(WrappingPathIndex idx) const;
    WrappingPath& updPath(WrappingPathIndex idx);

    SimTK_PIMPL_DOWNCAST(WrappingPathSubsystem, Subsystem);
    class Impl;
    Impl& updImpl();
    const Impl& getImpl() const;
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

    class Impl;

private:
    friend WrappingPathSubsystem::Impl;

    const Impl& getImpl() const;
    Impl& updImpl();

    std::shared_ptr<Impl> impl = nullptr;
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

            // TODO solver matrices
        };

        std::vector<WrapObstacle>& updObstacles() {return m_Obstacles;}
        const std::vector<WrapObstacle>& getObstacles() const {return m_Obstacles;}

        // Allocate state variables and cache entries.
        void realizeTopology(State& state);
        void realizePosition(const State& state) const;
        void realizeVelocity(const State& state) const;
        void realizeAcceleration(const State& state) const;
        void invalidateTopology()
        {   if (m_Subsystem) m_Subsystem->invalidateSubsystemTopologyCache(); }

        const PosInfo& getPosInfo(const State& state) const;

    private:

        PosInfo& updPosInfo(const State& state) const;

        void calcPosInfo(PosInfo& posInfo) const;

        void calcPath(State& state, bool preventLiftOff = false) const;
        void calcInitPath(State& state, std::function<Vec3(WrapObstacleIndex)> pointHints);

        std::shared_ptr<WrappingPathSubsystem> m_Subsystem = nullptr;
        WrappingPathIndex          m_Index {};

        std::vector<WrapObstacle> m_Obstacles {};

        // TOPOLOGY CACHE (set during realizeTopology())
        CacheEntryIndex       m_PosInfoIx;
        CacheEntryIndex       m_VelInfoIx;
        CacheEntryIndex       m_VizInfoIx;

        friend class WrappingPath;
};

//==============================================================================
//                         SUBSYSTEM :: IMPL
//==============================================================================
class WrappingPathSubsystem::Impl : public Subsystem::Guts {
    public:
        Impl() {}
        ~Impl() {}
        Impl* cloneImpl() const override 
        {   return new Impl(*this); }

        int getNumPaths() const {return paths.size();}

        const WrappingPath& getCablePath(WrappingPathIndex index) const 
        {   return paths[index]; }

        WrappingPath& updCablePath(WrappingPathIndex index) 
        {   return paths[index]; }

        // Add a cable path to the list, bumping the reference count.
        WrappingPathIndex adoptCablePath(WrappingPath& path) {
            paths.push_back(path);
            return WrappingPathIndex(paths.size()-1);
        }

        // Return the MultibodySystem which owns this WrappingPathSubsystem.
        const MultibodySystem& getMultibodySystem() const 
        {   return MultibodySystem::downcast(getSystem()); }

        // Return the SimbodyMatterSubsystem from which this WrappingPathSubsystem
        // gets the bodies to track.
        const SimbodyMatterSubsystem& getMatterSubsystem() const 
        {   return getMultibodySystem().getMatterSubsystem(); }

        // Allocate state variables.
        int realizeSubsystemTopologyImpl(State& state) const override {
            // Briefly allow writing into the Topology cache; after this the
            // Topology cache is const.
            Impl* wThis = const_cast<Impl*>(this);

            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                WrappingPath& path = wThis->updCablePath(ix);
                path.updImpl().realizeTopology(state);
            }

            return 0;
        }

        int realizeSubsystemInstanceImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizeInstance(state);
            }
            return 0;
        }

        int realizeSubsystemPositionImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizePosition(state);
            }
            return 0;
        }

        int realizeSubsystemVelocityImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizeVelocity(state);
            }
            return 0;
        }


        int realizeSubsystemAccelerationImpl(const State& state) const override {
            for (WrappingPathIndex ix(0); ix < paths.size(); ++ix) {
                const WrappingPath& path = getCablePath(ix);
                path.getImpl().realizeAcceleration(state);
            }
            return 0;
        }

        SimTK_DOWNCAST(Impl, Subsystem::Guts);

    private:
        // TOPOLOGY STATE
        Array_<WrappingPath, WrappingPathIndex> paths;
};

} // namespace SimTK

#endif
