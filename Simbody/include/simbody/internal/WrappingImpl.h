#ifndef SimTK_SIMBODY_WRAPPING_PATH_IMPL_H_
#define SimTK_SIMBODY_WRAPPING_PATH_IMPL_H_

#include "SimTKcommon/internal/State.h"
#include "SimTKmath.h"
#include "simbody/internal/Wrapping.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/common.h"
#include "simmath/internal/ContactGeometry.h"

#include <functional>
#include <stdexcept>

namespace SimTK
{

//==============================================================================
//                      SURFACE IMPL
//==============================================================================
class Surface::Impl {
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;
    using PointOnLineResult = ContactGeometry::PointOnLineResult;

public:
    Impl()                              = default;
    ~Impl() = default;
    Impl(Impl&&) noexcept            = default;
    Impl& operator=(Impl&&) noexcept = default;

    Impl(const Impl&)                = delete;
    Impl& operator=(const Impl&)     = delete;

    struct GeodesicInitialConditions
    {
        Vec3 x {NaN, NaN, NaN};
        Vec3 t {NaN, NaN, NaN};
        Real l = NaN;
    };

    struct LocalGeodesicInfo : Surface::LocalGeodesic, Surface::WrappingStatus
    {
        std::vector<Vec3> points {};
        double sHint = NaN;
    };

    Impl(
            WrappingPathSubsystem subsystem,
            ContactGeometry geometry,
            Vec3 initPointGuess
            ) : 
        m_Subsystem(subsystem),
        m_Geometry(geometry),
        m_InitPointGuess(initPointGuess)
    {}

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    // TODO requires all levels?
    void realizeInstance(const State& state) const;
    void realizePosition(const State& state) const;
    void realizeVelocity(const State& state) const;
    void realizeAcceleration(const State& state) const;
    // Requires invalidateTopology?
    /* void invalidateTopology() */
    /* {   m_Subsystem.invalidateSubsystemTopologyCache(); } */

    // The virtual methods:
    // 1. calcGeodesic()
    // 2. calcPointOnLineNearSurface
    // 3. isPointAboveSurface
    // 4. calcPathPoints

    // 1. calcGeodesic
    const LocalGeodesic& calcGeodesic(State& state, Vec3 x, Vec3 t, Real l) const;
    /* void calcGeodesic(GeodesicInitialConditions g0, Real sHint, LocalGeodesic& geodesic) const; */
    /* void calcGeodesic(Vec3 x, Vec3 t, Real l, Real sHint, LocalGeodesic& geodesic) const; */
    const LocalGeodesic& getGeodesic(const State& state) const {return getLocalGeodesicInfo(state);}
    void applyGeodesicCorrection(const State& state, const Correction& c) const;

    // TODO rename
    bool isActive(const State& state) const;

    // 4. calcPathPoints
    size_t calcPathPoints(const State& state, std::vector<Vec3>& points) const;

    // Detect liftoff or touchdown, and compute the point on line near the surface.
    const WrappingStatus& calcWrappingStatus(const State& state, Vec3 prev, Vec3 next, size_t maxIter, Real eps) const;

    // TODO should be in OpenSim?
    void setInitialPointGuess(Vec3 pointGuess);
    Vec3 getInitialPointGuess() const;

private:
    const LocalGeodesicInfo& getLocalGeodesicInfo(const State& state) const
    {
        realizePosition(state);
        return Value<LocalGeodesicInfo>::downcast(
                m_Subsystem.getDiscreteVarUpdateValue(state, m_GeodesicInfoIx)
                );
    }

    LocalGeodesicInfo& updLocalGeodesicInfo(const State& state) const
    {
        return Value<LocalGeodesicInfo>::updDowncast(
                m_Subsystem.updDiscreteVarUpdateValue(state, m_GeodesicInfoIx)
                );
    }

    const LocalGeodesicInfo& getPrevLocalGeodesicInfo(const State& state) const
    {
        return Value<LocalGeodesicInfo>::downcast(
                m_Subsystem.getDiscreteVariable(state, m_GeodesicInfoIx)
                );
    }

    LocalGeodesicInfo& updPrevLocalGeodesicInfo(State& state) const
    {
        return Value<LocalGeodesicInfo>::updDowncast(
                m_Subsystem.updDiscreteVariable(state, m_GeodesicInfoIx)
                );
    }

    void calcLocalGeodesicInfo(Vec3 x, Vec3 t, Real l, Real sHint,
            LocalGeodesicInfo& geodesic) const;

//------------------------------------------------------------------------------
    WrappingPathSubsystem m_Subsystem;

    ContactGeometry m_Geometry;

    Vec3 m_InitPointGuess;

    DiscreteVariableIndex m_GeodesicInfoIx;
};

//==============================================================================
//                                OBSTACLE IMPL
//==============================================================================
class WrapObstacle::Impl
{
    using FrenetFrame = ContactGeometry::FrenetFrame;
    using Variation = ContactGeometry::GeodesicVariation;
    using Correction = ContactGeometry::GeodesicCorrection;

private:
    Impl()                                      = default;

public:
    Impl(const Impl& source)            = default;
    Impl& operator=(const Impl& source) = default;
    ~Impl()                                     = default;

    /* Surface(const MobilizedBody& mobod, const Transform& X_BS, const ContactGeometry& geometry) */
    /*     : impl(std::make_shared<SurfaceImpl>(mobod, X_BS, geometry)) {} */
    /* Transform calcSurfaceToGroundTransform(const State& state) const {return impl->calcSurfaceToGroundTransform(state);} */

    Impl(
            WrappingPathSubsystem subsystem,
            const MobilizedBody& mobod,
            const Transform& X_BS,
            ContactGeometry geometry,
            Vec3 initPointGuess
            ) : 
        m_Subsystem(subsystem),
        m_Mobod(mobod),
        m_Offset(X_BS),
        m_Surface(subsystem, geometry, initPointGuess)
    {}

    enum class Status
    {
        Ok,
        Liftoff,
        Disabled,
    };

    // Ground frame solution.
    struct PosInfo
    {
        FrenetFrame KP {};
        FrenetFrame KQ {};

        Variation dKP {};
        Variation dKQ {};

        Real length = NaN;

        Status status = Status::Ok;
    };

    // Allocate state variables and cache entries.
    void realizeTopology(State& state);
    /* void realizeInstance(const State& state) const; */
    void realizePosition(const State& state) const;
    /* void realizeVelocity(const State& state) const; */
    /* void realizeAcceleration(const State& state) const; */
    /* void invalidateTopology(); */

    const PosInfo& getPosInfo(const State& state) const;

    bool isActive(const State& state) const;
    bool isDisabled(const State& state) const;

    PosInfo& updPosInfo(const State &state) const
    {
        return Value<PosInfo>::updDowncast(m_Subsystem.updCacheEntry(state, m_PosInfoIx));
    }

    size_t calcPathPoints(const State& state, std::vector<Vec3>& points) const;
    void applyGeodesicCorrection(const State& state, const ContactGeometry::GeodesicCorrection& c) const;

    void calcPosInfo(const State& state, Vec3 prev, Vec3 next, size_t maxIter, Real eps) const;
    void calcGeodesic(State& state, Vec3 x, Vec3 t, Real l) const;

private:
    void calcPosInfo(const State& state, PosInfo& posInfo) const;

    // TODO Required for accessing the cache variable?
    WrappingPathSubsystem m_Subsystem;

    MobilizedBody m_Mobod;
    Transform m_Offset;

    Surface m_Surface;

    // TOPOLOGY CACHE
    CacheEntryIndex       m_PosInfoIx;
};


//==============================================================================
//                         PATH :: IMPL
//==============================================================================
class WrappingPath::Impl {
    public:
        Impl(
                WrappingPathSubsystem subsystem,
                MobilizedBody originBody,
                Vec3 originPoint,
                MobilizedBody terminationBody,
                Vec3 terminationPoint):
            m_Subsystem(subsystem),
            m_OriginBody(originBody),
            m_OriginPoint(originPoint),
            m_TerminationBody(terminationBody),
            m_TerminationPoint(terminationPoint) {}

        struct PosInfo
        {
            Vec3 xO {NaN, NaN, NaN};
            Vec3 xI {NaN, NaN, NaN};

            Real l = NaN;

            std::vector<LineSegment> lines;

            Vector pathError {};
            Matrix pathMatrix {};
            Matrix pathErrorJacobian {};
            Vector pathCorrections {};

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
        {   m_Subsystem.invalidateSubsystemTopologyCache(); }

        const PosInfo& getPosInfo(const State& state) const;

        void calcInitZeroLengthGeodesic(State& state, std::function<Vec3(size_t)> GetInitPointGuess) const;

        WrapObstacle* findPrevActiveObstacle(const State& s, size_t obsIdx);
        WrapObstacle* findNextActiveObstacle(const State& s, size_t obsIdx);

        void callCurrentWithPrevAndNext(
                const State& s,
                std::function<void(size_t prevActiveIdx, size_t currentIdx, size_t nextActiveIdx)> f);

        void callCurrentWithPrevAndNext(
                const State& s,
                std::function<void(Vec3 x_O, const WrapObstacle& o, Vec3 x_I)> f) const;

    private:
        PosInfo& updPosInfo(const State& s) const;
        void calcPosInfo(const State& s, PosInfo& posInfo) const;

        WrappingPathSubsystem m_Subsystem;

        MobilizedBody m_OriginBody;
        Vec3 m_OriginPoint;

        MobilizedBody m_TerminationBody;
        Vec3 m_TerminationPoint;

        std::vector<WrapObstacle> m_Obstacles {};

        Real m_PathErrorBound = 0.1;
        size_t m_MaxIter = 10;

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
