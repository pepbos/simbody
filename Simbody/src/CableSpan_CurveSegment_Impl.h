#ifndef SimTK_SIMBODY_CABLE_SPAN_CURVE_SEGMENT_IMPL_H_
#define SimTK_SIMBODY_CABLE_SPAN_CURVE_SEGMENT_IMPL_H_

#include "simbody/internal/CableSpan.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

namespace SimTK
{

//==============================================================================
//                          CURVE SEGMENT IMPL
//==============================================================================
class CurveSegment::Impl final
{
public:
    // The frenet frame along the geodesic on the surface obstacle.
    using FrenetFrame = Transform;

    // Helper struct for logging the frames obtained during shooting of the
    // geodesic.
    struct LocalGeodesicSample final
    {
        LocalGeodesicSample(Real l, FrenetFrame X) :
            length(l), frame(std::move(X))
        {}

        Real length;
        FrenetFrame frame;
    };

    // Discrete auto update cache entry to hold the local (body fixed) geodesic
    // information.
    //
    // At the heart of computing the cable path lies a Newton iteration, which
    // greatly benefits from a good solution estimate. Using a discrete cache
    // entry allows for starting the solver from the previous solution.
    // Furthermore, storing the previous geodesic allows winding over obstacles
    // multiple times.
    // An auto update cache variable is used because the cable path over an
    // obstacle can be seen as a nonholonomic constraint. If we wish to redo
    // the cable path computation for some reason, we must be able to reset to
    // the previous geodesic.
    //
    // Following Scholz2015 we use the P and Q subscript to indicate the
    // curve's start and end contact point respectively. The S subscript is
    // used to indicate the surface frame. For example: point_SP would indicate
    // the initial contact point represented and measured from the surface's
    // origin frame.
    struct InstanceEntry final
    {
        bool isInContactWithSurface() const
        {
            return status == WrappingStatus::InContactWithSurface ||
                   status == WrappingStatus::InitialGuess;
        }

        // Frenet frame at the initial contact point w.r.t. the surface frame.
        FrenetFrame X_SP{};
        // Frenet frame at the final contact point w.r.t. the surface frame.
        FrenetFrame X_SQ{};

        // Length of this curve segment.
        Real length = NaN;

        // Curvatures at the contact points, with the tangential and binormal
        // direction as the first and last element respectively.
        Vec2 curvatures_P{NaN}; // initial contact point.
        Vec2 curvatures_Q{NaN}; // final contact point.

        // Geodesic torsion at the initial and final contact points.
        Real torsion_P = NaN;
        Real torsion_Q = NaN;

        // Jacobi scalars at the final contact point. Contains the
        // translational and rotational direction as the first and last
        // element (see Scholz2015).
        Vec2 jacobi_Q{NaN};
        Vec2 jacobiDot_Q{NaN};

        // Samples recorded from shooting the geodesic over the surface.
        // For an analytic contact geometry this container will be empty.
        // Otherwise, the first and last sample will contain X_SP and X_SQ
        // defined above.
        std::vector<LocalGeodesicSample> samples;

        // The initial integrator stepsize to try next time when shooting a
        // geodesic. This step size estimate will improve each time after
        // shooting a new geodesic.
        Real integratorInitialStepSize = NaN;

        // Given the line spanned by the point before and after this curve
        // segment, the tracking point lies on that line and is nearest to the
        // surface. This point is used to find the touchdown location when the
        // curve is not in contact with the obstacles's surface.
        Vec3 trackingPointOnLine_S{NaN, NaN, NaN};

        // This is a flag to indicate whether the cable comes into contact with
        // the obstacle at all. If the cable is not in contact with the
        // surface, only trackingPointOnLine_S will contain valid data. If the
        // cable is in contact with the surface, trackingPointOnLine_S will not
        // contain valid data.
        WrappingStatus status = WrappingStatus::InitialGuess;
    };

    // Stage::Position level cache entry for holding the curve information in
    // ground frame.
    struct PosEntry final
    {
        // Position and orientation of contact geometry w.r.t. ground.
        Transform X_GS{};

        // Frenet frame at the initial contact point w.r.t. ground.
        FrenetFrame X_GP{};
        // Frenet frame at the final contact point w.r.t. ground.
        FrenetFrame X_GQ{};
    };

//------------------------------------------------------------------------------

    // TODO refactor rule of five and constructor
    Impl()                       = delete;
    Impl(const Impl&)            = delete;
    Impl& operator=(const Impl&) = delete;

    ~Impl()                          = default;
    Impl(Impl&&) noexcept            = default;
    Impl& operator=(Impl&&) noexcept = default;

    Impl(
        CableSpan path,
        MobilizedBodyIndex body,
        const Transform& X_BS,
        ContactGeometry geometry,
        Vec3 initPointGuess);

//------------------------------------------------------------------------------

    // Allocate state variables and cache entries.
    void realizeTopology(State& s)
    {
        // Allocate an auto-update discrete variable for the last computed
        // geodesic.
        Value<InstanceEntry>* cache = new Value<InstanceEntry>();
        m_InstanceIx = updSubsystem().allocateAutoUpdateDiscreteVariable(
            s,
            Stage::Report,
            cache,
            Stage::Position);

        PosEntry posInfo{};
        m_PosIx = updSubsystem().allocateCacheEntry(
            s,
            Stage::Position,
            Stage::Infinity,
            new Value<PosEntry>(posInfo));
    }

    void realizePosition(const State& s, Vec3 prevPoint_G, Vec3 nextPoint_G)
        const;

    // Realize the CableSpan, that this segment is a part of, to
    // Stage::Position.
    void realizeCablePosition(const State& s) const;

    void invalidatePosEntry(const State& state) const
    {
        getSubsystem().markCacheValueNotRealized(state, m_PosIx);
    }

//------------------------------------------------------------------------------

    const CableSpan& getCable() const
    {
        return m_Path;
    }

    CurveSegmentIndex getIndex() const
    {
        return m_Index;
    }

    void setIndex(CurveSegmentIndex ix)
    {
        m_Index = ix;
    }

    // Set the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    void setContactPointHint(Vec3 contactPointHint_S)
    {
        m_ContactPointHint_S = contactPointHint_S;
    }

    // Get the user defined point that controls the initial wrapping path.
    // Point is in surface coordinates.
    Vec3 getContactPointHint() const
    {
        return m_ContactPointHint_S;
    }

    const ContactGeometry& getContactGeometry() const
    {
        return m_Geometry;
    }

    const DecorativeGeometry& getDecoration() const
    {
        return m_Decoration;
    }

    const MobilizedBody& getMobilizedBody() const;

    const Transform& getXformSurfaceToBody() const
    {
        return m_X_BS;
    }

    void setXformSurfaceToBody(Transform X_BS)
    {
        m_X_BS = std::move(X_BS);
    }

    const CableSubsystem& getSubsystem() const
    {
        return *m_Subsystem;
    }

    CableSubsystem& updSubsystem()
    {
        return *m_Subsystem;
    }

//------------------------------------------------------------------------------

    const InstanceEntry& getInstanceEntry(const State& s) const
    {
        const CableSubsystem& subsystem = getSubsystem();
        if (!subsystem.isDiscreteVarUpdateValueRealized(s, m_InstanceIx)) {
            updInstanceEntry(s) = getPrevInstanceEntry(s);
            subsystem.markDiscreteVarUpdateValueRealized(s, m_InstanceIx);
        }
        return Value<InstanceEntry>::downcast(
            subsystem.getDiscreteVarUpdateValue(s, m_InstanceIx));
    }

    // TODO rename posInfo to PosEntry
    const PosEntry& getPosInfo(const State& s) const
    {
        return Value<PosEntry>::downcast(
            getSubsystem().getCacheEntry(s, m_PosIx));
    }

//------------------------------------------------------------------------------

    Transform calcSurfaceFrameInGround(const State& s) const
    {
        return getMobilizedBody().getBodyTransform(s).compose(m_X_BS);
    }

    bool isInContactWithSurface(const State& s) const
    {
        return getInstanceEntry(s).isInContactWithSurface();
    }

    // See CurveSegment::calcPathPoints for description.
    int calcPathPoints(
        const State& state,
        std::function<void(Vec3 point_G)>& sink,
        int nSamples = 0) const;

    int calcPathPointsAndTangents(
        const State& state,
        std::function<void(Vec3 point_G, UnitVec3 tangent_G)>& sink,
        int nSamples = 0) const;

    // Compute a new geodesic from provided initial conditions.
    // This method will update the InstanceEntry cache, and invalidates the
    // PosEntry cache.
    void calcLocalGeodesic(
        const State& s,
        Vec3 point_S,
        Vec3 tangent_S,
        Real length,
        Real stepSizeHint) const;

    // Lift curve from surface, and start tracking the given point.
    // This method will update the InstanceEntry cache, and invalidates the
    // PosEntry cache.
    void liftCurveFromSurface(const State& s, Vec3 trackingPoint_S) const;

//------------------------------------------------------------------------------

private:
    InstanceEntry& updInstanceEntry(const State& state) const
    {
        return Value<InstanceEntry>::updDowncast(
            getSubsystem().updDiscreteVarUpdateValue(state, m_InstanceIx));
    }

    const InstanceEntry& getPrevInstanceEntry(const State& state) const
    {
        return Value<InstanceEntry>::downcast(
            getSubsystem().getDiscreteVariable(state, m_InstanceIx));
    }

    InstanceEntry& updPrevInstanceEntry(State& state) const
    {
        return Value<InstanceEntry>::updDowncast(
            getSubsystem().updDiscreteVariable(state, m_InstanceIx));
    }

    PosEntry& updPosInfo(const State& state) const
    {
        return Value<PosEntry>::updDowncast(
            getSubsystem().updCacheEntry(state, m_PosIx));
    }

//------------------------------------------------------------------------------

    // Subsystem info.
    CableSubsystem* m_Subsystem; // The subsystem this segment belongs to.
    CableSpan m_Path;            // The path this segment belongs to.
    CurveSegmentIndex m_Index;   // The index in its path.

    // MobilizedBody that surface is attached to.
    MobilizedBodyIndex m_Body;
    // Surface to body transform.
    Transform m_X_BS;

    // Obstacle surface.
    ContactGeometry m_Geometry;
    DecorativeGeometry m_Decoration;

    // Topology cache.
    CacheEntryIndex m_PosIx;
    DiscreteVariableIndex m_InstanceIx;

    // Initial contact point hint used to setup the initial path.
    Vec3 m_ContactPointHint_S{NaN, NaN, NaN};

    // Helper class for unit tests.
    friend CableSubsystemTestHelper;
};

} // namespace SimTK

#endif
