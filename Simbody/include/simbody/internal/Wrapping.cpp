#include "Wrapping.h"

#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKmath.h"
#include "WrappingImpl.h"
#include "simmath/internal/ContactGeometry.h"
#include <cstddef>
#include <stdexcept>

using namespace SimTK;

using Correction        = ContactGeometry::GeodesicCorrection;
using FrameVariation    = ContactGeometry::GeodesicFrameVariation;
using FrenetFrame       = ContactGeometry::FrenetFrame;
using GeodesicInfo      = CurveSegment::Impl::PosInfo;
using LocalGeodesicInfo = CurveSegment::Impl::LocalGeodesic::LocalGeodesicInfo;
using GeodesicInitialConditions =
    CurveSegment::Impl::LocalGeodesic::GeodesicInitialConditions;
using GeodesicJacobian = Vec4;
using PointVariation   = ContactGeometry::GeodesicPointVariation;
using Variation        = ContactGeometry::GeodesicVariation;
using LineSegment      = CableSpan::LineSegment;
using Status           = CurveSegment::Impl::Status;

//==============================================================================
//                                CONSTANTS
//==============================================================================
namespace
{
static const CoordinateAxis TangentAxis  = ContactGeometry::TangentAxis;
static const CoordinateAxis NormalAxis   = ContactGeometry::NormalAxis;
static const CoordinateAxis BinormalAxis = ContactGeometry::BinormalAxis;
static const int GeodesicDOF             = 4;
} // namespace

//==============================================================================
//                                SOLVER
//==============================================================================

namespace
{
class SolverDataCache;
SolverDataCache& findDataCache(size_t nActive);

struct SolverData
{
    std::vector<LineSegment> lineSegments;

    Matrix pathErrorJacobian;
    Vector pathCorrection;
    Vector pathError;
    Matrix mat;
    // TODO Cholesky decomposition...
    FactorLU matInv;
    Vector vec;
};

class SolverDataCache
{
private:
    SolverDataCache(size_t n)
    {
        static constexpr int Q = 4;
        static constexpr int C = 4;

        m_Data.lineSegments.resize(n + 1);
        m_Data.pathErrorJacobian = Matrix(C * n, Q * n, 0.);
        m_Data.pathCorrection    = Vector(Q * n, 0.);
        m_Data.pathError         = Vector(C * n, 0.);
        m_Data.mat               = Matrix(Q * n, Q * n, NaN);
        m_Data.vec               = Vector(Q * n, NaN);
    }

public:
    void lock()
    {
        m_Mutex.lock();
    }

    void unlock()
    {
        m_Mutex.unlock();
    }

    SolverData& updData()
    {
        return m_Data;
    }

private:
    SolverData m_Data;
    std::mutex m_Mutex{};

    friend SolverDataCache& findDataCache(size_t nActive);
};

SolverDataCache& findDataCache(size_t nActive)
{
    static std::vector<SolverDataCache> s_GlobalCache{};
    static std::mutex s_GlobalLock{};

    {
        std::lock_guard<std::mutex> lock(s_GlobalLock);

        for (size_t i = s_GlobalCache.size(); i < nActive - 1; ++i) {
            int n = i + 1;
            s_GlobalCache.emplace_back(nActive);
        }
    }

    return s_GlobalCache.at(nActive - 1);
}

const Correction* calcPathCorrections(SolverData& data)
{
    Real w = data.pathError.normInf();

    data.mat = data.pathErrorJacobian.transpose() * data.pathErrorJacobian;
    for (int i = 0; i < data.mat.nrow(); ++i) {
        data.mat[i][i] += w;
    }
    data.matInv = data.mat;
    data.vec    = data.pathErrorJacobian.transpose() * data.pathError;
    data.matInv.solve(data.vec, data.pathCorrection);

    static_assert(
        sizeof(Correction) == sizeof(double) * GeodesicDOF,
        "Invalid size of corrections vector");
    SimTK_ASSERT(
        data.pathCorrection.size() * sizeof(double) == n * sizeof(Correction),
        "Invalid size of path corrections vector");
    return reinterpret_cast<const Correction*>(&data.pathCorrection[0]);
}
} // namespace

//==============================================================================
//                          GEODESIC INITIAL CONDITIONS
//==============================================================================

GeodesicInitialConditions GeodesicInitialConditions::CreateCorrected(
    const FrenetFrame& KP,
    const Variation& dKP,
    Real l,
    const Correction& c)
{
    GeodesicInitialConditions g0;

    Vec3 v = dKP[1] * c;
    g0.x   = KP.p() + v;

    Vec3 w           = dKP[0] * c;
    const UnitVec3 t = KP.R().getAxisUnitVec(ContactGeometry::TangentAxis);
    g0.t             = t + cross(w, t);

    Real dl = c[3];
    g0.l    = l + dl;

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::CreateFromGroundInSurfaceFrame(
    const Transform& X_GS,
    Vec3 x_G,
    Vec3 t_G,
    Real l)
{
    GeodesicInitialConditions g0;

    g0.x = X_GS.shiftBaseStationToFrame(x_G);
    g0.t = X_GS.xformBaseVecToFrame(t_G);
    g0.l = l;

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::CreateZeroLengthGuess(
    Vec3 prev_QS,
    Vec3 xGuess_S)
{
    GeodesicInitialConditions g0;

    g0.x = xGuess_S;
    g0.t = xGuess_S - prev_QS;
    g0.l = 0.;

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::CreateAtTouchdown(
    Vec3 prev_QS,
    Vec3 next_PS,
    Vec3 trackingPointOnLine)
{
    GeodesicInitialConditions g0;

    g0.x = trackingPointOnLine;
    g0.t = next_PS - prev_QS;
    g0.l = 0.;

    return g0;
}

//==============================================================================
//                                SURFACE IMPL
//==============================================================================
using LocalGeodesic = CurveSegment::Impl::LocalGeodesic;

//------------------------------------------------------------------------------
//                               REALIZE CACHE
//------------------------------------------------------------------------------
void LocalGeodesic::realizeTopology(State& s)
{
    // Allocate an auto-update discrete variable for the last computed geodesic.
    CacheEntry cache{};
    m_CacheIx = m_Subsystem.allocateAutoUpdateDiscreteVariable(
        s,
        Stage::Velocity,
        new Value<CacheEntry>(cache),
        Stage::Position);
}

void LocalGeodesic::realizePosition(const State& s) const
{
    if (!m_Subsystem.isDiscreteVarUpdateValueRealized(s, m_CacheIx)) {
        updCacheEntry(s) = getPrevCacheEntry(s);
        m_Subsystem.markDiscreteVarUpdateValueRealized(s, m_CacheIx);
    }
}

//------------------------------------------------------------------------------
//                               PUBLIC METHODS
//------------------------------------------------------------------------------
const LocalGeodesicInfo& LocalGeodesic::calcInitialGeodesic(
    State& s,
    const GeodesicInitialConditions& g0) const
{
    // TODO is this correct?
    CacheEntry& cache = updPrevCacheEntry(s);
    calcGeodesic(g0, cache);
    updCacheEntry(s) = cache;
    m_Subsystem.markDiscreteVarUpdateValueRealized(s, m_CacheIx);
}

const LocalGeodesicInfo& LocalGeodesic::calcLocalGeodesicInfo(
    const State& s,
    Vec3 prev_QS,
    Vec3 next_PS) const
{
    realizePosition(s);
    CacheEntry& cache = updCacheEntry(s);
    calcCacheEntry(prev_QS, next_PS, cache);
    return cache;
}

void LocalGeodesic::applyGeodesicCorrection(const State& s, const Correction& c) const
{
    realizePosition(s);

    // Get the previous geodesic.
    const CacheEntry& g = getCacheEntry(s);

    // Get corrected initial conditions.
    const GeodesicInitialConditions g0 =
        GeodesicInitialConditions::CreateCorrected(g.KP, g.dKP, g.length, c);

    // Shoot the new geodesic.
    calcGeodesic(g0, updCacheEntry(s));
}

size_t LocalGeodesic::calcPathPoints(const State& s, std::vector<Vec3>& points) const
{
    realizePosition(s);

    size_t count            = 0;
    const CacheEntry& cache = getCacheEntry(s);
    for (Vec3 p : cache.points) {
        points.push_back(p);
        ++count;
    }
    throw std::runtime_error("NOTYETIMPLEMENTED for analytic");
    return count;
}

//------------------------------------------------------------------------------
//                               PRIVATE METHODS
//------------------------------------------------------------------------------
void LocalGeodesic::calcGeodesic(
    const GeodesicInitialConditions& g0,
    CacheEntry& cache) const
{
    // Compute geodesic start boundary frame and variation.
    m_Geometry.calcNearestFrenetFrameFast(g0.x, g0.t, cache.KP);
    m_Geometry.calcGeodesicStartFrameVariation(cache.KP, cache.dKP);

    // Compute geodesic end boundary frame amd variation (shoot new geodesic).
    m_Geometry.calcGeodesicEndFrameVariationImplicitly(
        cache.KP.p(),
        cache.KP.R().getAxisUnitVec(ContactGeometry::TangentAxis),
        g0.l,
        cache.sHint,
        cache.KQ,
        cache.dKQ,
        cache.points);

    // TODO  update step size.
    // TODO  update line tracking?
    throw std::runtime_error("NOTYETIMPLEMENTED");
}

void LocalGeodesic::calcLiftoffIfNeeded(
    const Vec3& prev_QS,
    const Vec3& next_PS,
    CacheEntry& cache) const
{
    // Only attempt liftoff when currently wrapping the surface.
    LocalGeodesicInfo& g = cache;
    if (g.status != Status::Ok) {
        return;
    }

    // The curve length must have shrunk completely before lifting off.
    if (g.length < 0.) {
        return;
    }

    // For a zero-length curve, trigger liftoff when the prev and next points
    // lie above the surface plane.
    if (
        dot(prev_QS - g.KP.p(), g.KP.R().getAxisUnitVec(NormalAxis)) <= 0.
        &&
        dot(next_PS - g.KP.p(), g.KP.R().getAxisUnitVec(NormalAxis)) <= 0.)
    {
        // No liftoff.
        return;
    }

    // Liftoff detected: update status.
    g.status = Status::Liftoff;
    // Initialize the tracking point from the last geodesic start point.
    cache.trackingPointOnLine = g.KP.p();

}

void LocalGeodesic::calcTouchdownIfNeeded(
    const Vec3& prev_QS,
    const Vec3& next_PS,
    CacheEntry& cache) const
{
    // Only attempt touchdown when liftoff.
    LocalGeodesicInfo& g = cache;
    if (g.status != Status::Liftoff) {
        return;
    }

    // Detect touchdown by computing the point on the line from x_QS to x_PS
    // that is nearest to the surface.
    if(!m_Geometry.calcNearestPointOnLine(
        prev_QS,
        next_PS,
        cache.trackingPointOnLine,
        m_TouchdownIter,
        m_TouchdownAccuracy))
    {
        // No touchdown detected.
        return;
    }

    // Touchdown detected: Remove the liftoff status flag.
    g.status = Status::Ok;
    // Shoot a zero length geodesic at the touchdown point.
    GeodesicInitialConditions g0 =
        GeodesicInitialConditions::CreateAtTouchdown(
                prev_QS,
                next_PS,
                cache.trackingPointOnLine);
    calcGeodesic(g0, cache);
}

void LocalGeodesic::assertSurfaceBounds(
    const Vec3& prev_QS,
    const Vec3& next_PS) const
{
    // Make sure that the previous point does not lie inside the surface.
    SimTK_ASSERT(m_Geometry.calcSurfaceValue(prev_QS) < 0.,
        "Unable to wrap over surface: Preceding point lies inside the surface");
    SimTK_ASSERT(m_Geometry.calcSurfaceValue(next_PS) < 0.,
        "Unable to wrap over surface: Next point lies inside the surface");
}

void LocalGeodesic::calcCacheEntry(
    const Vec3& prev_QS,
    const Vec3& next_PS,
    CacheEntry& cache) const
{
    LocalGeodesicInfo& g = cache;

    if (g.status == Status::Disabled) {
        return;
    }

    assertSurfaceBounds(prev_QS, next_PS);
    calcTouchdownIfNeeded(prev_QS, next_PS, cache);
    calcLiftoffIfNeeded(prev_QS, next_PS, cache);
}

//==============================================================================
//                               OBSTACLE
//==============================================================================

CurveSegment::Impl::Impl(
    CableSpan path,
    CurveSegmentIndex ix,
    const MobilizedBody& mobod,
    const Transform& X_BS,
    ContactGeometry geometry,
    Vec3 initPointGuess) :
    m_Subsystem(path.getImpl().getSubsystem()),
    m_Path(path), m_Index(ix), m_Mobod(mobod), m_Offset(X_BS),
    m_Surface(m_Subsystem, geometry, initPointGuess)
{}

void CurveSegment::Impl::realizeTopology(State& s)
{
    // Allocate position level cache.
    PosInfo posInfo{};
    m_PosInfoIx = m_Subsystem.allocateCacheEntry(
        s,
        Stage::Position,
        new Value<PosInfo>(posInfo));
}

void CurveSegment::Impl::realizePosition(const State& s) const
{
    if (!m_Subsystem.isCacheValueRealized(s, m_PosInfoIx)) {
        calcPosInfo(s, updPosInfo(s));
        m_Subsystem.markCacheValueRealized(s, m_PosInfoIx);
    }
}

void CurveSegment::Impl::calcInitZeroLengthGeodesic(State& s, Vec3 prev_QG)
    const
{
    // Get tramsform from local surface frame to ground.
    Transform X_GS = m_Mobod.getBodyTransform(s).compose(m_Offset);

    Vec3 prev_QS = X_GS.shiftBaseStationToFrame(prev_QG);
    Vec3 xGuess_S =
        m_Surface.getInitialPointGuess(); // TODO move into function call?

    GeodesicInitialConditions g0 =
        GeodesicInitialConditions::CreateZeroLengthGuess(
            prev_QS,
            xGuess_S);
    m_Surface.calcInitialGeodesic(s, g0);

    m_Subsystem.markCacheValueNotRealized(s, m_PosInfoIx);
}

void CurveSegment::Impl::applyGeodesicCorrection(
    const State& s,
    const CurveSegment::Impl::Correction& c) const
{
    // Apply correction to curve.
    m_Surface.applyGeodesicCorrection(s, c);

    // Invalidate position level cache.
    m_Subsystem.markCacheValueNotRealized(s, m_PosInfoIx);
}

size_t CurveSegment::Impl::calcPathPoints(
    const State& s,
    std::vector<Vec3>& points) const
{
    const Transform& X_GS = getPosInfo(s).X_GS;
    size_t n              = m_Surface.calcPathPoints(s, points);
    for (size_t i = points.size() - n; i < points.size(); ++i) {
        points.at(i) = X_GS.shiftFrameStationToBase(points.at(i));
    }
    return n;
}
void CurveSegment::Impl::calcGeodesicInGround(
        const LocalGeodesicInfo& geodesic_S,
        const Transform& X_GS,
        PosInfo& posInfo) const
{
    posInfo.X_GS = X_GS;

    // Store the the local geodesic in ground frame.
    posInfo.KP = X_GS.compose(geodesic_S.KP);
    posInfo.KQ = X_GS.compose(geodesic_S.KQ);

    posInfo.dKP[0] = X_GS.R() * geodesic_S.dKP[0];
    posInfo.dKP[1] = X_GS.R() * geodesic_S.dKP[1];

    posInfo.dKQ[0] = X_GS.R() * geodesic_S.dKQ[0];
    posInfo.dKQ[1] = X_GS.R() * geodesic_S.dKQ[1];

    // TODO use SpatialVec for variation.
    /* for (size_t i = 0; i < GeodesicDOF; ++i) { */
    /*     posInfo.dKP[i][0] = X_GS.xformFrameVecToBase(geodesic_S.dKP[i][0]) */
    /*     posInfo.dKP[i][1] = X_GS.xformFrameVecToBase(geodesic_S.dKP[i][1]) */

    /*     posInfo.dKQ[i][0] = X_GS.xformFrameVecToBase(geodesic_S.dKQ[i][0]) */
    /*     posInfo.dKQ[i][1] = X_GS.xformFrameVecToBase(geodesic_S.dKQ[i][1]) */
    /* } */

    posInfo.length = geodesic_S.length;

    throw std::runtime_error("NOTYETIMPLEMENTED: Check transformation order");
}

void CurveSegment::Impl::calcPosInfo(const State& s, PosInfo& posInfo) const
{
    if (m_Surface.getStatus(s) == Status::Disabled) {
        return;
    }

    // Compute tramsform from local surface frame to ground.
    const Transform& X_GS = posInfo.X_GS =
        m_Mobod.getBodyTransform(s).compose(m_Offset);

    // Get the path points before and after this segment.
    const Vec3 prev_G = m_Path.getImpl().findPrevPoint(s, m_Index);
    const Vec3 next_G = m_Path.getImpl().findNextPoint(s, m_Index);

    // Transform the prev and next path points to the surface frame.
    const Vec3 prev_S = X_GS.shiftBaseStationToFrame(prev_G);
    const Vec3 next_S = X_GS.shiftBaseStationToFrame(next_G);

    // Detect liftoff, touchdown and potential invalid configurations.
    // TODO this doesnt follow the regular invalidation scheme...
    // Grab the last geodesic that was computed.
    const LocalGeodesicInfo& geodesic_S =
        m_Surface.calcLocalGeodesicInfo(s, prev_S, next_S);

    // Store the the local geodesic in ground frame.
    calcGeodesicInGround(geodesic_S, X_GS, updPosInfo(s));
}

//==============================================================================
//                                PATH HELPERS
//==============================================================================

namespace
{

static const int N_PATH_CONSTRAINTS = 4;

// TODO this is awkward
void addBlock(const Vec4& values, Matrix& block)
{
    for (int i = 0; i < Vec4::size(); ++i) {
        block[0][i] = values[i];
    }
}

void addDirectionJacobian(
    const LineSegment& e,
    const UnitVec3& axis,
    const PointVariation& dx,
    Matrix& J,
    bool invert = false)
{
    Vec3 y = axis - e.d * dot(e.d, axis);
    y /= e.l * (invert ? 1. : -1);
    addBlock(~dx * y, J);
}

double calcPathError(
    const LineSegment& e,
    const Rotation& R,
    CoordinateAxis axis)
{
    return dot(e.d, R.getAxisUnitVec(axis));
}

GeodesicJacobian addPathErrorJacobian(
    const LineSegment& e,
    const UnitVec3& axis,
    const Variation& dK,
    Matrix& J,
    bool invertV = false)
{
    addDirectionJacobian(e, axis, dK[1], J, invertV);
    addBlock(~dK[0] * cross(axis, e.d), J);
}

} // namespace

//==============================================================================
//                                PATH IMPL
//==============================================================================

void CableSpan::Impl::realizeTopology(State& s)
{
    PosInfo posInfo{};
    m_PosInfoIx = m_Subsystem.allocateCacheEntry(
        s,
        Stage::Position,
        new Value<PosInfo>(posInfo));

    VelInfo velInfo{};
    m_VelInfoIx = m_Subsystem.allocateCacheEntry(
        s,
        Stage::Velocity,
        new Value<VelInfo>(velInfo));
}

void CableSpan::Impl::realizePosition(const State& s) const
{
    if (m_Subsystem.isCacheValueRealized(s, m_PosInfoIx)) {
        return;
    }
    calcPosInfo(s, updPosInfo(s));
    m_Subsystem.markCacheValueRealized(s, m_PosInfoIx);
}

void CableSpan::Impl::realizeVelocity(const State& s) const
{
    if (m_Subsystem.isCacheValueRealized(s, m_VelInfoIx)) {
        return;
    }
    calcVelInfo(s, updVelInfo(s));
    m_Subsystem.markCacheValueRealized(s, m_VelInfoIx);
}

const CableSpan::Impl::PosInfo& CableSpan::Impl::getPosInfo(
    const State& s) const
{
    realizePosition(s);
    return Value<PosInfo>::downcast(m_Subsystem.getCacheEntry(s, m_PosInfoIx));
}

CableSpan::Impl::PosInfo& CableSpan::Impl::updPosInfo(const State& s) const
{
    return Value<PosInfo>::updDowncast(
        m_Subsystem.updCacheEntry(s, m_PosInfoIx));
}

const CableSpan::Impl::VelInfo& CableSpan::Impl::getVelInfo(
    const State& s) const
{
    realizeVelocity(s);
    return Value<VelInfo>::downcast(m_Subsystem.getCacheEntry(s, m_VelInfoIx));
}

CableSpan::Impl::VelInfo& CableSpan::Impl::updVelInfo(const State& s) const
{
    return Value<VelInfo>::updDowncast(
        m_Subsystem.updCacheEntry(s, m_VelInfoIx));
}

void CableSpan::Impl::calcInitZeroLengthGeodesic(State& s) const
{
    Vec3 prev_QG =
        m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
    size_t i = 0;
    for (const CurveSegment& obstacle : m_CurveSegments) {
        if (obstacle.getImpl().getStatus(s) == Status::Disabled) {
            continue;
        }
        obstacle.getImpl().calcInitZeroLengthGeodesic(s, prev_QG);

        prev_QG = obstacle.getImpl().getPosInfo(s).KQ.p();
    }
}

const CurveSegment* CableSpan::Impl::findPrevActiveCurveSegment(
    const State& s,
    CurveSegmentIndex ix) const
{
    for (int i = ix - 1; i > 0; --i) {
        // Find the active segment before the current.
        if (m_CurveSegments.at(CurveSegmentIndex(i)).getImpl().isActive(s)) {
            return &m_CurveSegments.at(CurveSegmentIndex(i));
        }
    }
    return nullptr;
}

const CurveSegment* CableSpan::Impl::findNextActiveCurveSegment(
    const State& s,
    CurveSegmentIndex ix) const
{
    // Find the active segment after the current.
    for (int i = ix + 1; i < m_CurveSegments.size(); ++i) {
        if (m_CurveSegments.at(CurveSegmentIndex(i)).getImpl().isActive(s)) {
            return &m_CurveSegments.at(CurveSegmentIndex(i));
        }
    }
    return nullptr;
}

Vec3 CableSpan::Impl::findPrevPoint(const State& s, CurveSegmentIndex ix) const
{
    const CurveSegment* segment = findPrevActiveCurveSegment(s, ix);
    return segment ? segment->getImpl().getPosInfo(s).KQ.p()
                   : m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(
                         m_OriginPoint);
}

Vec3 CableSpan::Impl::findNextPoint(const State& s, CurveSegmentIndex ix) const
{
    const CurveSegment* segment = findNextActiveCurveSegment(s, ix);
    return segment
               ? segment->getImpl().getPosInfo(s).KP.p()
               : m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(
                     m_TerminationPoint);
}

template <size_t N>
void CableSpan::Impl::calcPathErrorVector(
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Vector& pathError) const
{
    size_t lineIx = 0;
    ptrdiff_t row = -1;

    for (int i = 0; i < getNumCurveSegments(); ++i) {
        const CurveSegment& segment = getCurveSegment(CurveSegmentIndex(i));
        if (!segment.getImpl().isActive(s)) {
            continue;
        }

        const GeodesicInfo& g = segment.getImpl().getPosInfo(s);
        for (CoordinateAxis axis : axes) {
            pathError(++row) = calcPathError(lines.at(lineIx), g.KP.R(), axis);
        }
        ++lineIx;
        for (CoordinateAxis axis : axes) {
            pathError(++row) = calcPathError(lines.at(lineIx), g.KQ.R(), axis);
        }
    }
}

template <size_t N>
void CableSpan::Impl::calcPathErrorJacobian(
    const State& s,
    const std::vector<LineSegment>& lines,
    std::array<CoordinateAxis, N> axes,
    Matrix& J) const
{
    constexpr size_t Nq = GeodesicDOF;

    // TODO perhaps just not make method static.
    const size_t n = lines.size() - 1;

    SimTK_ASSERT(
        J.rows() == n * N,
        "Invalid number of rows in jacobian matrix");
    SimTK_ASSERT(
        J.cols() == n * Nq,
        "Invalid number of columns in jacobian matrix");

    size_t row      = 0;
    size_t col      = 0;
    size_t activeIx = 0;
    for (const CurveSegment& segment : m_CurveSegments) {
        if (!segment.getImpl().isActive(s)) {
            continue;
        }
        const GeodesicInfo& g = segment.getImpl().getPosInfo(s);

        const LineSegment& l_P = lines.at(activeIx);
        const LineSegment& l_Q = lines.at(activeIx + 1);

        const CurveSegmentIndex ix = segment.getImpl().getIndex();
        const CurveSegment* prev   = findPrevActiveCurveSegment(s, ix);
        const CurveSegment* next   = findNextActiveCurveSegment(s, ix);

        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_P    = g.KP.R().getAxisUnitVec(axis);
            const Variation& dK_P = g.dKP;

            addPathErrorJacobian(l_P, a_P, dK_P, J.block(row, col, 1, Nq));

            if (prev) {
                const Variation& prev_dK_Q = prev->getImpl().getPosInfo(s).dKQ;
                addDirectionJacobian(
                    l_P,
                    a_P,
                    prev_dK_Q[1],
                    J.block(row, col - Nq, 1, Nq),
                    true);
            }
            ++row;
        }

        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_Q    = g.KQ.R().getAxisUnitVec(axis);
            const Variation& dK_Q = g.dKQ;

            addPathErrorJacobian(
                l_Q,
                a_Q,
                dK_Q,
                J.block(row, col, 1, Nq),
                true);

            if (next) {
                const Variation& next_dK_P = next->getImpl().getPosInfo(s).dKP;
                addDirectionJacobian(
                    l_Q,
                    a_Q,
                    next_dK_P[1],
                    J.block(row, col + Nq, 1, Nq));
            }
            ++row;
        }

        col += Nq;
    };
}

double CableSpan::Impl::calcPathLength(
    const State& s,
    const std::vector<LineSegment>& lines) const
{
    double lTot = 0.;
    for (const LineSegment& line : lines) {
        // TODO spell out as length.
        lTot += line.l;
    }

    for (const CurveSegment& segment : m_CurveSegments) {
        if (!segment.getImpl().isActive(s)) {
            continue;
        }
        lTot += segment.getImpl().getPosInfo(s).length;
    }
    return lTot;
}

void CableSpan::Impl::calcLineSegments(
    const State& s,
    Vec3 p_O,
    Vec3 p_I,
    std::vector<LineSegment>& lines) const
{
    lines.resize(m_CurveSegments.size() + 1);
    lines.clear();

    Vec3 lineStart = std::move(p_O);
    for (const CurveSegment& segment : m_CurveSegments) {
        if (!segment.getImpl().isActive(s)) {
            continue;
        }

        const GeodesicInfo& g = segment.getImpl().getPosInfo(s);
        const Vec3 lineEnd    = g.KP.p();
        lines.emplace_back(lineStart, lineEnd);

        lineStart = g.KQ.p();
    }
    lines.emplace_back(lineStart, p_I);
}

size_t CableSpan::Impl::countActive(const State& s) const
{
    size_t count = 0;
    for (const CurveSegment& segment : m_CurveSegments) {
        if (segment.getImpl().isActive(s)) {
            ++count;
        }
    }
    return count;
}

void CableSpan::Impl::calcPosInfo(const State& s, PosInfo& posInfo) const
{
    // Path origin and termination point.
    const Vec3 x_O =
        m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
    const Vec3 x_I =
        m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(
            m_TerminationPoint);

    const std::array<CoordinateAxis, 2> axes{NormalAxis, BinormalAxis};

    for (posInfo.loopIter = 0; posInfo.loopIter < m_PathMaxIter;
         ++posInfo.loopIter) {
        const size_t nActive = countActive(s);

        // Grab the shared data cache for computing the matrices, and lock it.
        SolverDataCache& dataCache = findDataCache(nActive);
        dataCache.lock();
        SolverData& data = dataCache.updData();

        // Compute the straight-line segments.
        calcLineSegments(s, x_O, x_I, data.lineSegments);

        // Evaluate path error, and stop when converged.
        calcPathErrorVector<2>(s, data.lineSegments, axes, data.pathError);
        const Real maxPathError = data.pathError.normInf();
        if (maxPathError < m_PathErrorBound) {
            return;
        }

        // Evaluate the path error jacobian.
        calcPathErrorJacobian<2>(
            s,
            data.lineSegments,
            axes,
            data.pathErrorJacobian);

        // Compute path corrections.
        const Correction* corrIt = calcPathCorrections(data);

        // Apply path corrections.
        for (const CurveSegment& obstacle : m_CurveSegments) {
            if (!obstacle.getImpl().isActive(s)) {
                continue;
            }
            obstacle.getImpl().applyGeodesicCorrection(s, *corrIt);
            ++corrIt;
        }

        // Release the lock on the shared data.
        dataCache.unlock();

        // Path has changed: invalidate each segment's cache.
        for (const CurveSegment& obstacle : m_CurveSegments) {
            obstacle.getImpl().invalidatePositionLevelCache(s);
        }
    }

    throw std::runtime_error("Failed to converge");
}

void CableSpan::Impl::calcVelInfo(const State& s, VelInfo& velInfo) const
{
    const PosInfo& pos = getPosInfo(s);

    Real lengthDot = 0.;

    Vec3 v_GQ = m_OriginBody.findStationVelocityInGround(s, m_OriginPoint);
    const CurveSegment* lastActive = nullptr;
    for (const CurveSegment& obstacle : m_CurveSegments) {
        if (!obstacle.getImpl().isActive(s)) {
            continue;
        }

        const GeodesicInfo& g = obstacle.getImpl().getPosInfo(s);
        const UnitVec3 e_G    = g.KP.R().getAxisUnitVec(TangentAxis);

        Vec3 next_v_GQ;
        Vec3 v_GP;
        obstacle.getImpl().calcContactPointVelocitiesInGround(
            s,
            v_GP,
            next_v_GQ);

        lengthDot += dot(e_G, v_GP - v_GQ);

        v_GQ       = next_v_GQ;
        lastActive = &obstacle;
    }

    const Vec3 v_GP =
        m_TerminationBody.findStationVelocityInGround(s, m_TerminationPoint);
    const UnitVec3 e_G =
        lastActive ? lastActive->getImpl().getPosInfo(s).KQ.R().getAxisUnitVec(
                         TangentAxis)
                   : UnitVec3(pos.xI - pos.xO);

    lengthDot += dot(e_G, v_GP - v_GQ);
}
void CableSpan::Impl::applyBodyForces(
    const State& state,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    throw std::runtime_error("NOTYETIMPLEMENTED");
}

//==============================================================================
//            SUBSYSTEM
//==============================================================================

bool CableSubsystem::isInstanceOf(const Subsystem& s)
{
    return Impl::isA(s.getSubsystemGuts());
}

const CableSubsystem& CableSubsystem::downcast(const Subsystem& s)
{
    assert(isInstanceOf(s));
    return static_cast<const CableSubsystem&>(s);
}
CableSubsystem& CableSubsystem::updDowncast(Subsystem& s)
{
    assert(isInstanceOf(s));
    return static_cast<CableSubsystem&>(s);
}

const CableSubsystem::Impl& CableSubsystem::getImpl() const
{
    return SimTK_DYNAMIC_CAST_DEBUG<const Impl&>(getSubsystemGuts());
}
CableSubsystem::Impl& CableSubsystem::updImpl()
{
    return SimTK_DYNAMIC_CAST_DEBUG<Impl&>(updSubsystemGuts());
}

// Create Subsystem but don't associate it with any System. This isn't much use
// except for making std::vectors, which require a default constructor to be
// available.
CableSubsystem::CableSubsystem()
{
    adoptSubsystemGuts(new Impl());
}

CableSubsystem::CableSubsystem(MultibodySystem& mbs)
{
    adoptSubsystemGuts(new Impl());
    mbs.adoptSubsystem(*this);
} // steal ownership

int CableSubsystem::getNumPaths() const
{
    return getImpl().getNumPaths();
}

const CableSpan& CableSubsystem::getPath(WrappingPathIndex cableIx) const
{
    return getImpl().getCablePath(cableIx);
}

CableSpan& CableSubsystem::updPath(WrappingPathIndex cableIx)
{
    return updImpl().updCablePath(cableIx);
}
