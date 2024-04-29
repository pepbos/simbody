#include "Wrapping.h"

#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKmath.h"
#include "WrappingImpl.h"
#include "simmath/internal/ContactGeometry.h"
#include <cstddef>
#include <memory>
#include <stdexcept>

using namespace SimTK;

using Correction        = ContactGeometry::GeodesicCorrection;
using FrameVariation    = ContactGeometry::GeodesicFrameVariation;
using FrenetFrame       = ContactGeometry::FrenetFrame;
using GeodesicInfo      = CurveSegment::Impl::PosInfo;
using LocalGeodesicInfo = CurveSegment::Impl::LocalGeodesic::LocalGeodesicInfo;
using LocalGeodesicSample =
    CurveSegment::Impl::LocalGeodesic::LocalGeodesicSample;
using GeodesicInitialConditions =
    CurveSegment::Impl::LocalGeodesic::GeodesicInitialConditions;
using GeodesicJacobian = Vec4;
using PointVariation   = ContactGeometry::GeodesicPointVariation;
using Variation        = ContactGeometry::GeodesicVariation;
using LineSegment      = CableSpan::LineSegment;
using Status           = CurveSegment::Status;

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
//                                GEODESIC MATH TODO move to contact geometry
//==============================================================================

namespace
{

struct ImplicitGeodesicState
{
    ImplicitGeodesicState() = default;

    explicit ImplicitGeodesicState(
        Vec<10, Real>&& implicitGeodesicStateAsVector)
    {
        asVecMut() = implicitGeodesicStateAsVector;
    }

    ImplicitGeodesicState(Vec3 point, Vec3 tangent) : x(point), t(tangent){};

    const Vec<10, Real>& asVec() const
    {
        return reinterpret_cast<const Vec<10, Real>&>(x[0]);
    }

    Vec<10, Real>& asVecMut()
    {
        return reinterpret_cast<Vec<10, Real>&>(x[0]);
    }

    Vec<10, Real> calcDerivativeVector(
        const Vec3& acceleration,
        Real gaussianCurvature) const
    {
        ImplicitGeodesicState dy;
        dy.x    = t;
        dy.t    = acceleration;
        dy.a    = aDot;
        dy.aDot = -a * gaussianCurvature;
        dy.r    = rDot;
        dy.rDot = -r * gaussianCurvature;
        return {dy.asVec()};
        ;
    }

    Vec3 x    = {NaN, NaN, NaN};
    Vec3 t    = {NaN, NaN, NaN};
    Real a    = 1.;
    Real aDot = 0.;
    Real r    = 0.;
    Real rDot = 1.;
};

void calcSurfaceProjectionFast(
    const ContactGeometry& geometry,
    Vec3& x,
    Vec3& t,
    size_t maxIter,
    Real eps)
{
    size_t it = 0;
    for (; it < maxIter; ++it) {
        const Real c = geometry.calcSurfaceValue(x);

        if (std::abs(c) < eps) {
            break;
        }

        const Vec3 g = geometry.calcSurfaceGradient(x);
        x += -g * c / dot(g, g);
    }

    SimTK_ASSERT(
        it < maxIter,
        "Surface projection failed: Reached max iterations");

    UnitVec3 n(geometry.calcSurfaceGradient(x));
    t         = t - dot(n, t) * n;
    Real norm = t.norm();
    SimTK_ASSERT(!isNaN(norm), "Surface projection failed: Detected NaN");
    SimTK_ASSERT(
        norm > 1e-13,
        "Surface projection failed: Tangent guess is parallel to surface "
        "normal");
    t = t / norm;
}

Real calcPointOnLineNearOriginAsFactor(Vec3 a, Vec3 b)
{
    const Vec3 e = b - a;
    Real c       = -dot(a, e) / dot(e, e);
    return std::max(0., std::min(1., c));
};

Real calcPointOnLineNearPointAsFactor(Vec3 a, Vec3 b, Vec3 point)
{
    return calcPointOnLineNearOriginAsFactor(a - point, b - point);
};

bool calcNearestPointOnLineImplicitly(
    const ContactGeometry& geometry,
    Vec3 a,
    Vec3 b,
    Vec3& point,
    size_t maxIter,
    double eps)
{
    // Initial guess.
    double alpha = calcPointOnLineNearPointAsFactor(a, b, point);
    size_t iter  = 0;

    for (; iter < maxIter; ++iter) {
        // Touchdown point on line.
        const Vec3 d  = b - a;
        const Vec3 pl = a + (b - a) * alpha;

        // Constraint evaluation at touchdown point.
        const double c = geometry.calcSurfaceValue(pl);

        // Break on touchdown, TODO or not?
        if (std::abs(c) < eps)
            break;

        // Gradient at point on line.
        const Vec3 g  = geometry.calcSurfaceGradient(pl);
        const Mat33 H = geometry.calcSurfaceHessian(pl);

        // Add a weight to the newton step to avoid large steps.
        constexpr double w = 0.5;

        // Update alpha.
        const double step = dot(g, d) / (dot(d, H * d) + w);

        // Stop when converged.
        if (std::abs(step) < eps)
            break;

        // Clamp the stepsize.
        constexpr double maxStep = 0.25;
        alpha -= std::min(std::max(-maxStep, step), maxStep);

        // Stop when leaving bounds.
        if (alpha < 0. || alpha > 1.)
            break;
    }

    // Write the point on line nearest the surface.
    alpha = std::max(std::min(1., alpha), 0.);
    point = a + (b - a) * alpha;

    // Assumes a negative constraint evaluation means touchdown.
    const bool contact = geometry.calcSurfaceValue(point) < eps;

    // TODO handle here?
    if (iter >= maxIter) {
        std::cout << "a = " << a << "\n";
        std::cout << "b = " << b << "\n";
        std::cout << "p = " << point << "\n";
        std::cout << "c = " << alpha << "\n";
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Failed to compute point on line nearest "
                                 "surface: Reached max iterations");
    }

    // Return number of iterations required.
    return contact;
}

ImplicitGeodesicState operator-(
    const ImplicitGeodesicState& lhs,
    const ImplicitGeodesicState& rhs)
{
    ImplicitGeodesicState out;
    out.asVecMut() = lhs.asVec() - rhs.asVec();
    return out;
}

ImplicitGeodesicState operator+(
    const ImplicitGeodesicState& lhs,
    const Vec<10, Real>& rhs)
{
    ImplicitGeodesicState out;
    out.asVecMut() = lhs.asVec() + rhs;
    return out;
}

Real calcInfNorm(const ImplicitGeodesicState& q)
{
    Real infNorm = 0.;
    for (size_t r = 0; r < q.asVec().nrow(); ++r) {
        infNorm = std::max(infNorm, q.asVec()[r]);
    }
    return infNorm;
}

class RKM
{
    using Y  = ImplicitGeodesicState;
    using DY = Vec<10, Real>;

public:
    RKM() = default;

    // Integrate y0, populating y1, y2.
    Real step(Real h, std::function<DY(const Y&)>& f);

    struct Sample
    {
        Real l;
        FrenetFrame K;
    };

    Y stepTo(
        Y y0,
        Real x1,
        Real h0,
        std::function<DY(const Y&)>& f,         // Dynamics
        std::function<void(Y&)>& g,             // Surface projection
        std::function<void(Real, const Y&)>& m, // State log
        Real accuracy);

    size_t getNumberOfFailedSteps() const
    {
        return _failedCount;
    }

    size_t getInitStepSize() const
    {
        return _h0;
    }

private:
    static constexpr size_t ORDER = 5;

    std::array<DY, ORDER> _k{};
    // [yInit, ySave, yStep]
    std::array<Y, 3> _y{};

    Real _hMin = 1e-10;
    Real _hMax = 1e-1;
    Real _h    = NaN;
    Real _h0   = NaN;
    Real _e    = NaN;

    size_t _failedCount = 0;
};

Real RKM::step(Real h, std::function<DY(const Y&)>& f)
{
    Y& yk = _y.at(1);

    for (size_t i = 0; i < 5; ++i) {
        yk = _y.at(0) + (h / 3.) * _k.at(0);
    }

    // k1
    _k.at(0) = f(_y.at(0));

    // k2
    {
        yk       = _y.at(0) + (h / 3.) * _k.at(0);
        _k.at(1) = f(yk);
    }

    // k3
    {
        yk       = _y.at(0) + (h / 6.) * _k.at(0) + (h / 6.) * _k.at(1);
        _k.at(2) = f(yk);
    }

    // k4
    {
        yk = _y.at(0) + (1. / 8. * h) * _k.at(0) + (3. / 8. * h) * _k.at(2);
        _k.at(3) = f(yk);
    }

    // k5
    {
        yk = _y.at(0) + (1. / 2. * h) * _k.at(0) + (-3. / 2. * h) * _k.at(2) +
             (2. * h) * _k.at(3);
        _k.at(4) = f(yk);
    }

    // y1: Auxiliary --> Already updated in k5 computation.

    // y2: Final state.
    _y.at(2) = _y.at(0) + (1. / 6. * h) * _k.at(0) + (2. / 3. * h) * _k.at(3) +
               (1. / 6. * h) * _k.at(4);

    return calcInfNorm(_y.at(1) - _y.at(2)) * 0.2;
}

RKM::Y RKM::stepTo(
    Y y0,
    Real x1,
    Real h0,
    std::function<DY(const Y&)>& f,
    std::function<void(Y&)>& g,
    std::function<void(Real, const Y&)>& m,
    Real accuracy)
{
    g(y0);
    m(0., y0);

    _y.at(0)     = std::move(y0);
    _h           = h0;
    _e           = 0.;
    _failedCount = 0;

    Real x = 0.;
    while (x < x1 - 1e-13) {
        const bool init = x == 0.;

        _h = x + _h > x1 ? x1 - x : _h;

        // Attempt step.
        Real err = step(_h, f);

        // Reject if accuracy was not met.
        if (err > accuracy) { // Rejected
            // Decrease stepsize.
            _h /= 2.;
            _y.at(1) = _y.at(0);
            _y.at(2) = _y.at(0);
            ++_failedCount;
        } else {         // Accepted
            g(_y.at(2)); // Enforce constraints.
            _y.at(0) = _y.at(2);
            _y.at(1) = _y.at(2);
            x += _h;
            m(x, _y.at(0));
        }

        // Potentially increase stepsize.
        if (err < accuracy / 64.) {
            _h *= 2.;
        }

        _e = std::max(_e, err);

        SimTK_ASSERT(
            _h > _hMin,
            "Geodesic Integrator failed: Reached very small stepsize");
        SimTK_ASSERT(
            _h < _hMax,
            "Geodesic Integrator failed: Reached very large stepsize");

        if (init) {
            _h0 = _h;
        }
    }
    return _y.at(0);
}

Real calcNormalCurvature(
    const ContactGeometry& geometry,
    Vec3 point,
    Vec3 tangent)
{
    const Vec3& p  = point;
    const Vec3& v  = tangent;
    const Vec3 g   = geometry.calcSurfaceGradient(p);
    const Vec3 h_v = geometry.calcSurfaceHessian(p) * v;
    // Sign flipped compared to thesis: kn = negative, see eq 3.63
    return -dot(v, h_v) / g.norm();
}

Real calcGeodesicTorsion(
    const ContactGeometry& geometry,
    Vec3 point,
    Vec3 tangent)
{
    // TODO verify this!
    const Vec3& p  = point;
    const Vec3& v  = tangent;
    const Vec3 g   = geometry.calcSurfaceGradient(p);
    const Vec3 h_v = geometry.calcSurfaceHessian(p) * v;
    const Vec3 gxv = cross(g, v);
    return -dot(h_v, gxv) / dot(g, g);
}

UnitVec3 calcSurfaceNormal(const ContactGeometry& geometry, Vec3 point)
{
    const Vec3& p       = point;
    const Vec3 gradient = geometry.calcSurfaceGradient(p);
    return UnitVec3{gradient};
}

Vec3 calcAcceleration(const ContactGeometry& geometry, Vec3 point, Vec3 tangent)
{
    // TODO Writing it out saves a root, but optimizers are smart.
    // Sign flipped compared to thesis: kn = negative, see eq 3.63
    return calcNormalCurvature(geometry, point, std::move(tangent)) *
           calcSurfaceNormal(geometry, point);
}

Mat33 calcAdjoint(const Mat33& mat)
{
    Real fxx = mat(0, 0);
    Real fyy = mat(1, 1);
    Real fzz = mat(2, 2);

    Real fxy = mat(0, 1);
    Real fxz = mat(0, 2);
    Real fyz = mat(1, 2);

    std::array<Real, 9> elements = {
        fyy * fzz - fyz * fyz,
        fyz * fxz - fxy * fzz,
        fxy * fyz - fyy * fxz,
        fxz * fyz - fxy * fzz,
        fxx * fzz - fxz * fxz,
        fxy * fxz - fxx * fyz,
        fxy * fyz - fxz * fyy,
        fxy * fxz - fxx * fyz,
        fxx * fyy - fxy * fxy};
    Mat33 adj;
    size_t i = 0;
    for (size_t r = 0; r < 3; ++r) {
        for (size_t c = 0; c < 3; ++c) {
            adj(r, c) = elements[i];
            ++i;
        }
    }
    return adj;
}

Real calcGaussianCurvature(const ContactGeometry& geometry, Vec3 point)
{
    const Vec3& p = point;
    Vec3 g        = geometry.calcSurfaceGradient(p);
    Real gDotg    = dot(g, g);
    Mat33 adj     = calcAdjoint(geometry.calcSurfaceHessian(p));

    if (gDotg * gDotg < 1e-13) {
        throw std::runtime_error(
            "Gaussian curvature inaccurate: are we normal to surface?");
    }

    return (dot(g, adj * g)) / (gDotg * gDotg);
}

void calcFrenetFrame(
    const ContactGeometry& geometry,
    const ImplicitGeodesicState& q,
    FrenetFrame& K)
{
    K.setP(q.x);
    K.updR().setRotationFromTwoAxes(
        calcSurfaceNormal(geometry, q.x),
        NormalAxis,
        q.t,
        TangentAxis);
}

void calcGeodesicBoundaryState(
    const ContactGeometry& geometry,
    const ImplicitGeodesicState& q,
    bool isEnd,
    FrenetFrame& K,
    Mat34& v,
    Mat34& w)
{
    calcFrenetFrame(geometry, q, K);

    const UnitVec3& t = K.R().getAxisUnitVec(TangentAxis);
    const UnitVec3& n = K.R().getAxisUnitVec(NormalAxis);
    const UnitVec3& b = K.R().getAxisUnitVec(BinormalAxis);

    // TODO remove fourth element?
    v.col(0) = t;
    v.col(1) = b * q.a;
    v.col(2) = b * q.r;
    v.col(3) = isEnd ? v.col(0) : Vec3{0.};

    const Real tau_g   = geometry.calcGeodesicTorsion(q.x, t);
    const Real kappa_n = geometry.calcNormalCurvature(q.x, t);
    const Real kappa_a = geometry.calcNormalCurvature(q.x, b);

    w.col(0) = tau_g * t + kappa_n * b;
    w.col(1) = -q.a * kappa_a * t - q.aDot * n - q.a * tau_g * b;
    w.col(2) = -q.r * kappa_a * t - q.rDot * n - q.r * tau_g * b;
    w.col(3) = isEnd ? w.col(0) : Vec3{0.};
}

void calcGeodesicAndVariationImplicitly(
    const ContactGeometry& geometry,
    Vec3 x,
    Vec3 t,
    Real l,
    Real& ds,
    FrenetFrame& K_P,
    Variation& dK_P,
    FrenetFrame& K_Q,
    Variation& dK_Q,
    Real accuracy,
    size_t prjMaxIter,
    Real prjAccuracy,
    std::vector<LocalGeodesicSample>& log)
{
    using Y  = ImplicitGeodesicState;
    using DY = Vec<10, Real>;

    std::function<DY(const Y&)> f = [&](const Y& q) -> DY {
        return q.calcDerivativeVector(
            calcAcceleration(geometry, q.x, q.t),
            calcGaussianCurvature(geometry, q.x));
    };
    std::function<void(Y&)> g = [&](Y& q) {
        calcSurfaceProjectionFast(geometry, q.x, q.t, prjMaxIter, prjAccuracy);
    };
    std::function<void(Real, const Y&)> m = [&](Real l, const Y& q) {
        FrenetFrame frame;
        calcFrenetFrame(geometry, q, frame);
        log.emplace_back(l, frame);
    };

    Y y0(x, t);

    RKM rkm{};

    Y y1 = rkm.stepTo(y0, l, ds, f, g, m, accuracy);

    SimTK_ASSERT(log.size() > 0, "Failed to integrate geodesic: Log is empty");
    calcGeodesicBoundaryState(geometry, y0, false, K_P, dK_P[1], dK_P[0]);
    calcGeodesicBoundaryState(geometry, y1, true, K_Q, dK_Q[1], dK_Q[0]);

    ds = rkm.getInitStepSize();
}

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
        sizeof(Correction) == sizeof(Real) * GeodesicDOF,
        "Invalid size of corrections vector");
    SimTK_ASSERT(
        data.pathCorrection.size() * sizeof(Real) == n * sizeof(Correction),
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

    // Take the length correction, and add to the current length.
    Real dl = c[3]; // Length increment is the last correction element.
    g0.l    = std::max(l + dl, 0.); // Clamp length to be nonnegative.

    return g0;
}

GeodesicInitialConditions GeodesicInitialConditions::
    CreateFromGroundInSurfaceFrame(
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
    shootNewGeodesic(g0, cache);
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

void LocalGeodesic::applyGeodesicCorrection(const State& s, const Correction& c)
    const
{
    realizePosition(s);

    // Get the previous geodesic.
    const CacheEntry& g = getCacheEntry(s);

    // Get corrected initial conditions.
    const GeodesicInitialConditions g0 =
        GeodesicInitialConditions::CreateCorrected(g.K_P, g.dK_P, g.length, c);

    // Shoot the new geodesic.
    shootNewGeodesic(g0, updCacheEntry(s));
}

void LocalGeodesic::calcPathPoints(const State& s, std::vector<Vec3>& points)
    const
{
    realizePosition(s);
    const CacheEntry& cache = getCacheEntry(s);

    if (analyticFormAvailable()) {
        throw std::runtime_error("NOTYETIMPLEMENTED");
        m_Geometry.resampleGeodesicPointsAnalytically(
            cache.K_P,
            cache.K_Q,
            cache.length,
            m_NumberOfAnalyticPoints,
            points);
    } else {
        for (const LocalGeodesicSample& sample : cache.samples) {
            points.push_back(sample.frame.p());
        }
    }
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
    if (dot(prev_QS - g.K_P.p(), g.K_P.R().getAxisUnitVec(NormalAxis)) <= 0. &&
        dot(next_PS - g.K_P.p(), g.K_P.R().getAxisUnitVec(NormalAxis)) <= 0.) {
        // No liftoff.
        return;
    }

    // Liftoff detected: update status.
    g.status = Status::Liftoff;
    // Initialize the tracking point from the last geodesic start point.
    cache.trackingPointOnLine = g.K_P.p();
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
    bool touchdownDetected;
    if (m_Geometry.analyticFormAvailable()) {
        throw std::runtime_error("NOTYETIMPLEMENTED");
        /* touchdownDetected = m_Geometry.calcNearestPointOnLineAnalytically( */
        /*     prev_QS, */
        /*     next_PS, */
        /*     cache.trackingPointOnLine); */
    } else {
        touchdownDetected = calcNearestPointOnLineImplicitly(
            m_Geometry,
            prev_QS,
            next_PS,
            cache.trackingPointOnLine,
            m_TouchdownIter,
            m_TouchdownAccuracy);
    }
    if (!touchdownDetected) {
        return;
    }

    // Touchdown detected: Remove the liftoff status flag.
    g.status = Status::Ok;
    // Shoot a zero length geodesic at the touchdown point.
    GeodesicInitialConditions g0 = GeodesicInitialConditions::CreateAtTouchdown(
        prev_QS,
        next_PS,
        cache.trackingPointOnLine);
    shootNewGeodesic(g0, cache);
}

void LocalGeodesic::assertSurfaceBounds(
    const Vec3& prev_QS,
    const Vec3& next_PS) const
{
    // Make sure that the previous point does not lie inside the surface.
    SimTK_ASSERT(
        m_Geometry.calcSurfaceValue(prev_QS) < 0.,
        "Unable to wrap over surface: Preceding point lies inside the surface");
    SimTK_ASSERT(
        m_Geometry.calcSurfaceValue(next_PS) < 0.,
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

void LocalGeodesic::shootNewGeodesic(
    const GeodesicInitialConditions& g0,
    CacheEntry& cache) const
{
    calcGeodesicAndVariationImplicitly(
        m_Geometry,
        g0.x,
        g0.t,
        g0.l,
        cache.sHint,
        cache.K_P,
        cache.dK_P,
        cache.K_Q,
        cache.dK_Q,
        m_IntegratorAccuracy,
        m_ProjectionMaxIter,
        m_ProjectionAccuracy,
        cache.samples);

    /* using Y  = ImplicitGeodesicState; */
    /* if (m_Geometry.analyticFormAvailable()) { */

    /*     m_Geometry.calcGeodesicWithVariationAnalytically( */
    /*         g0.x, */
    /*         g0.t, */
    /*         g0.l, */
    /*         cache.KP, */
    /*         cache.dKP, */
    /*         cache.KQ, */
    /*         cache.dKQ); */
    /* } else { */
    /*     // Compute geodesic start boundary frame and variation. */
    /*     m_Geometry.calcNearestFrenetFrameImplicitlyFast( */
    /*         g0.x, */
    /*         g0.t, */
    /*         cache.KP, */
    /*         m_ProjectionMaxIter, */
    /*         m_ProjectionRequiredAccuracy); */
    /*     m_Geometry.calcGeodesicStartFrameVariationImplicitly( */
    /*         cache.KP, */
    /*         cache.dKP); */
    /*     // Compute geodesic end boundary frame amd variation (shoot new */
    /*     // geodesic). */
    /*     m_Geometry.calcGeodesicEndFrameVariationImplicitly( */
    /*         cache.KP, */
    /*         g0.l, */
    /*         cache.KQ, */
    /*         cache.dKQ, */
    /*         cache.sHint, */
    /*         m_IntegratorAccuracy, */
    /*         cache.frames); */
    /* } */
}

//==============================================================================
//                               CURVE SEGMENT
//==============================================================================

CurveSegment::CurveSegment(
    CableSpan cable,
    const MobilizedBody& mobod,
    Transform X_BS,
    const ContactGeometry& geometry,
    Vec3 xHint) :
    m_Impl(std::shared_ptr<CurveSegment::Impl>(
        new CurveSegment::Impl(cable, mobod, X_BS, geometry, xHint)))
{
    // TODO bit awkward to set the index later.
    updImpl().setIndex(cable.adoptSegment(*this));
    /* CurveSegmentIndex ix = cable.adoptSegment(*this); */
    /* m_Impl = std::shared_ptr<CurveSegment::Impl>(new
     * CurveSegment::Impl(cable, ix, mobod, X_BS, geometry, xHint)); */
}

const CableSpan& CurveSegment::getCable() const
{
    return getImpl().getCable();
}

Real CurveSegment::getSegmentLength(const State& s)
{
    getImpl().realizeCablePosition(s);
    return getImpl().getPosInfo(s).length;
}

CurveSegment::Status CurveSegment::getStatus(const State& s) const
{
    getImpl().realizeCablePosition(s);
    return getImpl().getStatus(s);
}

//==============================================================================
//                          CURVE SEGMENT IMPL
//==============================================================================

CurveSegment::Impl::Impl(
    CableSpan path,
    const MobilizedBody& mobod,
    const Transform& X_BS,
    ContactGeometry geometry,
    Vec3 initPointGuess) :
    m_Subsystem(path.getImpl().getSubsystem()),
    m_Path(path),
    m_Index(-1), // TODO what to do with this index, and when
    m_Mobod(mobod),
    m_Offset(X_BS),
    m_Geodesic(m_Subsystem, geometry, initPointGuess)
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
        m_Geodesic.getInitialPointGuess(); // TODO move into function call?

    GeodesicInitialConditions g0 =
        GeodesicInitialConditions::CreateZeroLengthGuess(prev_QS, xGuess_S);
    m_Geodesic.calcInitialGeodesic(s, g0);

    m_Subsystem.markCacheValueNotRealized(s, m_PosInfoIx);
}

void CurveSegment::Impl::applyGeodesicCorrection(
    const State& s,
    const CurveSegment::Impl::Correction& c) const
{
    // Apply correction to curve.
    m_Geodesic.applyGeodesicCorrection(s, c);

    // Invalidate position level cache.
    m_Subsystem.markCacheValueNotRealized(s, m_PosInfoIx);
}

void CurveSegment::Impl::calcPathPoints(
    const State& s,
    std::vector<Vec3>& points) const
{
    const Transform& X_GS = getPosInfo(s).X_GS;
    size_t initSize       = points.size();
    m_Geodesic.calcPathPoints(s, points);
    for (size_t i = points.size() - initSize; i < points.size(); ++i) {
        points.at(i) = X_GS.shiftFrameStationToBase(points.at(i));
    }
}

namespace
{
void xformSurfaceGeodesicToGround(
    const LocalGeodesicInfo& geodesic_S,
    const Transform& X_GS,
    CurveSegment::Impl::PosInfo& geodesic_G)
{
    geodesic_G.X_GS = X_GS;

    // Store the the local geodesic in ground frame.
    geodesic_G.KP = X_GS.compose(geodesic_S.K_P);
    geodesic_G.KQ = X_GS.compose(geodesic_S.K_Q);

    geodesic_G.dKP[0] = X_GS.R() * geodesic_S.dK_P[0];
    geodesic_G.dKP[1] = X_GS.R() * geodesic_S.dK_P[1];

    geodesic_G.dKQ[0] = X_GS.R() * geodesic_S.dK_Q[0];
    geodesic_G.dKQ[1] = X_GS.R() * geodesic_S.dK_Q[1];

    // TODO use SpatialVec for variation.
    /* for (size_t i = 0; i < GeodesicDOF; ++i) { */
    /*     geodesic_G.dKP[i][0] = X_GS.xformFrameVecToBase(geodesic_S.dKP[i][0])
     */
    /*     geodesic_G.dKP[i][1] = X_GS.xformFrameVecToBase(geodesic_S.dKP[i][1])
     */

    /*     geodesic_G.dKQ[i][0] = X_GS.xformFrameVecToBase(geodesic_S.dKQ[i][0])
     */
    /*     geodesic_G.dKQ[i][1] = X_GS.xformFrameVecToBase(geodesic_S.dKQ[i][1])
     */
    /* } */

    geodesic_G.length = geodesic_S.length;

    throw std::runtime_error("NOTYETIMPLEMENTED: Check transformation order");
}
} // namespace

void CurveSegment::Impl::calcPosInfo(const State& s, PosInfo& posInfo) const
{
    if (m_Geodesic.getStatus(s) == Status::Disabled) {
        return;
    }

    // Compute tramsform from local surface frame to ground.
    const Transform& X_GS = calcSurfaceFrameInGround(s);

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
        m_Geodesic.calcLocalGeodesicInfo(s, prev_S, next_S);

    // Store the the local geodesic in ground frame.
    xformSurfaceGeodesicToGround(geodesic_S, X_GS, updPosInfo(s));
}

SpatialVec CurveSegment::Impl::calcAppliedWrenchInGround(
    const State& s,
    Real tension) const
{
    const PosInfo& posInfo = getPosInfo(s);

    const UnitVec3& t_P = posInfo.KP.R().getAxisUnitVec(TangentAxis);
    const UnitVec3& t_Q = posInfo.KQ.R().getAxisUnitVec(TangentAxis);
    const Vec3 F_G      = tension * (t_Q - t_P);

    const Vec3 x_GS = posInfo.X_GS.p();
    const Vec3& r_P = posInfo.KP.p() - x_GS;
    const Vec3& r_Q = posInfo.KQ.p() - x_GS;
    const Vec3 M_G  = tension * (r_Q % t_Q - r_P % t_P);
    return {M_G, F_G};
}

void CurveSegment::Impl::applyBodyForce(
    const State& s,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    // TODO why?
    if (tension <= 0) {
        return;
    }

    if (getStatus(s) != Status::Ok) {
        return;
    }

    m_Mobod.applyBodyForce(
        s,
        calcAppliedWrenchInGround(s, tension),
        bodyForcesInG);
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

Real calcPathError(const LineSegment& e, const Rotation& R, CoordinateAxis axis)
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
//                                CABLE SPAN
//==============================================================================

CableSpan::CableSpan(
    CableSubsystem& subsystem,
    const MobilizedBody& originBody,
    const Vec3& defaultOriginPoint,
    const MobilizedBody& terminationBody,
    const Vec3& defaultTerminationPoint) :
    m_Impl(std::shared_ptr<Impl>(new Impl(
        subsystem,
        originBody,
        defaultOriginPoint,
        terminationBody,
        defaultTerminationPoint)))
{}

CurveSegmentIndex CableSpan::adoptSegment(const CurveSegment& segment)
{
    return updImpl().adoptSegment(segment);
}

int CableSpan::getNumCurveSegments() const
{
    return getImpl().getNumCurveSegments();
}

const CurveSegment& CableSpan::getCurveSegment(CurveSegmentIndex ix) const
{
    return getImpl().getCurveSegment(ix);
}

Real CableSpan::getLength(const State& s) const
{
    return getImpl().getPosInfo(s).l;
}

Real CableSpan::getLengthDot(const State& s) const
{
    return getImpl().getVelInfo(s).lengthDot;
}
void CableSpan::applyBodyForces(
    const State& s,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    return getImpl().applyBodyForces(s, tension, bodyForcesInG);
}

//==============================================================================
//                              CABLE SPAN IMPL
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

Real CableSpan::Impl::calcPathLength(
    const State& s,
    const std::vector<LineSegment>& lines) const
{
    Real lTot = 0.;
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
        try {
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

        } catch (const std::exception& e) {
            dataCache.unlock();
            throw e;
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
    const State& s,
    Real tension,
    Vector_<SpatialVec>& bodyForcesInG) const
{
    // TODO why?
    if (tension <= 0) {
        return;
    }

    realizePosition(s);

    for (const CurveSegment& segment : m_CurveSegments) {
        segment.getImpl().applyBodyForce(s, tension, bodyForcesInG);
    }
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

const CableSpan& CableSubsystem::getPath(CableSpanIndex ix) const
{
    return getImpl().getCablePath(ix);
}

CableSpan& CableSubsystem::updPath(CableSpanIndex ix)
{
    return updImpl().updCablePath(ix);
}
