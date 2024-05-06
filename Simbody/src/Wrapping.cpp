#include "SimTKcommon/internal/CoordinateAxis.h"
#include "SimTKcommon/internal/ExceptionMacros.h"
#include "SimTKmath.h"
#include "WrappingImpl.h"
#include "simmath/internal/ContactGeometry.h"
#include <cstddef>
#include <memory>
#include <stdexcept>

using namespace SimTK;

using Correction          = ContactGeometry::GeodesicCorrection;
using FrameVariation      = ContactGeometry::GeodesicFrameVariation;
using FrenetFrame         = ContactGeometry::FrenetFrame;
using GeodesicInfo        = CurveSegment::Impl::PosInfo;
using GeodesicJacobian    = Vec4;
using LineSegment         = CableSpan::LineSegment;
using LocalGeodesicInfo   = CurveSegment::Impl::InstanceEntry;
using LocalGeodesicSample = CurveSegment::Impl::LocalGeodesicSample;
using PointVariation      = ContactGeometry::GeodesicPointVariation;
using SolverData          = CableSubsystem::Impl::SolverData;
using Status              = CurveSegment::Status;
using Variation           = ContactGeometry::GeodesicVariation;

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
    updImpl().setIndex(cable.updImpl().adoptSegment(*this));
}

const CableSpan& CurveSegment::getCable() const
{
    return getImpl().getCable();
}

Real CurveSegment::getSegmentLength(const State& s) const
{
    getImpl().realizeCablePosition(s);
    return getImpl().getInstanceEntry(s).length;
}

const ContactGeometry& CurveSegment::getContactGeometry() const
{
    return getImpl().getContactGeometry();
}

const Mobod& CurveSegment::getMobilizedBody() const
{
    return getImpl().getMobilizedBody();
}

void CurveSegment::setXformSurfaceToBody(Transform X_BS)
{
    updImpl().setXformSurfaceToBody(std::move(X_BS));
}

const Transform& CurveSegment::getXformSurfaceToBody() const
{
    return getImpl().getXformSurfaceToBody();
}

CurveSegment::Status CurveSegment::getStatus(const State& s) const
{
    getImpl().realizeCablePosition(s);
    return getImpl().getInstanceEntry(s).status;
}

const Transform& CurveSegment::getFrenetFrameStart(const State& state) const
{
    getImpl().realizeCablePosition(state);
    return getImpl().getPosInfo(state).KP;
}

const Transform& CurveSegment::getFrenetFrameEnd(const State& state) const
{
    getImpl().realizeCablePosition(state);
    return getImpl().getPosInfo(state).KQ;
}

int CurveSegment::getNumberOfIntegratorStepsTaken(const State& state)
{
    getImpl().realizeCablePosition(state);
    return getImpl().getInstanceEntry(state).samples.size();
}

Real CurveSegment::getInitialIntegratorStepSize(const State& state)
{
    getImpl().realizeCablePosition(state);
    return getImpl().getInstanceEntry(state).sHint;
}

void CurveSegment::calcUnitForce(const State& state, SpatialVec& unitForce_G)
    const
{
    getImpl().realizeCablePosition(state);
    getImpl().calcUnitForce(state, unitForce_G);
}

int CurveSegment::calcPoints(
    const State& state,
    std::vector<Vec3>& points_G,
    int nPoints) const
{
    getImpl().realizeCablePosition(state);
    return getImpl().calcPathPoints(state, points_G, nPoints);
}

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
{
    subsystem.updImpl().adoptCablePath(*this);
}

void CableSpan::adoptWrappingObstacle(
    const MobilizedBody& mobod,
    Transform X_BS,
    const ContactGeometry& geometry,
    Vec3 contactPointHint)
{
    CurveSegment(*this, mobod, X_BS, geometry, contactPointHint);
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

int CableSpan::calcPoints(
    const State& state,
    std::vector<Vec3>& points_G,
    int nPointsPerCurveSegment) const
{
    return getImpl().calcPathPoints(state, points_G, nPointsPerCurveSegment);
}

void CableSpan::calcUnitForceAtOrigin(
    const State& state,
    SpatialVec& unitForce_G) const
{
    getImpl().calcUnitForceAtOrigin(state, unitForce_G);
}

void CableSpan::calcUnitForceAtTermination(
    const State& state,
    SpatialVec& unitForce_G) const
{
    getImpl().calcUnitForceAtTermination(state, unitForce_G);
}

Real CableSpan::calcCablePower(const State& state, Real tension) const
{
    return getImpl().calcCablePower(state, tension);
}

//==============================================================================
//                                CONSTANTS
//==============================================================================
namespace
{
static const int GeodesicDOF = 4;
} // namespace

//==============================================================================
//                                GEODESIC MATH TODO move to contact geometry
//==============================================================================

namespace
{

// TODO Sign of implicit function was flipped.
Real calcSurfaceConstraintValue(const ContactGeometry& geometry, Vec3 point)
{
    return -geometry.calcSurfaceValue(point);
}

// TODO Sign of implicit function was flipped.
Vec3 calcSurfaceConstraintGradient(const ContactGeometry& geometry, Vec3 point)
{
    return -geometry.calcSurfaceGradient(point);
}

// TODO Sign of implicit function was flipped.
Mat33 calcSurfaceConstraintHessian(const ContactGeometry& geometry, Vec3 point)
{
    return -geometry.calcSurfaceHessian(point);
}

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
        const Real c = calcSurfaceConstraintValue(geometry, x);

        if (std::abs(c) < eps) {
            break;
        }

        const Vec3 g = calcSurfaceConstraintGradient(geometry, x);
        x += -g * c / dot(g, g);
    }

    if (it >= maxIter) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "Surface projection failed: Reached max iterations");
    }

    if (std::abs(geometry.calcSurfaceValue(x)) > eps) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "Surface projection failed: no longer on surface");
    }

    UnitVec3 n(calcSurfaceConstraintGradient(geometry, x));
    t         = t - dot(n, t) * n;
    Real norm = t.norm();
    if (isNaN(norm))
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Surface projection failed: Detected NaN");
    if (norm < 1e-13)
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Surface projection failed: Tangent guess is "
                                 "parallel to surface normal");

    t = t / norm;
}

bool calcNearestPointOnLineImplicitly(
    const ContactGeometry& geometry,
    Vec3 a,
    Vec3 b,
    Vec3& point,
    size_t maxIter,
    Real eps)
{
    // Initial guess.
    const Vec3 d = b - a;
    Real alpha   = -dot(d, a - point) / dot(d, d);
    alpha        = std::max(0., std::min(alpha, 1.));

    size_t iter = 0;

    for (; iter < maxIter; ++iter) {
        // Touchdown point on line.
        const Vec3 pl = a + d * alpha;

        // Constraint evaluation at touchdown point.
        const Real c = calcSurfaceConstraintValue(geometry, pl);

        // Break on touchdown, TODO or not?
        if (std::abs(c) < eps)
            break;

        // Gradient at point on line.
        const Vec3 g  = calcSurfaceConstraintGradient(geometry, pl);
        const Mat33 H = calcSurfaceConstraintHessian(geometry, pl);

        // Add a weight to the newton step to avoid large steps.
        constexpr Real w = 0.5;

        // Update alpha.
        const Real step = dot(g, d) / (dot(d, H * d) + w);

        // Stop when converged.
        if (std::abs(step) < eps)
            break;

        // Clamp the stepsize.
        constexpr Real maxStep = 0.25;
        alpha -= std::min(std::max(-maxStep, step), maxStep);

        // Stop when leaving bounds.
        if (alpha < 0. || alpha > 1.)
            break;
    }

    // Write the point on line nearest the surface.
    alpha = std::max(0., std::min(alpha, 1.));
    point = a + d * alpha;

    // Assumes a negative constraint evaluation means touchdown.
    const bool contact = calcSurfaceConstraintValue(geometry, point) < eps;

    // TODO handle here?
    if (iter >= maxIter) {
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
    Real infNorm           = 0.;
    const Vec<10, Real>& v = q.asVec();
    for (size_t r = 0; r < v.nrow(); ++r) {
        infNorm = std::max(infNorm, std::abs(v[r]));
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
        Real& h0,
        std::function<DY(const Y&)>& f,         // Dynamics
        std::function<void(Y&)>& g,             // Surface projection
        std::function<void(Real, const Y&)>& m, // State log
        Real accuracy);

    size_t getNumberOfFailedSteps() const
    {
        return _failedCount;
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
    Real& h0,
    std::function<DY(const Y&)>& f,
    std::function<void(Y&)>& g,
    std::function<void(Real, const Y&)>& m,
    Real accuracy)
{
    g(y0);
    m(0., y0);

    _y.at(0)     = std::move(y0);
    _h           = std::min(h0 > _hMin ? h0 : _hMin, _hMax);
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

        if (_h < _hMin) {
            // TODO use SimTK_ASSERT
            throw std::runtime_error(
                "Geodesic Integrator failed: Reached very small stepsize");
        }

        if (init) {
            h0 = _h;
        }
    }
    if (std::abs(x - x1) > 1e-13) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("failed to integrate");
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
    const Vec3 g   = calcSurfaceConstraintGradient(geometry, p);
    const Vec3 h_v = calcSurfaceConstraintHessian(geometry, p) * v;
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
    const Vec3 g   = calcSurfaceConstraintGradient(geometry, p);
    const Vec3 h_v = calcSurfaceConstraintHessian(geometry, p) * v;
    const Vec3 gxv = cross(g, v);
    return -dot(h_v, gxv) / dot(g, g);
}

UnitVec3 calcSurfaceNormal(const ContactGeometry& geometry, Vec3 point)
{
    const Vec3& p       = point;
    const Vec3 gradient = calcSurfaceConstraintGradient(geometry, p);
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
    Vec3 g        = calcSurfaceConstraintGradient(geometry, p);
    Real gDotg    = dot(g, g);
    Mat33 adj     = calcAdjoint(calcSurfaceConstraintHessian(geometry, p));

    if (gDotg * gDotg < 1e-13) {
        // TODO use SimTK_ASSERT
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
    v.col(0) = Vec3{t};
    v.col(1) = b * q.a;
    v.col(2) = b * q.r;
    v.col(3) = isEnd ? v.col(0) : Vec3{0.};

    const Real tau_g   = calcGeodesicTorsion(geometry, q.x, t);
    const Real kappa_n = calcNormalCurvature(geometry, q.x, t);
    const Real kappa_a = calcNormalCurvature(geometry, q.x, b);

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
        log.push_back(LocalGeodesicSample(l, frame));
    };

    Y y0(x, t);
    g(y0);

    RKM rkm{};

    Y y1 = rkm.stepTo(y0, l, ds, f, g, m, accuracy);

    SimTK_ASSERT(log.size() > 0, "Failed to integrate geodesic: Log is empty");
    calcGeodesicBoundaryState(geometry, y0, false, K_P, dK_P[1], dK_P[0]);
    calcGeodesicBoundaryState(geometry, y1, true, K_Q, dK_Q[1], dK_Q[0]);
}

} // namespace

//==============================================================================
//                                Resampling geodesic
//==============================================================================

namespace
{
Vec3 calcHermiteInterpolation(
    Real x0,
    const Vec3& y0,
    const Vec3& y0Dot,
    Real x1,
    const Vec3& y1,
    const Vec3& y1Dot,
    Real x)
{
    const Real dx    = x1 - x0;
    const Vec3 dy    = y1 - y0;
    const Vec3 dyDot = y1Dot - y0Dot;

    const Vec3& c0 = y0;
    const Vec3& c1 = y0Dot;
    const Vec3 c3 =
        -2. * (dy - y0Dot * dx - 0.5 * dx * dyDot) / std::pow(dx, 3);
    const Vec3 c2 = (dyDot / dx - 3. * c3 * dx) / 2.;

    const Real h = x - x0;
    return c0 + h * (c1 + h * (c2 + h * c3));
}

Vec3 calcHermiteInterpolation(
    const LocalGeodesicSample& a,
    const LocalGeodesicSample& b,
    Real l)
{
    return calcHermiteInterpolation(
        a.length,
        a.frame.p(),
        a.frame.R().getAxisUnitVec(TangentAxis),
        b.length,
        b.frame.p(),
        b.frame.R().getAxisUnitVec(TangentAxis),
        l);
}

size_t calcResampledGeodesicPoints(
    const std::vector<LocalGeodesicSample>& geodesic,
    const Transform& X_GS,
    int nSamples,
    std::vector<Vec3>& interpolatedSamples)
{
    // Some sanity checks.
    if (geodesic.empty()) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "Resampling of geodesic failed: Provided geodesic is empty.");
    }
    if (geodesic.front().length != 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: First frame "
                                 "must be at length = zero");
    }
    if (geodesic.front().length < 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: Last frame "
                                 "must be at length > zero");
    }
    if (nSamples == 1 && geodesic.size() != 1) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Resampling of geodesic failed: Requested "
                                 "number of samples must be unequal to 1");
    }

    // Capture the start of the geodesic.
    interpolatedSamples.push_back(
        X_GS.shiftFrameStationToBase(geodesic.front().frame.p()));

    // If there is but one sample in the geodesic, write that sample and exit.
    if (geodesic.size() == 1) {
        return 1;
    }

    // Seperate the interpolation points by equal length increments.
    const Real dl = geodesic.back().length / static_cast<Real>(nSamples - 1);

    // Compute the interpolated points from the geodesic.
    auto itGeodesic = geodesic.begin();
    // We can skip the first and last samples, because these are pushed
    // manually before and after this loop respectively (we start at i=1
    // and stop at i < nSamples-1).
    for (size_t i = 1; i < nSamples - 1; ++i) {

        // Length at the current interpolation point.
        const Real length = dl * static_cast<Real>(i);

        // Find the two samples (lhs, rhs) of the geodesic such that the
        // length of the interpolation point lies between them.
        // i.e. find: lhs.length <= length < rhs.length
        while (true) {
            // Sanity check: We should stay within range.
            if ((itGeodesic + 1) == geodesic.end()) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Resampling of geodesic failed: Attempted to read out of "
                    "array range");
            }

            // The candidate samples to use for interpolation.
            const LocalGeodesicSample& lhs = *itGeodesic;
            const LocalGeodesicSample& rhs = *(itGeodesic + 1);

            // Sanity check: Samples are assumed to be monotonically increasing
            // in length.
            if (lhs.length > rhs.length) {
                // TODO use SimTK_ASSERT
                throw std::runtime_error(
                    "Resampling of geodesic failed: Samples are not "
                    "monotonically increasing in length.");
            }

            // Check that the interpolation point lies between these samples:
            // lhs.length <= length < rhs.length
            if (length >= rhs.length) {
                // Try the next two samples.
                ++itGeodesic;
                continue;
            }

            // Do the interpolation, and write to the output buffer.
            const Vec3 point_S = calcHermiteInterpolation(lhs, rhs, length);
            // Transform to ground frame.
            const Vec3 point_G = X_GS.shiftFrameStationToBase(point_S);
            // Write interpolated point to the output buffer.
            interpolatedSamples.push_back(point_G);

            break;
        }
    }

    // Capture the last point of the geodesic.
    interpolatedSamples.push_back(
        X_GS.shiftFrameStationToBase(geodesic.back().frame.p()));

    return nSamples;
}
} // namespace

int CurveSegment::Impl::calcPathPoints(
    const State& s,
    std::vector<Vec3>& points,
    int nSamples) const
{
    const Transform& X_GS           = getPosInfo(s).X_GS;
    const InstanceEntry& geodesic_S = getInstanceEntry(s);
    if (!geodesic_S.isActive()) {
        return 0;
    }

    // Do not do any resampling if nSamples==0, simply write the points from the
    // integrator to the output buffer.
    if (nSamples == 0) {
        for (const LocalGeodesicSample& sample : geodesic_S.samples) {
            points.push_back(X_GS.shiftFrameStationToBase(sample.frame.p()));
        }
    }

    // Resample the points from the integrator by interpolating at equal
    // intervals.
    return calcResampledGeodesicPoints(
        geodesic_S.samples,
        X_GS,
        nSamples,
        points);
}

//==============================================================================
//                                ???
//==============================================================================

namespace
{

const Correction* calcPathCorrections(SolverData& data)
{
    Real w = data.pathError.normInf();

    data.mat = data.pathErrorJacobian.transpose() * data.pathErrorJacobian;
    for (int i = 0; i < data.mat.nrow(); ++i) {
        data.mat[i][i] += w + 1e-3;
    }
    data.matInv = data.mat;
    data.vec    = data.pathErrorJacobian.transpose() * (data.pathError * (-1.));
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
//                                SURFACE IMPL
//==============================================================================

void CurveSegment::Impl::calcLiftoffIfNeeded(
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S,
    CurveSegment::Impl::InstanceEntry& cache) const
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
    if (dot(prevPoint_S - g.K_P.p(), g.K_P.R().getAxisUnitVec(NormalAxis)) <=
            0. ||
        dot(nextPoint_S - g.K_P.p(), g.K_P.R().getAxisUnitVec(NormalAxis)) <=
            0.) {
        // No Lifted.
        return;
    }

    // Liftoff detected: update status.
    g.status = Status::Lifted;
    // Initialize the tracking point from the last geodesic start point.
    cache.trackingPointOnLine = g.K_P.p();
}

void CurveSegment::Impl::calcTouchdownIfNeeded(
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S,
    CurveSegment::Impl::InstanceEntry& cache) const
{
    // Only attempt touchdown when lifted.
    LocalGeodesicInfo& g = cache;
    if (g.status != Status::Lifted) {
        return;
    }

    // Detect touchdown by computing the point on the line from x_QS to x_PS
    // that is nearest to the surface.
    const bool touchdownDetected = calcNearestPointOnLineImplicitly(
        m_Geometry,
        prevPoint_S,
        nextPoint_S,
        cache.trackingPointOnLine,
        m_TouchdownIter,
        m_TouchdownAccuracy);
    if (!touchdownDetected) {
        return;
    }

    // Touchdown detected: Remove the Lifted status flag.
    g.status = Status::Ok;
    // Shoot a zero length geodesic at the touchdown point.
    shootNewGeodesic(
        cache.trackingPointOnLine,
        nextPoint_S - prevPoint_S,
        0.,
        cache.sHint,
        cache);
}

void CurveSegment::Impl::assertSurfaceBounds(
    const Vec3& prevPoint_S,
    const Vec3& nextPoint_S) const
{
    // Make sure that the previous point does not lie inside the surface.
    if (calcSurfaceConstraintValue(m_Geometry, prevPoint_S) < 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("Unable to wrap over surface: Preceding point "
                                 "lies inside the surface");
    }
    if (calcSurfaceConstraintValue(m_Geometry, nextPoint_S) < 0.) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error(
            "Unable to wrap over surface: Next point lies inside the surface");
    }
}

void CurveSegment::Impl::shootNewGeodesic(
    Vec3 x,
    Vec3 t,
    Real l,
    Real sHint,
    InstanceEntry& cache) const
{
    cache.samples.clear();
    calcGeodesicAndVariationImplicitly(
        m_Geometry,
        x,
        t,
        l,
        sHint,
        cache.K_P,
        cache.dK_P,
        cache.K_Q,
        cache.dK_Q,
        m_IntegratorAccuracy,
        m_ProjectionMaxIter,
        m_ProjectionAccuracy,
        cache.samples);
    cache.length = l;
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
    m_Subsystem(&path.updImpl().updSubsystem()),
    m_Path(path), m_Index(-1), // TODO what to do with this index, and when
    m_Mobod(mobod), m_X_BS(X_BS), m_Geometry(geometry),
    m_ContactPointHint_S(initPointGuess),
    m_Decoration(geometry.createDecorativeGeometry()
                     .setColor(Orange)
                     .setOpacity(.75)
                     .setResolution(3))
{}

void CurveSegment::Impl::realizeCablePosition(const State& s) const
{
    m_Path.getImpl().realizePosition(s);
}

//==============================================================================
//                                PATH HELPERS
//==============================================================================

namespace
{

static const int N_PATH_CONSTRAINTS = 4;

void addDirectionJacobian(
    const LineSegment& e,
    const UnitVec3& axis,
    const PointVariation& dx,
    std::function<void(const Vec4&)>& AddBlock,
    bool invert = false)
{
    Vec3 y = axis - e.d * dot(e.d, axis);
    y /= e.l * (invert ? -1. : 1);
    AddBlock((~dx) * y);
}

Real calcPathError(const LineSegment& e, const Rotation& R, CoordinateAxis axis)
{
    return dot(e.d, R.getAxisUnitVec(axis));
}

void addPathErrorJacobian(
    const LineSegment& e,
    const UnitVec3& axis,
    const Variation& dK,
    std::function<void(const Vec4&)>& AddBlock,
    bool invertV = false)
{
    addDirectionJacobian(e, axis, dK[1], AddBlock, invertV);
    AddBlock((~dK[0]) * cross(axis, e.d));
}

} // namespace

//==============================================================================
//                              CABLE SPAN IMPL
//==============================================================================

void CableSpan::Impl::realizeTopology(State& s)
{
    for (CurveSegment segment : m_CurveSegments) {
        segment.updImpl().realizeTopology(s);
    }

    PosInfo posInfo{};
    m_PosInfoIx = updSubsystem().allocateCacheEntry(
        s,
        Stage::Position,
        Stage::Infinity,
        new Value<PosInfo>(posInfo));

    getSubsystem().markCacheValueNotRealized(s, m_PosInfoIx);

    VelInfo velInfo{};
    m_VelInfoIx = updSubsystem().allocateCacheEntry(
        s,
        Stage::Velocity,
        Stage::Infinity,
        new Value<VelInfo>(velInfo));

    getSubsystem().markCacheValueNotRealized(s, m_VelInfoIx);
}

void CableSpan::Impl::realizePosition(const State& s) const
{
    if (getSubsystem().isCacheValueRealized(s, m_PosInfoIx)) {
        return;
    }
    calcPosInfo(s, updPosInfo(s));
    getSubsystem().markCacheValueRealized(s, m_PosInfoIx);
}

void CableSpan::Impl::realizeVelocity(const State& s) const
{
    realizePosition(s);
    if (getSubsystem().isCacheValueRealized(s, m_VelInfoIx)) {
        return;
    }
    calcVelInfo(s, updVelInfo(s));
    getSubsystem().markCacheValueRealized(s, m_VelInfoIx);
}

const CableSpan::Impl::PosInfo& CableSpan::Impl::getPosInfo(
    const State& s) const
{
    realizePosition(s);
    return Value<PosInfo>::downcast(
        getSubsystem().getCacheEntry(s, m_PosInfoIx));
}

CableSpan::Impl::PosInfo& CableSpan::Impl::updPosInfo(const State& s) const
{
    return Value<PosInfo>::updDowncast(
        getSubsystem().updCacheEntry(s, m_PosInfoIx));
}

const CableSpan::Impl::VelInfo& CableSpan::Impl::getVelInfo(
    const State& s) const
{
    realizeVelocity(s);
    return Value<VelInfo>::downcast(
        getSubsystem().getCacheEntry(s, m_VelInfoIx));
}

CableSpan::Impl::VelInfo& CableSpan::Impl::updVelInfo(const State& s) const
{
    return Value<VelInfo>::updDowncast(
        getSubsystem().updCacheEntry(s, m_VelInfoIx));
}

const CurveSegment* CableSpan::Impl::findPrevActiveCurveSegment(
    const State& s,
    CurveSegmentIndex ix) const
{
    for (int i = ix - 1; i >= 0; --i) {
        // Find the active segment before the current.
        if (m_CurveSegments.at(CurveSegmentIndex(i))
                .getImpl()
                .getInstanceEntry(s)
                .isActive()) {
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
        if (m_CurveSegments.at(CurveSegmentIndex(i))
                .getImpl()
                .getInstanceEntry(s)
                .isActive()) {
            return &m_CurveSegments.at(CurveSegmentIndex(i));
        }
    }
    return nullptr;
}

Vec3 CableSpan::Impl::findPrevPoint(const State& s, CurveSegmentIndex ix) const
{
    const CurveSegment* segment = findPrevActiveCurveSegment(s, ix);
    return segment ? segment->getImpl().calcFinalContactPoint(s)
                   : m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(
                         m_OriginPoint);
}

Vec3 CableSpan::Impl::findNextPoint(const State& s, CurveSegmentIndex ix) const
{
    const CurveSegment* segment = findNextActiveCurveSegment(s, ix);
    return segment
               ? segment->getImpl().calcInitialContactPoint(s)
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
    pathError *= 0;

    for (const CurveSegment& segment : m_CurveSegments) {
        if (!segment.getImpl().getInstanceEntry(s).isActive()) {
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
    J *= 0.;

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
        if (!segment.getImpl().getInstanceEntry(s).isActive()) {
            continue;
        }
        const GeodesicInfo& g = segment.getImpl().getPosInfo(s);

        const LineSegment& l_P = lines.at(activeIx);
        const LineSegment& l_Q = lines.at(++activeIx);

        const CurveSegmentIndex ix = segment.getImpl().getIndex();
        const CurveSegment* prev   = findPrevActiveCurveSegment(s, ix);
        const CurveSegment* next   = findNextActiveCurveSegment(s, ix);

        int blkCol                                = col;
        std::function<void(const Vec4&)> AddBlock = [&](const Vec4& block) {
            for (int ix = 0; ix < 4; ++ix) {
                J.row(row).col(blkCol + ix) += block[ix];
            }
        };

        const Variation& dK_P = g.dKP;
        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_P = g.KP.R().getAxisUnitVec(axis);

            blkCol = col;
            addPathErrorJacobian(l_P, a_P, dK_P, AddBlock, false);

            if (prev) {
                const Variation& prev_dK_Q = prev->getImpl().getPosInfo(s).dKQ;
                blkCol                     = col - Nq;
                addDirectionJacobian(l_P, a_P, prev_dK_Q[1], AddBlock, true);
            }
            ++row;
        }

        const Variation& dK_Q = g.dKQ;
        for (CoordinateAxis axis : axes) {
            const UnitVec3 a_Q = g.KQ.R().getAxisUnitVec(axis);

            blkCol = col;
            addPathErrorJacobian(l_Q, a_Q, dK_Q, AddBlock, true);

            if (next) {
                const Variation& next_dK_P = next->getImpl().getPosInfo(s).dKP;
                blkCol                     = col + Nq;
                addDirectionJacobian(l_Q, a_Q, next_dK_P[1], AddBlock, false);
            }
            ++row;
        }

        col += Nq;
    };
}

Real CableSpan::Impl::calcCableLength(
    const State& s,
    const std::vector<LineSegment>& lines) const
{
    Real lTot = 0.;
    for (const LineSegment& line : lines) {
        // TODO spell out as length.
        lTot += line.l;
    }

    for (const CurveSegment& segment : m_CurveSegments) {
        if (!segment.getImpl().getInstanceEntry(s).isActive()) {
            continue;
        }
        lTot += segment.getImpl().getInstanceEntry(s).length;
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
        if (!segment.getImpl().getInstanceEntry(s).isActive()) {
            continue;
        }

        const GeodesicInfo& g = segment.getImpl().getPosInfo(s);
        const Vec3 lineEnd    = g.KP.p();
        lines.push_back(LineSegment(lineStart, lineEnd));

        lineStart = g.KQ.p();
    }
    lines.emplace_back(lineStart, p_I);
}

size_t CableSpan::countActive(const State& s) const
{
    return getImpl().countActive(s);
}

size_t CableSpan::Impl::countActive(const State& s) const
{
    size_t count = 0;
    for (const CurveSegment& segment : m_CurveSegments) {
        if (segment.getImpl().getInstanceEntry(s).isActive()) {
            ++count;
        }
    }
    return count;
}

/* template <size_t N> */
void CableSpan::calcPathErrorJacobian(
    const State& s,
    /* const std::array<CoordinateAxis, N>& axes, */
    Vector& e,
    Matrix& J) const
{
    getImpl().calcPathErrorJacobianUtility(s, e, J);
}

void CableSpan::applyCorrection(const State& s, const Vector& c) const
{
    getImpl().applyCorrection(s, c);
}

void CableSpan::Impl::applyCorrection(const State& s, const Vector& c) const
{
    const size_t nActive = countActive(s);
    if (nActive * GeodesicDOF != c.size()) {
        // TODO use SimTK_ASSERT
        throw std::runtime_error("invalid size of corrections vector");
    }

    const Correction* corrIt = reinterpret_cast<const Correction*>(&c[0]);
    for (const CurveSegment& curve : m_CurveSegments) {
        if (!curve.getImpl().getInstanceEntry(s).isActive()) {
            continue;
        }
        curve.getImpl().applyGeodesicCorrection(s, *corrIt);
        ++corrIt;
    }

    // Path has changed: invalidate each segment's cache.
    invalidatePositionLevelCache(s);
}

void CableSpan::Impl::invalidatePositionLevelCache(const State& s) const
{
    for (const CurveSegment& curve : m_CurveSegments) {
        curve.getImpl().invalidatePosEntry(s);
    }
    getSubsystem().markCacheValueNotRealized(s, m_PosInfoIx);
}

/* template<size_t N> */
void CableSpan::Impl::calcPathErrorJacobianUtility(
    const State& s,
    /* const std::array<CoordinateAxis, N>& axes, */
    Vector& e,
    Matrix& J) const
{
    const std::array<CoordinateAxis, 2> axes = {NormalAxis, BinormalAxis};
    // Force re-realizing the cache.
    invalidatePositionLevelCache(s);
    for (const CurveSegment& curve : m_CurveSegments) {
        curve.getImpl().realizePosition(
            s,
            findPrevPoint(s, curve),
            findNextPoint(s, curve));
    }

    const size_t nActive = countActive(s);

    if (nActive == 0) {
        J.resize(0, 0);
        return;
    }

    // Grab the shared data cache for computing the matrices, and lock it.
    SolverData& data =
        getSubsystem().getImpl().updCachedScratchboard(s).updOrInsert(nActive);

    // Compute the straight-line segments.
    const Vec3 x_O =
        m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
    const Vec3 x_I =
        m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(
            m_TerminationPoint);
    calcLineSegments(s, x_O, x_I, data.lineSegments);

    // Evaluate path error.
    calcPathErrorVector<2>(s, data.lineSegments, axes, data.pathError);

    // Evaluate the path error jacobian.
    calcPathErrorJacobian<2>(
        s,
        data.lineSegments,
        axes,
        data.pathErrorJacobian);

    e = data.pathError;
    J = data.pathErrorJacobian;
}

void CableSpan::Impl::calcPosInfo(const State& s, PosInfo& posInfo) const
{
    // Path origin and termination point.
    const Vec3 x_O =
        m_OriginBody.getBodyTransform(s).shiftFrameStationToBase(m_OriginPoint);
    const Vec3 x_I =
        m_TerminationBody.getBodyTransform(s).shiftFrameStationToBase(
            m_TerminationPoint);

    posInfo.xO = x_O;
    posInfo.xI = x_I;

    // Axes considered when computing the path error.
    const std::array<CoordinateAxis, 2> axes{NormalAxis, BinormalAxis};

    posInfo.loopIter = 0;
    while (true) {
        // Make sure all curve segments are realized to position stage.
        // This will transform all last computed geodesics to Ground frame, and
        // will update each curve's Status.
        for (const CurveSegment& curve : m_CurveSegments) {
            curve.getImpl().realizePosition(
                s,
                findPrevPoint(s, curve),
                findNextPoint(s, curve));
        }

        // Count the number of active curve segments.
        const size_t nActive = countActive(s);

        // If the path contains no curved segments it is a straight line.
        if (nActive == 0) {
            // Update the path length, and exit.
            posInfo.l = (x_I - x_O).norm();
            break;
        }

        // Grab the shared data cache for helping with computing the path corrections.
        // This data is only used as an intermediate variable, and will be
        // discarded after each iteration.
        SolverData& data =
            getSubsystem().getImpl().updCachedScratchboard(s).updOrInsert(
                nActive);

        // Compute the straight-line segments of this cable span. Note that
        // there is one more straight line segment, than there are active curve
        // segments.
        calcLineSegments(s, x_O, x_I, data.lineSegments);

        // Evaluate path error as the misalignment of the straight line
        // segments with the curve segment's tangent vectors at the contact
        // points.
        calcPathErrorVector<2>(s, data.lineSegments, axes, data.pathError);
        const Real maxPathError = data.pathError.normInf();

        // Stop iterating if max path error is small, or max iterations is reached.
        if (maxPathError < m_PathErrorBound || posInfo.loopIter >= m_PathMaxIter) {
            posInfo.l = calcCableLength(s, data.lineSegments);
            break;
        }

        // Evaluate the path error jacobian to the natural geodesic corrections
        // of each curve segment.
        calcPathErrorJacobian<2>(
            s,
            data.lineSegments,
            axes,
            data.pathErrorJacobian);

        // Compute the geodesic corrections for each curve segment.
        const Correction* corrIt = calcPathCorrections(data);

        // Apply corrections to the curve segments.
        for (const CurveSegment& curve : m_CurveSegments) {
            // The corrections were computed for the active segments: Skip any
            // non-active here.
            if (!curve.getImpl().getInstanceEntry(s).isActive()) {
                continue;
            }
            curve.getImpl().applyGeodesicCorrection(s, *corrIt);
            ++corrIt;
        }

        // Applying the corrections changes the path: invalidate each segment's
        // cache.
        for (const CurveSegment& curve : m_CurveSegments) {
            // Also invalidate non-active segments: They might touchdown again.
            curve.getImpl().invalidatePosEntry(s);
        }

         ++posInfo.loopIter;
    }

    // TODO throw?
    if (posInfo.loopIter >= m_PathMaxIter) {
        std::cout << "CableSpan::calcPosInfo(): Reached max iterations!\n";
    }
}

void CableSpan::Impl::calcVelInfo(const State& s, VelInfo& velInfo) const
{
    auto CalcPointVelocityInGround = [&](const MobilizedBody& mobod,
                                         const Vec3& point_G) -> Vec3 {
        // Not using MobilizedBody::findStationVelocityInGround because the
        // point_G is in ground frame. The following computation is the same
        // though (minus transforming the point to the ground frame).

        // Get body kinematics in ground frame.
        const Vec3& x_BG = mobod.getBodyOriginLocation(s);
        const Vec3& w_BG = mobod.getBodyAngularVelocity(s);
        const Vec3& v_BG = mobod.getBodyOriginVelocity(s);

        // Compute surface point velocity in ground frame.
        return v_BG + w_BG % (point_G - x_BG);
    };

    Real& lengthDot = (velInfo.lengthDot = 0.);

    Vec3 v_GQ = m_OriginBody.findStationVelocityInGround(s, m_OriginPoint);
    const CurveSegment* lastActive = nullptr;
    for (const CurveSegment& curve : m_CurveSegments) {
        if (!curve.getImpl().getInstanceEntry(s).isActive()) {
            continue;
        }

        const MobilizedBody& mobod = curve.getImpl().getMobilizedBody();
        // TODO odd name: "g"
        const CurveSegment::Impl::PosInfo& g = curve.getImpl().getPosInfo(s);
        const UnitVec3 e_G = g.KP.R().getAxisUnitVec(TangentAxis);

        const Vec3 v_GP = CalcPointVelocityInGround(mobod, g.KP.p());

        lengthDot += dot(e_G, v_GP - v_GQ);

        v_GQ = CalcPointVelocityInGround(mobod, g.KQ.p());

        lastActive = &curve;
    }

    const Vec3 v_GP =
        m_TerminationBody.findStationVelocityInGround(s, m_TerminationPoint);

    const PosInfo& pos = getPosInfo(s);
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
    if (tension < 0.) {
        // TODO throw? or skip?
        throw std::runtime_error("Cable tension should be nonnegative.");
    }

    realizePosition(s);
    SpatialVec unitForce_G;

    {
        calcUnitForceAtOrigin(s, unitForce_G);
        getOriginBody().applyBodyForce(s, unitForce_G * tension, bodyForcesInG);
    }

    for (const CurveSegment& curve : m_CurveSegments) {
        if (!curve.isActive(s)) {
            return;
        }

        curve.calcUnitForce(s, unitForce_G);
        curve.getMobilizedBody().applyBodyForce(
            s,
            unitForce_G * tension,
            bodyForcesInG);
    }

    {
        calcUnitForceAtTermination(s, unitForce_G);
        getTerminationBody().applyBodyForce(
            s,
            unitForce_G * tension,
            bodyForcesInG);
    }
}

int CableSpan::Impl::calcDecorativeGeometryAndAppend(
    const State& s,
    Stage stage,
    Array_<DecorativeGeometry>& decorations) const
{
    const Vec3 color = Green; // Red Purple

    const PosInfo& ppe = getPosInfo(s);

    // Draw point at origin and termination.
    decorations.push_back(DecorativePoint(ppe.xO).setColor(Green));
    decorations.push_back(DecorativePoint(ppe.xI).setColor(Red));

    if (countActive(s) == 0) {
        decorations.push_back(DecorativeLine(ppe.xO, ppe.xI)
                .setColor(Purple)
                .setLineThickness(3));
    }

    for (const CurveSegment& curveSegment : m_CurveSegments) {
        const CurveSegment::Impl& curve = curveSegment.getImpl();

        const Transform X_GS   = curve.calcSurfaceFrameInGround(s);
        DecorativeGeometry geo = curve.getDecoration();
        const Transform& X_SD  = geo.getTransform(); // TODO getTransform?

        // Inactive surfaces are dimmed.
        if (!curve.getInstanceEntry(s).isActive()) {
            decorations.push_back(geo.setTransform(X_GS * X_SD)
                                      .setColor(Real(0.75) * geo.getColor()));
            continue;
        }

        decorations.push_back(geo.setTransform(X_GS * X_SD));

        {
            const Vec3 prevPoint = findPrevPoint(s, curveSegment);
            const Vec3 x_P       = curve.getPosInfo(s).KP.p();

            decorations.push_back(DecorativeLine(prevPoint, x_P)
                                      .setColor(Orange)
                                      .setLineThickness(2));
        }

        Transform K_P                            = curve.getPosInfo(s).KP;
        Transform K_Q                            = curve.getPosInfo(s).KQ;
        const std::array<CoordinateAxis, 3> axes = {
            TangentAxis,
            NormalAxis,
            BinormalAxis};
        const std::array<Vec3, 3> colors = {Red, Green, Blue};

        for (size_t i = 0; i < 3; ++i) {
            decorations.push_back(
                DecorativeLine(
                    K_P.p(),
                    K_P.p() + 0.5 * K_P.R().getAxisUnitVec(axes.at(i)))
                    .setColor(colors.at(i))
                    .setLineThickness(4));
            decorations.push_back(
                DecorativeLine(
                    K_Q.p(),
                    K_Q.p() + 0.5 * K_Q.R().getAxisUnitVec(axes.at(i)))
                    .setColor(colors.at(i))
                    .setLineThickness(4));
        }

        {
            const Vec3 nextPoint = findNextPoint(s, curveSegment);
            const Vec3 x_Q       = curve.getPosInfo(s).KQ.p();

            decorations.push_back(DecorativeLine(nextPoint, x_Q)
                                      .setColor(Gray)
                                      .setLineThickness(2));
        }

        curveSegment.getImpl().calcDecorativeGeometryAndAppend(s, decorations);
    }

    /* decorations.push_back(DecorativeLine(lastCurvePoint, ppe.xI) */
    /*         .setColor(Purple) */
    /*         .setLineThickness(3)); */
    return 0;
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
