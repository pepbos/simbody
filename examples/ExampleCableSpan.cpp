#include "Simbody.h"
#include "simbody/internal/Wrapping.h"
#include <cassert>
#include <iostream>
using std::cout;
using std::endl;

using namespace SimTK;

// This force element implements an elastic cable of a given nominal length,
// and a stiffness k that generates a k*x force opposing stretch beyond
// the slack length. There is also a dissipation term (k*x)*c*xdot. We keep
// track of dissipated power here so we can use conservation of energy to check
// that the cable and force element aren't obviously broken.
class MyCableSpringImpl : public Force::Custom::Implementation
{
public:
    MyCableSpringImpl(
        const GeneralForceSubsystem& forces,
        const CableSpan& cable,
        Real stiffness,
        Real nominal,
        Real damping) :
        forces(forces),
        m_Cable(cable), k(stiffness), x0(nominal), c(damping)
    {
        assert(stiffness >= 0 && nominal >= 0 && damping >= 0);
    }

    const CableSpan& getCable() const
    {
        return m_Cable;
    }

    // Must be at stage Velocity. Evalutes tension if necessary.
    Real getTension(const State& s) const
    {
        ensureTensionCalculated(s);
        return Value<Real>::downcast(forces.getCacheEntry(s, tensionx));
    }

    // Must be at stage Velocity.
    Real getPowerDissipation(const State& s) const
    {
        const Real stretch = calcStretch(s);
        if (stretch == 0)
            return 0;
        const Real rate = m_Cable.getLengthDot(s);
        return k * stretch * std::max(c * rate, -1.) * rate;
    }

    // This integral is always available.
    Real getDissipatedEnergy(const State& s) const
    {
        return forces.getZ(s)[workx];
    }

    //--------------------------------------------------------------------------
    //                       Custom force virtuals

    // Ask the cable to apply body forces given the tension calculated here.
    void calcForce(
        const State& s,
        Vector_<SpatialVec>& bodyForces,
        Vector_<Vec3>& particleForces,
        Vector& mobilityForces) const override
    {
        m_Cable.applyBodyForces(s, getTension(s), bodyForces);
    }

    // Return the potential energy currently stored by the stretch of the cable.
    Real calcPotentialEnergy(const State& s) const override
    {
        const Real stretch = calcStretch(s);
        if (stretch == 0)
            return 0;
        return k * square(stretch) / 2;
    }

    // Allocate the s variable for tracking dissipated energy, and a
    // cache entry to hold the calculated tension.
    void realizeTopology(State& s) const override
    {
        Vector initWork(1, 0.);
        workx    = forces.allocateZ(s, initWork);
        tensionx = forces.allocateLazyCacheEntry(
            s,
            Stage::Velocity,
            new Value<Real>(NaN));
    }

    // Report power dissipation as the derivative for the work variable.
    void realizeAcceleration(const State& s) const override
    {
        Real& workDot = forces.updZDot(s)[workx];
        workDot       = getPowerDissipation(s);
    }
    //--------------------------------------------------------------------------

private:
    // Return the amount by which the cable is stretched beyond its nominal
    // length or zero if the cable is slack. Must be at stage Position.
    Real calcStretch(const State& s) const
    {
        const Real stretch = m_Cable.getLength(s) - x0;
        return std::max(stretch, 0.);
    }

    // Must be at stage Velocity to calculate tension.
    Real calcTension(const State& s) const
    {
        const Real stretch = calcStretch(s);
        if (stretch == 0)
            return 0;
        const Real rate = m_Cable.getLengthDot(s);
        if (c * rate < -1)
            cout << "c*rate=" << c * rate << "; limited to -1\n";
        const Real tension = k * stretch * (1 + std::max(c * rate, -1.));
        return tension;
    }

    // If s is at stage Velocity, we can calculate and store tension
    // in the cache if it hasn't already been calculated.
    void ensureTensionCalculated(const State& s) const
    {
        if (forces.isCacheValueRealized(s, tensionx))
            return;
        Value<Real>::updDowncast(forces.updCacheEntry(s, tensionx)) =
            calcTension(s);
        forces.markCacheValueRealized(s, tensionx);
    }

    const GeneralForceSubsystem& forces;
    CableSpan m_Cable;
    Real k, x0, c;
    mutable ZIndex workx;
    mutable CacheEntryIndex tensionx;
};

// A nice handle to hide most of the cable spring implementation. This defines
// a user's API.
class MyCableSpring : public Force::Custom
{
public:
    MyCableSpring(
        GeneralForceSubsystem& forces,
        const CableSpan& cable,
        Real stiffness,
        Real nominal,
        Real damping) :
        Force::Custom(
            forces,
            new MyCableSpringImpl(forces, cable, stiffness, nominal, damping))
    {}

    // Expose some useful methods.
    const CableSpan& getCable() const
    {
        return getImpl().getCable();
    }
    Real getTension(const State& s) const
    {
        return getImpl().getTension(s);
    }
    Real getPowerDissipation(const State& s) const
    {
        return getImpl().getPowerDissipation(s);
    }
    Real getDissipatedEnergy(const State& s) const
    {
        return getImpl().getDissipatedEnergy(s);
    }

private:
    const MyCableSpringImpl& getImpl() const
    {
        return dynamic_cast<const MyCableSpringImpl&>(getImplementation());
    }
};

// This gets called periodically to dump out interesting things about
// the cables and the system as a whole. It also saves states so that we
// can play back at the end.
static Array_<State> saveStates;
class ShowStuff : public PeriodicEventReporter
{
public:
    ShowStuff(
        const MultibodySystem& mbs,
        const MyCableSpring& cable1,
        /* const MyCableSpring& cable2, */
        Real interval) :
        PeriodicEventReporter(interval),
        mbs(mbs), cable1(cable1)
    {}

    static void showHeading(std::ostream& o)
    {
        printf(
            "%8s %10s %10s %10s %10s %10s %10s %10s %10s %12s\n",
            "time",
            "length",
            "rate",
            "integ-rate",
            "unitpow",
            "tension",
            "disswork",
            "KE",
            "PE",
            "KE+PE-W");
    }

    /** This is the implementation of the EventReporter virtual. **/
    void handleEvent(const State& s) const override
    {
        const CableSpan& path1 = cable1.getCable();
        /* const CableSpan& path2 = cable2.getCable(); */
        printf(
            "%8g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g CPU=%g\n",
            s.getTime(),
            path1.getLength(s),
            path1.getLengthDot(s),
            path1.calcCablePower(s, 1), // unit power
            cable1.getTension(s),
            cable1.getDissipatedEnergy(s),
            cpuTime());
        /* printf( */
        /*     "%8s %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g %10.4g " */
        /*     "%12.6g\n", */
        /*     "", */
        /*     path2.getLength(s), */
        /*     path2.getLengthDot(s), */
        /*     path2.calcCablePower(s, 1), // unit power */
        /*     cable2.getTension(s), */
        /*     cable2.getDissipatedEnergy(s), */
        /*     mbs.calcKineticEnergy(s), */
        /*     mbs.calcPotentialEnergy(s), */
        /*     mbs.calcEnergy(s) + cable1.getDissipatedEnergy(s) + */
        /*         cable2.getDissipatedEnergy(s)); */
        saveStates.push_back(s);
    }

private:
    const MultibodySystem& mbs;
    MyCableSpring cable1;
};

int main()
{
    try {
        // Create the system.
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        CableSubsystem cables(system);
        GeneralForceSubsystem forces(system);

        system.setUseUniformBackground(true); // no ground plane in display

        Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
        // Force::GlobalDamper(forces, matter, 5);

        Body::Rigid someBody(MassProperties(1.0, Vec3(0), Inertia(1)));
        const Real Rad = .25;

        Body::Rigid biggerBody(MassProperties(1.0, Vec3(0), Inertia(1)));
        const Real BiggerRad = .5;

        const Vec3 radii(.4, .25, .15);
        Body::Rigid ellipsoidBody(
            MassProperties(1.0, Vec3(0), 1. * UnitInertia::ellipsoid(radii)));

        const Real CylRad = .3, HalfLen = .5;
        Body::Rigid cylinderBody(MassProperties(
            1.0,
            Vec3(0),
            1. * UnitInertia::cylinderAlongX(Rad, HalfLen)));

        Body::Rigid fancyBody = biggerBody; // NOT USING ELLIPSOID

        MobilizedBody Ground = matter.Ground();

        MobilizedBody::Ball body1(
            Ground,
            Transform(Vec3(0)),
            someBody,
            Transform(Vec3(0, 1, 0)));
        MobilizedBody::Ball body2(
            body1,
            Transform(Vec3(0)),
            someBody,
            Transform(Vec3(0, 1, 0)));
        MobilizedBody::Ball body3(
            body2,
            Transform(Vec3(0)),
            someBody,
            Transform(Vec3(0, 1, 0)));
        MobilizedBody::Ball body4(
            body3,
            Transform(Vec3(0)),
            fancyBody,
            Transform(Vec3(0, 1, 0)));
        MobilizedBody::Ball body5(
            body4,
            Transform(Vec3(0)),
            someBody,
            Transform(Vec3(0, 1, 0)));

        CableSpan path1(
            cables,
            body1,
            Vec3(-2*Rad, 0, 0), // origin
            body5,
            Vec3(2*Rad, 0, 0)); // termination

        /* CableObstacle::ViaPoint p1(path1, body2, Rad * UnitVec3(1, 1, 0)); */

        // obs4
        path1.adoptWrappingObstacle(
            body3,
            Transform(),
            ContactGeometry::Sphere(Rad),
            {0., 1., 0.});

        // obs5
        path1.adoptWrappingObstacle(
            body4,
            Transform(),
            ContactGeometry::Sphere(BiggerRad),
            {-1., 0., -2.});

        MyCableSpring cable1(forces, path1, 100., 3.5, 0.1);

        /* CableSpan path2( */
        /*     cables, */
        /*     body3, */
        /*     2 * Rad * UnitVec3(1, 1, 1), */
        /*     Ground, */
        /*     Vec3(-2.5, 1, 0)); */
        /* MyCableSpring cable2(forces, path2, 100., 2, 0.1); */

        // obs1.setPathPreferencePoint(Vec3(2,3,4));
        // obs1.setDecorativeGeometry(DecorativeSphere(0.25).setOpacity(.5));

        Visualizer viz(system);
        viz.setShowFrameNumber(true);
        system.addEventReporter(new Visualizer::Reporter(viz, 0.1 * 1. / 30));
        system.addEventReporter(
            new ShowStuff(system, cable1, 0.1 * 0.1));
        // Initialize the system and s.

        system.realizeTopology();
        State s = system.getDefaultState();
        Random::Gaussian random;
        for (int i = 0; i < s.getNQ(); ++i)
            s.updQ()[i] = random.getValue();
        for (int i = 0; i < s.getNU(); ++i)
            s.updU()[i] = 0.1 * random.getValue();

        system.realize(s, Stage::Position);
        viz.report(s);
        cout << "path1 init length=" << path1.getLength(s) << endl;
        /* cout << "path2 init length=" << path2.getLength(s) << endl; */
        cout << "Hit ENTER ...";
        getchar();

        // Simulate it.
        saveStates.clear();
        saveStates.reserve(2000);

        RungeKuttaMersonIntegrator integ(system);
        // CPodesIntegrator integ(system);
        integ.setAccuracy(1e-3);
        // integ.setAccuracy(1e-6);
        TimeStepper ts(system, integ);
        ts.initialize(s);
        ShowStuff::showHeading(cout);

        const Real finalTime   = 1.;
        const double startTime = realTime();
        ts.stepTo(finalTime);
        cout << "DONE with " << finalTime << "s simulated in "
             << realTime() - startTime << "s elapsed.\n";

        while (true) {
            cout << "Hit ENTER FOR REPLAY, Q to quit ...";
            const char ch = getchar();
            if (ch == 'q' || ch == 'Q')
                break;
            for (unsigned i = 0; i < saveStates.size(); ++i)
                viz.report(saveStates[i]);
        }

    } catch (const std::exception& e) {
        cout << "EXCEPTION: " << e.what() << "\n";
    }
}
