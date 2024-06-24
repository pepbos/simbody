/*-----------------------------------------------------------------------------
                Simbody(tm) Example: Cable Over Smooth Surfaces
-------------------------------------------------------------------------------
 Copyright (c) 2024 Authors.
 Authors: Pepijn van den Bos
 Contributors:

 Licensed under the Apache License, Version 2.0 (the "License"); you may
 not use this file except in compliance with the License. You may obtain a
 copy of the License at http://www.apache.org/licenses/LICENSE-2.0.

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 ----------------------------------------------------------------------------*/

/* This example is for experimenting with a CableSpan over obstacle surfaces.
The cable endpoints are manually repositioned to control the next path
solving problem. */

#include "Simbody.h"
#include "simbody/internal/CableSpan.h"
#include <iostream>

using namespace SimTK;

// A helper class for printing interesting cable outputs.
class ShowStuff : public PeriodicEventReporter {
public:
    ShowStuff(const CableSubsystem& cables, Real interval) :
        PeriodicEventReporter(interval), m_cables(cables)
    {}

    static void showHeading(std::ostream& o)
    {
        printf(
            "%10s %10s %10s %10s",
            "Cable index",
            "iterations",
            "path error",
            "length");
    }

    /** This is the implementation of the EventReporter virtual. **/
    void handleEvent(const State& s) const override
    {
        for (CableSpanIndex cableIx(0); cableIx < m_cables.getNumCables();
             ++cableIx) {
            const CableSpan& cable = m_cables.getCable(cableIx);
            printf(
                "%d %d %g %g",
                static_cast<int>(cableIx),
                cable.getNumSolverIter(s),
                cable.getMaxPathError(s),
                cable.getLength(s));
        }
    }

private:
    const CableSubsystem& m_cables;
};

/* A helper class for drawing interesting things of a CableSpan. */
class CableDecorator : public SimTK::DecorationGenerator {
public:
    CableDecorator(MultibodySystem& mbs, const CableSpan& cable) :
        m_mbs(&mbs), m_cable(cable)
    {
        for (CableSpanObstacleIndex ix(0);
             ix < m_cable.getNumSurfaceObstacles();
             ++ix) {
            m_obstacleDecorations.push_back(
                m_cable.getObstacleContactGeometry(ix)
                    .createDecorativeGeometry()
                    .setResolution(3));
        }
    }

    virtual void generateDecorations(
        const State& state,
        Array_<DecorativeGeometry>& decorations) override
    {
        for (CableSpanObstacleIndex ix(0);
             ix < m_cable.getNumSurfaceObstacles();
             ++ix) {
            // Draw the obstacle surface.
            // Green if cable is in contact with surface, grey otheriwse.
            const ContactGeometry& geometry =
                m_cable.getObstacleContactGeometry(ix);
            const bool isInContactWithSurface =
                m_cable.isInContactWithObstacle(state, ix);
            const Vec3 color = isInContactWithSurface ? Yellow : Gray;
            const Real opacity = isInContactWithSurface ? 1. : 0.25;
            Transform X_GB =
                m_mbs->getMatterSubsystem()
                    .getMobilizedBody(m_cable.getObstacleMobilizedBodyIndex(ix))
                    .getBodyTransform(state);
            Transform X_GS =
                X_GB.compose(m_cable.getObstacleXformSurfaceToBody(ix));
            decorations.push_back(m_obstacleDecorations.at(ix)
                                      .setTransform(X_GS)
                                      .setColor(color)
                                      .setOpacity(opacity));

            // Draw the initial contact point hints (these are user-defined) as
            // an orange line with a point.
            const Vec3 x_PS = m_cable.getObstacleInitialContactPointHint(ix);
            decorations.push_back(
                DecorativeLine(X_GS.p(), X_GS.shiftFrameStationToBase(x_PS))
                    .setColor(Orange)
                    .setLineThickness(3));
            decorations.push_back(
                DecorativePoint(X_GS.shiftFrameStationToBase(x_PS))
                    .setColor(Orange));

            if (!isInContactWithSurface) {
                continue;
            }

            // Draw the frenet frames at the geodesic boundary points.

            // Get the geodesic knots at the boundary points.
            bool isInitialKnot = true;
            std::array<ContactGeometry::GeodesicKnotPoint, 2>
                geodesicBoundaryKnots;
            m_cable.calcCurveSegmentKnots(
                state,
                ix,
                [&](const ContactGeometry::GeodesicKnotPoint& q,
                    const Transform&)
                {
                    geodesicBoundaryKnots.at(isInitialKnot ? 0 : 1) = q;
                    isInitialKnot                                   = false;
                });
            // Draw the frenet frames at the boundary points.
            for (const ContactGeometry::GeodesicKnotPoint& q :
                 geodesicBoundaryKnots) {
                // Calculate the frenet frame with tangent along x-axis and
                // normal along y-axis:
                const Vec3 point_G = X_GS.shiftFrameStationToBase(q.point);
                const UnitVec3 normal_G(X_GS.xformBaseVecToFrame(
                    geometry.calcSurfaceUnitNormal(q.point)));
                const UnitVec3 tangent_G(X_GS.xformBaseVecToFrame(q.tangent));
                Transform frenetFrame;
                frenetFrame.updP() = point_G;
                frenetFrame.updR().setRotationFromTwoAxes(
                    normal_G,
                    CoordinateAxis::YCoordinateAxis(),
                    tangent_G,
                    CoordinateAxis::YCoordinateAxis());
                // Add to list of decorations.
                decorations.push_back(
                    DecorativeFrame(0.2).setTransform(frenetFrame));
            };
        }
    }

    MultibodySystem* m_mbs;
    CableSpan m_cable;
    int m_lineThickness = 3;

    Array_<DecorativeGeometry, CableSpanObstacleIndex> m_obstacleDecorations;
    Array_<Transform, CableSpanObstacleIndex> m_obstacleDecorationsOffsets;
};

int main()
{
    try {
        // Create the system.
        MultibodySystem system;
        SimbodyMatterSubsystem matter(system);
        CableSubsystem cables(system);

        // A dummy body.
        Body::Rigid aBody(MassProperties(1.0, Vec3(0), Inertia(1)));

        // Mobilizer for path origin.
        MobilizedBody::Translation cableOriginBody(
            matter.Ground(),
            Vec3(-8., 0.1, 0.),
            aBody,
            Transform());

        // Mobilizer for path termination.
        MobilizedBody::Translation cableTerminationBody(
            matter.Ground(),
            Transform(Vec3(5., 1.0, -1.)),
            aBody,
            Transform());

        // Construct a new cable.
        CableSpan cable(
            cables,
            cableOriginBody,
            Vec3{0.},
            cableTerminationBody,
            Vec3{0.});

        // Add some obstacles to the cable.

        // Add torus obstacle.
        cable.addSurfaceObstacle(
            matter.Ground(), // Obstacle mobilizer body.
            Transform(
                Rotation(0.5 * Pi, YAxis),
                Vec3{-4., 0., 0.}), // Surface to body transform.
            std::shared_ptr<ContactGeometry>(
                new ContactGeometry::Torus(1., 0.2)), // Obstacle geometry.
            {0.1, 0.2, 0.} // Initial contact point hint.
        );

        // Add ellipsoid obstacle.
        cable.addSurfaceObstacle(
            matter.Ground(),
            Transform(Vec3{-2., 0., 0.}),
            std::shared_ptr<ContactGeometry>(
                new ContactGeometry::Ellipsoid({1.5, 2.6, 1.})),
            {0.0, 0., 1.1});

        // Add sphere obstacle.
        cable.addSurfaceObstacle(
            matter.Ground(),
            Transform(Vec3{2., 0., 0.}),
            std::shared_ptr<ContactGeometry>(new ContactGeometry::Sphere(1.)),
            {0.1, 1.1, 0.}
        );

        // Visaulize the system.
        system.setUseUniformBackground(true); // no ground plane in display
        Visualizer viz(system);
        viz.addDecorationGenerator(new CableDecorator(system, cable));
        ShowStuff showStuff(cables, 1e-3);

        // Initialize the system and s.
        system.realizeTopology();
        State s = system.getDefaultState();
        system.realize(s, Stage::Position);

        system.realize(s, Stage::Report);
        viz.report(s);
        showStuff.handleEvent(s);

        std::cout << "Hit ENTER ..., or q\n";
        const char ch = getchar();
        if (ch == 'Q' || ch == 'q') {
            return 0;
        }

        Real angle = 0.;
        while (true) {
            system.realize(s, Stage::Position);
            viz.report(s);
            cable.storeCurrentPath(s);

            // Move the cable origin.
            angle += 0.01;
            cableOriginBody.setQ(
                s,
                Vec3(sin(angle), 5. * sin(angle * 1.5), 5. * sin(angle * 2.)));
        }

    } catch (const std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << "\n";
    }
}
