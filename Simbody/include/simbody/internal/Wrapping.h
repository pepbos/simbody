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

    enum class Status
    {
        Ok,
        Liftoff,
        Disabled,
    };

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

    const Impl& getImpl() const { return *impl; }
    Impl& updImpl() { return *impl; }

    std::shared_ptr<Impl> impl = nullptr;

    friend WrapObstacle::Impl;
    friend WrappingPathSubsystem;
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
