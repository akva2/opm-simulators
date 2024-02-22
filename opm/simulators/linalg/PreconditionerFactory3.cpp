#include <config.h>

#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

namespace Opm {

INSTANCE_PF(double,3)

#if FLOW_INSTANCE_FLOAT
INSTANCE_PF(float,3)
#endif

}
