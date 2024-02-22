#include <config.h>

#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

namespace Opm {

INSTANCE_PF(double,5)

#if FLOW_INSTANCE_FLOAT
INSTANCE_PF(float,5)
#endif

}
