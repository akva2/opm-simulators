/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED
#define OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED

#include <functional>
#include <vector>

namespace Opm
{

class DeferredLogger;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class StandardWellConnections
{
public:
    StandardWellConnections(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well);

    void computeConnectionPressureDelta();

    // TODO: not total sure whether it is a good idea to put this function here
    // the major reason to put here is to avoid the usage of Wells struct
    void computeConnectionDensities(const std::vector<Scalar>& perfComponentRates,
                                    const std::vector<Scalar>& b_perf,
                                    const std::vector<Scalar>& rsmax_perf,
                                    const std::vector<Scalar>& rvmax_perf,
                                    const std::vector<Scalar>& rvwmax_perf,
                                    const std::vector<Scalar>& surf_dens_perf,
                                    DeferredLogger& deferred_logger);

    void computePropertiesForWellConnectionPressures(const WellState& well_state,
                                                     const std::function<Scalar(int,int)>& getTemperature,
                                                     const std::function<Scalar(int)>& getSaltConcentration,
                                                     const std::function<int(int)>& pvtRegionIdx,
                                                     const std::function<Scalar(int)>& solventInverseFormationVolumeFactor,
                                                     const std::function<Scalar(int)>& solventRefDensity,
                                                     std::vector<Scalar>& b_perf,
                                                     std::vector<Scalar>& rsmax_perf,
                                                     std::vector<Scalar>& rvmax_perf,
                                                     std::vector<Scalar>& rvwmax_perf,
                                                     std::vector<Scalar>& surf_dens_perf) const;

    void computeWellConnectionDensitesPressures(const WellState& well_state,
                                                const std::function<Scalar(int,int)>& invB,
                                                const std::function<Scalar(int,int)>& mobility,
                                                const std::function<Scalar(int)>& solventInverseFormationVolumeFactor,
                                                const std::function<Scalar(int)>& solventMobility,
                                                const std::vector<double>& b_perf,
                                                const std::vector<double>& rsmax_perf,
                                                const std::vector<double>& rvmax_perf,
                                                const std::vector<double>& rvwmax_perf,
                                                const std::vector<double>& surf_dens_perf,
                                                DeferredLogger& deferred_logger);

    Scalar getRho() const
    {
        return this->perf_densities_.empty() ? 0.0 : perf_densities_[0];
    }

    // densities of the fluid in each perforation
    std::vector<Scalar> perf_densities_;
    // pressure drop between different perforations
    std::vector<Scalar> perf_pressure_diffs_;

private:
    // Base interface reference
    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_;
};

}

#endif // OPM_STANDARDWELL_CONNECTIONS_HEADER_INCLUDED
