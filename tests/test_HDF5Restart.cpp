/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <config.h>

#define BOOST_TEST_MODULE TestHDF5Restart
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <ebos/ebos.hh>
#include <ebos/hdf5serializer.hh>
#include <ebos/eclproblem.hh>
#include <ebos/ecltracermodel.hh>
#include <ebos/ecltransmissibility.hh>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ASTNode.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>

#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilDefaultIndexTraits.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>

#include <boost/date_time/gregorian/gregorian.hpp>

#define TEST_FOR_TYPE(T) \
    BOOST_AUTO_TEST_CASE(T) \
    { \
        auto data_out = Opm::T::serializationTestObject(); \
        { \
            Opm::HDF5Serializer writer("hdf5_test_" #T ".hdf5", Opm::HDF5File::OpenMode::WRITE); \
            writer.write(data_out, "/0", "data"); \
        } \
        Opm::HDF5Serializer reader("hdf5_test_" #T ".hdf5", Opm::HDF5File::OpenMode::READ); \
        decltype(data_out) data_in; \
        reader.read(data_in, "/0", "data"); \
        BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized " #T " differ"); \
    }

TEST_FOR_TYPE(SimulatorTimer)
TEST_FOR_TYPE(HardcodedTimeStepControl)
TEST_FOR_TYPE(PIDAndIterationCountTimeStepControl)
TEST_FOR_TYPE(PIDTimeStepControl)
TEST_FOR_TYPE(Schedule)
TEST_FOR_TYPE(SimpleIterationCountTimeStepControl)

namespace Opm { using ATE = AdaptiveTimeSteppingEbos<Properties::TTag::EbosTypeTag>; }
TEST_FOR_TYPE(ATE)

BOOST_AUTO_TEST_CASE(EclTransmissibility)
{
    using Transmissibility = Opm::EclTransmissibility<Dune::CpGrid,
                                                      Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                                      Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                                      Dune::CartesianIndexMapper<Dune::CpGrid>,
                                                      double>;
    auto eclState = Opm::EclipseState::serializationTestObject();
    Dune::CpGrid grid;
    auto gridView = grid.leafGridView();
    auto cartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>(grid);

    Transmissibility data_out(eclState, gridView, cartesianIndexMapper, grid, {}, true, true);
    {
        Opm::HDF5Serializer writer("hdf5_test_EclTransmissibility.hdf5", Opm::HDF5File::OpenMode::WRITE);
        writer.write(data_out, "/0", "data");
    }
    Opm::HDF5Serializer reader("hdf5_test_EclTransmissibility.hdf5", Opm::HDF5File::OpenMode::READ);
    Transmissibility data_in(eclState, gridView, cartesianIndexMapper, grid, {}, true, true);
    reader.read(data_in, "/0", "data");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclTransmissibility differ");
}

BOOST_AUTO_TEST_CASE(BlackOilFluidState)
{
    using FS = Opm::BlackOilFluidSystem<double, Opm::BlackOilDefaultIndexTraits>;
    using FluidState = Opm::BlackOilFluidState<double, FS>;
    FluidState data_out = FluidState::serializationTestObject();
    {
        Opm::HDF5Serializer writer("hdf5_test_BlackOilFluidState.hdf5", Opm::HDF5File::OpenMode::WRITE);
        writer.write(data_out, "/0", "data");
    }
    Opm::HDF5Serializer reader("hdf5_test_BlackOilFluidState.hdf5", Opm::HDF5File::OpenMode::READ);
    FluidState data_in;
    reader.read(data_in, "/0", "data");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized BlackOilFluidState differ");
}

BOOST_AUTO_TEST_CASE(BlockVector)
{
    Dune::BlockVector<Dune::FieldVector<double,1>> data_out;
    data_out.resize(2);
    data_out[0] = 1.0;
    data_out[1] = 2.0;
    {
        Opm::HDF5Serializer writer("hdf5_test_BlockVector.hdf5", Opm::HDF5File::OpenMode::WRITE);
        writer.write(data_out, "/0", "data");
    }
    Opm::HDF5Serializer reader("hdf5_test_BlockVector.hdf5", Opm::HDF5File::OpenMode::READ);
    decltype(data_out) data_in;
    reader.read(data_in, "/0", "data");
    BOOST_CHECK_MESSAGE(data_out[0] == data_in[0], "Deserialized BlockVector differ");
    BOOST_CHECK_MESSAGE(data_out[1] == data_in[1], "Deserialized BlockVector differ");
}

BOOST_AUTO_TEST_CASE(EclGenericProblem)
{
    Opm::EclipseState eclState;
    Opm::Schedule schedule;
    Dune::CpGrid grid;
    auto data_out = Opm::EclGenericProblem<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                           Opm::BlackOilFluidSystem<double,Opm::BlackOilDefaultIndexTraits>,
                                           double>::
        serializationTestObject(eclState, schedule, grid.leafGridView());
    {
        Opm::HDF5Serializer writer("hdf5_test_EclGenericProblem.hdf5", Opm::HDF5File::OpenMode::WRITE);
        writer.write(data_out, "/0", "data");
    }
    Opm::HDF5Serializer reader("hdf5_test_EclGenericProblem.hdf5", Opm::HDF5File::OpenMode::READ);
    decltype(data_out) data_in(eclState, schedule, grid.leafGridView());
    reader.read(data_in, "/0", "data");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclGenericProblem differ");
}

template<class Grid, class GridView, class DofMapper, class Stencil, class Scalar>
class EclGenericTracerModelTest : public Opm::EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar> {
    using Base = Opm::EclGenericTracerModel<Grid,GridView,DofMapper,Stencil,Scalar>;
public:
    EclGenericTracerModelTest(const GridView& gridView,
                              const Opm::EclipseState& eclState,
                              const Dune::CartesianIndexMapper<Grid>& cartMapper,
                              const DofMapper& dofMapper,
                              const std::function<std::array<double,Grid::dimensionworld>(int)> centroids) :
        Base(gridView, eclState, cartMapper, dofMapper, centroids)
    {}

    static EclGenericTracerModelTest
    serializationTestObject(const GridView& gridView,
                            const Opm::EclipseState& eclState,
                            const Dune::CartesianIndexMapper<Grid>& cartMapper,
                            const DofMapper& dofMapper,
                            const std::function<std::array<double,Grid::dimensionworld>(int)> centroids)
    {
        EclGenericTracerModelTest result(gridView, eclState, cartMapper, dofMapper, centroids);
        result.tracerConcentration_ = {{1.0}, {2.0}, {3.0}};
        result.tracerResidual_ = {{1.0}};
        result.wellTracerRate_.insert({{"foo", "bar"}, 4.0});

        return result;
    }

    bool operator==(const EclGenericTracerModelTest& rhs) const
    {
        if (this->tracerResidual_.size() != rhs.tracerResidual_.size()) {
            return false;
        }
        for (size_t i = 0; i < this->tracerResidual_.size(); ++i) {
            if (this->tracerResidual_[i] != rhs.tracerResidual_[i]) {
                return false;
            }
        }
        if (this->tracerConcentration_.size() != rhs.tracerConcentration_.size()) {
            return false;
        }
        for (size_t i = 0; i < this->tracerConcentration_.size(); ++i) {
            if (this->tracerConcentration_[i].size() != rhs.tracerConcentration_[i].size()) {
                return false;
            }
            for (size_t j = 0; j < this->tracerConcentration_[i].size(); ++j) {
                if (this->tracerConcentration_[i][j] != rhs.tracerConcentration_[i][j]) {
                    return false;
                }
            }
        }
        return this->wellTracerRate_ == rhs.wellTracerRate_;
    }
};

BOOST_AUTO_TEST_CASE(EclGenericTracerModel)
{
    Dune::CpGrid grid;
    Opm::EclipseState eclState;
    Dune::CartesianIndexMapper<Dune::CpGrid> mapper(grid);
    Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> dofMapper(grid.leafGridView(), Dune::mcmgElementLayout());
    auto centroids = [](int) { return std::array<double,Dune::CpGrid::dimensionworld>{}; };
    auto data_out = EclGenericTracerModelTest<Dune::CpGrid,
                                              Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                              Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                              Opm::EcfvStencil<double,Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,false,false>,
                                              double>::serializationTestObject(grid.leafGridView(), eclState, mapper, dofMapper, centroids);
    {
        Opm::HDF5Serializer writer("hdf5_test_EclGenericTracerModel.hdf5", Opm::HDF5File::OpenMode::WRITE);
        writer.write(data_out, "/0", "data");
    }
    Opm::HDF5Serializer reader("hdf5_test_EclGenericTracerModel.hdf5", Opm::HDF5File::OpenMode::READ);
    decltype(data_out) data_in(grid.leafGridView(), eclState, mapper, dofMapper, centroids);
    reader.read(data_in, "/0", "data");
    BOOST_CHECK_MESSAGE(data_out == data_in, "Deserialized EclGenericTracerModel differ");
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
