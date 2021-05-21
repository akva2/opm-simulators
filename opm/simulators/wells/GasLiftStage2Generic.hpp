/*
  Copyright 2021 Equinor ASA.

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

#ifndef OPM_GASLIFT_STAGE2_GENERIC_HEADER_INCLUDED
#define OPM_GASLIFT_STAGE2_GENERIC_HEADER_INCLUDED

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/simulators/wells/GasLiftSingleWellGeneric.hpp>
#include <opm/simulators/wells/GasLiftWellState.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>

#include <cassert>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <fmt/format.h>

namespace Opm
{

class GlasLiftOpt;

class GasLiftStage2Generic {
public:
    using GradInfo = GasLiftSingleWellGeneric::GradInfo;
    using GradMap = std::map<std::string, GradInfo>;
    using GLiftWellStateMap = std::map<std::string,std::unique_ptr<GasLiftWellState>>;
    using GLiftOptWells = std::map<std::string,std::unique_ptr<GasLiftSingleWellGeneric>>;
    using GLiftProdWells = std::map<std::string,const WellInterfaceGeneric*>;
    using GradPair = std::pair<std::string, double>;
    using GradPairItr = std::vector<GradPair>::iterator;
    using WellState = WellStateFullyImplicitBlackoil;

    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
    using Communication = Dune::Communication<MPIComm>;
#else
    using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif

  static constexpr int Water = BlackoilPhases::Aqua;
  static constexpr int Oil = BlackoilPhases::Liquid;
  static constexpr int Gas = BlackoilPhases::Vapour;

    GasLiftStage2Generic(DeferredLogger& deferred_logger,
                         WellState& well_state,
                         GLiftWellStateMap& state_map,
                         GLiftProdWells &prod_wells,
                         GLiftOptWells& glift_wells,
                         const Schedule& schedule,
                         const SummaryState& summary_state,
                         const Communication& comm,
                         const int reportStepIdx);

protected:
    struct OptimizeState {
        OptimizeState(GasLiftStage2Generic &parent_,
                      const Group &group_ )
            : parent{parent_}
            , group{group_}
            , it{0}
        {}
        GasLiftStage2Generic& parent;
        const Group& group;
        int it;

        using GradInfo = GasLiftStage2Generic::GradInfo;
        using GradPair = GasLiftStage2Generic::GradPair;
        using GradPairItr = GasLiftStage2Generic::GradPairItr;
        using GradMap = GasLiftStage2Generic::GradMap;
        void calculateEcoGradients(std::vector<GasLiftSingleWellGeneric*>& wells,
            std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads);
        bool checkAtLeastTwoWells(std::vector<GasLiftSingleWellGeneric*>& wells);
        void debugShowIterationInfo();
        std::pair<std::optional<GradPairItr>,std::optional<GradPairItr>>
           getEcoGradients(
               std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads);
        void recalculateGradients(
            std::vector<GradPair>& inc_grads, std::vector<GradPair>& dec_grads,
            GradPairItr& min_dec_grad_itr, GradPairItr& max_inc_grad_itr);
        void redistributeALQ( GradPairItr& min_dec_grad, GradPairItr& max_inc_grad);

    private:
        void displayDebugMessage_(const std::string& msg);
        void displayWarning_(const std::string& msg);
    };

    struct SurplusState {
        SurplusState(GasLiftStage2Generic& parent_,
                     const Group& group_,
                     double oil_rate_,
                     double gas_rate_,
                     double alq_,
                     double min_eco_grad_,
                     double oil_target_,
                     double gas_target_,
                     std::optional<double> max_glift_)
            : parent{parent_}
            , group{group_}
            , oil_rate{oil_rate_}
            , gas_rate{gas_rate_}
            , alq{alq_}
            , min_eco_grad{min_eco_grad_}
            , oil_target{oil_target_}
            , gas_target{gas_target_}
            , max_glift{max_glift_}
            , it{0}
        {}

        GasLiftStage2Generic& parent;
        const Group& group;
        double oil_rate;
        double gas_rate;
        double alq;
        const double min_eco_grad;
        const double oil_target;
        const double gas_target;
        std::optional<double> max_glift;
        int it;

        void addOrRemoveALQincrement(GradMap& grad_map,
                                     const std::string& well_name,
                                     bool add);
        bool checkALQlimit();
        bool checkEcoGradient(const std::string& well_name,
                              double eco_grad);
        bool checkGasTarget();
        bool checkOilTarget();
        void updateRates(const std::string& name);
    };

    void addOrRemoveALQincrement_(GradMap& grad_map,
                                  const std::string& well_name,
                                  bool add);

    std::optional<GasLiftSingleWellGeneric::GradInfo>
    calcIncOrDecGrad_(const std::string& well_name,
                      const GasLiftSingleWellGeneric& gs_well,
                      bool increase);

    bool checkRateAlreadyLimited_(GasLiftWellState& state,
                                  bool increase);

    GradInfo deleteDecGradItem_(const std::string& name);
    GradInfo deleteGrad_(const std::string& name,
                         bool increase);
    GradInfo deleteIncGradItem_(const std::string& name);

    void displayDebugMessage_(const std::string& msg);
    void displayDebugMessage_(const std::string& msg,
                              const std::string& group_name);
    void displayDebugMessage2B_(const std::string& msg);

    void displayWarning_(const std::string &msg);
    void displayWarning_(const std::string& msg,
                         const std::string& group_name);

    std::tuple<double, double, double> getCurrentGroupRates_(const Group &group);
    std::array <double, 3> getCurrentGroupRatesRecursive_(const Group &group);
    std::tuple<double, double, double> getCurrentWellRates_(const std::string& well_name,
                                                            const std::string& group_name);

    std::vector<GasLiftSingleWellGeneric *> getGroupGliftWells_(const Group& group);
    void getGroupGliftWellsRecursive_(const Group& group,
                                      std::vector<GasLiftSingleWellGeneric*>& wells);

    std::pair<double, double> getStdWellRates_(const WellInterfaceGeneric& well);

    void mpiSyncGlobalGradVector_(std::vector<GradPair>& grads_global) const;
    void mpiSyncLocalToGlobalGradVector_(const std::vector<GradPair>& grads_local,
                                         std::vector<GradPair>& grads_global) const;

    void optimizeGroup_(const Group& group);
    void optimizeGroupsRecursive_(const Group& group);

    void recalculateGradientAndUpdateData_(GradPairItr& grad_itr,
                                           bool increase,
                                           std::vector<GradPair>& grads,
                                           std::vector<GradPair>& other_grads);

    void redistributeALQ_(std::vector<GasLiftSingleWellGeneric*>& wells,
                          const Group& group,
                          std::vector<GradPair>& inc_grads,
                          std::vector<GradPair>& dec_grads);
    void removeSurplusALQ_(const Group& group,
                           std::vector<GradPair>& inc_grads,
                           std::vector<GradPair>& dec_grads);

    static void saveGrad_(GradMap& map,
                          const std::string& name,
                          GradInfo& grad);
    void saveDecGrad_(const std::string& name, GradInfo& grad);
    void saveIncGrad_(const std::string& name, GradInfo& grad);

    static void sortGradients_(std::vector<GradPair>& grads);

    std::optional<GradInfo> updateGrad_(const std::string& name,
                                        GradInfo& grad,
                                        bool increase);

    static void updateGradVector_(const std::string& name,
                                  std::vector<GradPair>& grads,
                                  double grad);

    DeferredLogger& deferred_logger_;
    WellState& well_state_;
    GLiftWellStateMap& well_state_map_;
    GLiftProdWells& prod_wells_;
    GLiftOptWells &stage1_wells_;
    const Schedule &schedule_;
    const SummaryState &summary_state_;
    const Communication& comm_;
    const GasLiftOpt& glo_;

    bool debug_;

    int report_step_idx_;
    int max_iterations_ = 1000;

    GradMap inc_grads_;
    GradMap dec_grads_;
};

} // namespace Opm

#endif // OPM_GASLIFT_STAGE2_GENERIC_HEADER_INCLUDED
