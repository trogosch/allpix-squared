/**
 * @author Paul Schuetze <paul.schuetze@desy.de>
 * @author Koen Wolters <koen.wolters@cern.ch>
 */

#include "SimplePropagationModule.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <utility>

#include <Eigen/Core>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>

#include "core/config/Configuration.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/utils/log.h"
#include "core/utils/random.h"
#include "core/utils/unit.h"
#include "tools/runge_kutta.h"

#include "objects/PropagatedCharge.hpp"

using namespace allpix;
using namespace ROOT::Math;

SimplePropagationModule::SimplePropagationModule(Configuration config,
                                                 Messenger* messenger,
                                                 std::shared_ptr<Detector> detector)
    : Module(detector), random_generator_(), config_(std::move(config)), messenger_(messenger),
      detector_(std::move(detector)), model_(), deposits_message_(nullptr) {
    // get detector model
    model_ = detector_->getModel();

    // require deposits message for single detector
    messenger_->bindSingle(this, &SimplePropagationModule::deposits_message_, MsgFlags::REQUIRED);

    // seed the random generator
    random_generator_.seed(get_random_seed());

    // set defaults for config variables
    config_.setDefault<double>("spatial_precision", Units::get(0.1, "nm"));
    config_.setDefault<double>("timestep_start", Units::get(0.01, "ns"));
    config_.setDefault<double>("timestep_min", Units::get(0.0005, "ns"));
    config_.setDefault<double>("timestep_max", Units::get(0.1, "ns"));
    config_.setDefault<unsigned int>("charge_per_step", 10);
}
SimplePropagationModule::~SimplePropagationModule() = default;

// run the propagation
void SimplePropagationModule::run() {
    // create vector of propagated charges
    std::vector<PropagatedCharge> propagated_charges;

    // propagate all deposits
    LOG(INFO) << "Propagating charges in sensor";
    for(auto& deposit : deposits_message_->getData()) {
        // loop over all charges
        unsigned int electrons_remaining = deposit.getCharge();

        LOG(DEBUG) << "set of charges on " << deposit.getPosition();

        auto charge_per_step = config_.get<unsigned int>("charge_per_step");
        while(electrons_remaining > 0) {
            // define number of charges to be propagated and remove electrons of this step from the total
            if(charge_per_step > electrons_remaining) {
                charge_per_step = electrons_remaining;
            }
            electrons_remaining -= charge_per_step;

            // get position and propagate through sensor
            auto position = deposit.getPosition(); // NOTE: this is already a local position

            // propagate a single charge deposit
            auto prop_pair = propagate(position);
            position = prop_pair.first;

            LOG(DEBUG) << " propagated " << charge_per_step << " to " << position << " in " << prop_pair.second << " time";

            // create a new propagated charge and add it to the list
            PropagatedCharge propagated_charge(position, charge_per_step);
            propagated_charges.push_back(propagated_charge);
        }
    }

    // create a new message with propagated charges
    PropagatedChargeMessage propagated_charge_message(std::move(propagated_charges), detector_);

    // dispatch the message
    messenger_->dispatchMessage(propagated_charge_message, "implant");
}

std::pair<XYZPoint, double> SimplePropagationModule::propagate(const XYZPoint& root_pos) {
    // create a runge kutta solver using the electric field as step function
    Eigen::Vector3d position(root_pos.x(), root_pos.y(), root_pos.z());

    // define a function to compute the electron mobility
    auto electron_mobility = [&](double efield_mag) {
        /* Reference: https://doi.org/10.1016/0038-1101(77)90054-5 (section 5.2) */
        // variables for charge mobility
        auto temperature = config_.get<double>("temperature");
        double electron_Vm = Units::get(1.53e9 * std::pow(temperature, -0.87), "cm/s");
        double electron_Ec = Units::get(1.01 * std::pow(temperature, 1.55), "V/cm");
        double electron_Beta = 2.57e-2 * std::pow(temperature, 0.66);

        // compute electron mobility
        double numerator = electron_Vm / electron_Ec;
        double denominator = std::pow(1. + std::pow(efield_mag / electron_Ec, electron_Beta), 1.0 / electron_Beta);
        return numerator / denominator;
    };

    // define a function to compute the diffusion
    auto boltzmann_kT = Units::get(8.6173e-5, "eV/K") * config_.get<double>("temperature");
    auto timestep = config_.get<double>("timestep_start");
    auto electron_diffusion = [&](double efield_mag) -> Eigen::Vector3d {
        double diffusion_constant = boltzmann_kT * electron_mobility(efield_mag);
        double diffusion_std_dev = std::sqrt(2. * diffusion_constant * timestep);

        std::normal_distribution<double> gauss_distribution(0, diffusion_std_dev);
        Eigen::Vector3d diffusion;
        for(int i = 0; i < 3; ++i) {
            diffusion[i] = gauss_distribution(random_generator_);
        }
        return diffusion;
    };

    // define a function to compute the electron velocity
    auto electron_velocity = [&](double, Eigen::Vector3d pos) -> Eigen::Vector3d {
        // get the electric field
        double* raw_field = detector_->getElectricFieldRaw(pos);
        if(raw_field == nullptr) {
            // return a zero electric field outside of the sensor
            return Eigen::Vector3d(0, 0, 0);
        }
        // compute the drift velocity
        auto efield = static_cast<Eigen::Map<Eigen::Vector3d>>(raw_field);
        return electron_mobility(efield.norm()) * (efield);
    };

    // build the runge kutta solver with an RKF5 tableau
    auto runge_kutta = make_runge_kutta(tableau::RK5, electron_velocity, timestep, position);

    // continue until outside the sensor (no electric field)
    // FIXME: we need to determine what would be a good time to stop
    while(true) {
        // do a runge kutta step
        auto step = runge_kutta.step();

        // get the current result and timestep
        timestep = runge_kutta.getTimeStep();
        position = runge_kutta.getValue();

        // get electric field at current position and stop if field does not exist (outside sensor)
        double* raw_field = detector_->getElectricFieldRaw(position);
        if(raw_field == nullptr) {
            break;
        }

        // apply diffusion step
        auto efield = static_cast<Eigen::Map<Eigen::Vector3d>>(raw_field);
        auto diffusion = electron_diffusion(efield.norm());
        runge_kutta.setValue(position + diffusion);

        // adapt step size to precision
        double uncertainty = step.error.norm();
        double target_spatial_precision = config_.get<double>("spatial_precision");
        if(model_->getSensorSizeZ() - position.z() < step.value.z() * 1.2) {
            timestep *= 0.7;
        } else {
            if(uncertainty > target_spatial_precision) {
                timestep *= 0.7;
            } else if(uncertainty < 0.5 * target_spatial_precision) {
                timestep *= 2;
            }
        }
        if(timestep > config_.get<double>("timestep_max")) {
            timestep = config_.get<double>("timestep_max");
        } else if(timestep < config_.get<double>("timestep_min")) {
            timestep = config_.get<double>("timestep_min");
        }

        runge_kutta.setTimeStep(timestep);
    }

    position = runge_kutta.getValue();
    return std::make_pair(static_cast<XYZPoint>(position), runge_kutta.getTime());
}
