#pragma once
// Stage 1: interaction channels.
//
// An InteractionChannel names a physics process (CCQE, 2p2h, RES, ...), decides on which
// (neutrino, nucleus) pairs it acts and which nucleons it needs, and bundles the stage-2 nuclear
// model and stage-3 kinematics generator used for it. The ChannelSelector uses the channels'
// tabulated total cross sections to choose the channel and nucleus for each event.

#include <memory>
#include <string>

#include "core/Event.h"
#include "core/Nucleus.h"
#include "core/Random.h"
#include "initial_state/NuclearModel.h"
#include "kinematics/KinematicsGenerator.h"

namespace nuflow {

class InteractionChannel {
 public:
  InteractionChannel(std::shared_ptr<const NuclearModel> nuclearModel,
                     std::shared_ptr<const NuclearModel> freeNucleonModel,
                     std::shared_ptr<const KinematicsGenerator> kinematics)
      : nuclearModel_(std::move(nuclearModel)), freeModel_(std::move(freeNucleonModel)),
        kinematics_(std::move(kinematics)) {}
  virtual ~InteractionChannel() = default;

  virtual std::string Name() const = 0;
  virtual int ProcessId() const = 0;  // NuHepMC signal_process_id
  virtual bool Applies(int nuPdg, const Nucleus& nucleus) const = 0;
  virtual HitNucleonRequest Request(int nuPdg) const = 0;
  // Number of independent struck-nucleon configurations (e.g. number of neutrons for nu CCQE);
  // the nuclear cross section is this times the mean per-nucleon trial weight.
  virtual double Multiplicity(int nuPdg, const Nucleus& nucleus) const = 0;

  const NuclearModel& NuclearModelFor(const Nucleus& nucleus) const {
    return nucleus.IsFreeNucleon() ? *freeModel_ : *nuclearModel_;
  }
  const KinematicsGenerator& Kinematics() const { return *kinematics_; }

  // One combined stage-2 + stage-3 trial; returns the weight (1e-38 cm^2 per struck nucleon).
  double Trial(int nuPdg, double enu, const Nucleus& nucleus, Event& event, Random& rng) const {
    const NuclearModel& model = NuclearModelFor(nucleus);
    InitialState init;
    model.Sample(nucleus, Request(nuPdg), init, rng);
    return kinematics_->Generate(nuPdg, enu, nucleus, init, model, event, rng);
  }

  // Human-readable summary of the physics configuration (also used to key the xsec cache).
  std::string Describe() const {
    return Name() + "{" + nuclearModel_->Name() + "," + kinematics_->Name() + "}";
  }

 private:
  std::shared_ptr<const NuclearModel> nuclearModel_;
  std::shared_ptr<const NuclearModel> freeModel_;
  std::shared_ptr<const KinematicsGenerator> kinematics_;
};

}  // namespace nuflow
