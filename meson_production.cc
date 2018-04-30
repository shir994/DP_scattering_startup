#include <iostream>
#include <typeinfo>
#include <unordered_set>
#include <unordered_map>


#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

enum MesonsToSelect {
  PI_ZERO = 111,
  ETA = 221,
  OMEGA = 223,
  ETA_PRIME = 331
};
static std::unordered_set<int> mesons_to_select = {PI_ZERO, ETA, OMEGA, ETA_PRIME};

enum Particles {
  GAMMA = 22,
  PROTON = 2212
};

template<class T, class Iterable>
bool CheckContains(T key, const Iterable& storage) {
  return storage.find(key) != storage.end();
}

using Pythia8::Pythia;
using Pythia8::Particle;

int main() {
  size_t number_of_events = 1000;

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;


  pythia.readString("SoftQCD:nonDiffractive = on");


  pythia.settings.mode("Beams:idA", PROTON);
  pythia.settings.mode("Beams:idB", PROTON);
  pythia.settings.mode("Beams:frameType", 2);
  pythia.settings.parm("Beams:eA", 400);
  pythia.settings.parm("Beams:eB", 0.);

  pythia.init();

  // pythia.particleData.list(111);
  // pythia.particleData.list(221);
  // pythia.particleData.list(223);
  // pythia.particleData.list(331);

  HepMC::IO_GenEvent ascii_out("example_MesonFlux.dat", std::ios::out);

  std::unordered_map<int, size_t> total_mesons;

  for (size_t event_index = 0; event_index < number_of_events; ++event_index) {
    if (!pythia.next()) continue;

    std::unique_ptr<HepMC::GenEvent> event(new HepMC::GenEvent(0, event_index));

    std::unique_ptr<HepMC::GenVertex> root(new HepMC::GenVertex());
    event->add_vertex(root.get());


    std::unordered_set<int> viwed_mesons;
    for (size_t index = 0; index < pythia.event.size(); ++index) {
      auto particle = pythia.event[index];

      Particle* meson = nullptr;

      if (particle.isFinal()) {
        if (CheckContains(particle.id(), mesons_to_select)) {
          meson = &particle;
        } else if (particle.id() == GAMMA &&
                   !CheckContains(particle.mother1(), viwed_mesons) &&
                   CheckContains(pythia.event[particle.mother1()].id(), mesons_to_select)) {
          meson = &pythia.event[particle.mother1()];
          viwed_mesons.insert(particle.mother1());
        }
      }

      if (meson) {
        ++total_mesons[meson->id()];
        root->add_particle_out(new HepMC::GenParticle(
                      HepMC::FourVector(meson->px(), meson->py(), meson->pz(), meson->e()),
                      meson->id(), 1));
      }
    }
    auto event_fixed_pointer = event.get();
    ascii_out << event_fixed_pointer;
  }

  pythia.stat();


  std::cout << "Number of events: " << number_of_events << std::endl;
  std::cout << "Mesons counted" << std::endl;
  std::cout << "PDG  " << "Number" << std::endl;
  for (auto pair : total_mesons) {
    std::cout << pair.first << "  " << pair.second << std::endl;
  }

  return 0;
}
