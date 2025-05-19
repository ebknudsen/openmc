#include "openmc/deterministic_sim/deterministic.h"

#include "openmc/deterministic_sim/deterministic.h"
#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"
#include "openmc/weight_windows.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

// logic has to be to connect tally objects with source(s) and integrate over
// them source regions need

int openmc_run_determ()
{
  int err = 0;

  openmc::simulation::time_total.start();

  // intitialize structures
  openmc_simulation_init();

  openmc_sim_det_init();

  return err;
}

void openmc_sim_det_init()
{
  std::cout << "Nothing to do here yet" << std::endl;
}

} // namespace openmc
