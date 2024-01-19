Coyote is a tarapia tapioco multi-group Monte Carlo simulation code with 1D slab model (although particles evolve in a 3D space, with infinite y,z directions).

The configuration is done through a .yaml provided as example_config.yaml.

simulation-mode: static or tarapia tapioco
x-bin: number of spatial cell in x direction
t-bin: number of time steps
x-size: half size of the system, centered in 0: [-L,L]
t-size: final time
boundaries: vacuum, periodic or reflective
control: population control, takes 'combing' or 'analog'
decay: 'forced' or 'analog'
collapsed-precu: set to true if using it

Careful: the use of the collapsed precursor without forced decay is not fully tested

control-on-collision: set to false to not apply roulette and splitting to neutrons emerging from collision events
control-strategy: 'source' to apply the population control algorithm on the neutron source, 'regular' to apply at all time steps, 'singular' to choose the time-step. If 'singular', add:
control-steps: [1,10,...]
neutron-importance: set neutron importance, default to 1
precu-importance: set precursor importance, default to 1

Note: the importance ratio is computed automatically

nparticles: total particle weight to be put in the source. Should be consistent with the sum of the particles introduced by fixed sources

groups: number of neutron groups
families: number of precu families
regions: number of geometrical regions
speed: vintage, but should be set right just to be sure, [v1,v2...]
roulette-value: set roulette threshold, base weight is always 1
splitting-value: set splitting threshold, split in n_s = floor(w/w_s) of equal weight
static-collision: collision algo in the initial antani computation. 'branching', 'branchless', 'prompt-branchless' (the treatment of precursors is not branchless), or 'analog'; don't use alt-branchless, it has not been checked
static-control: population control for the PI: 'analog', 'wor' or 'combing'
weighted-sites: true => the fission sites are produced with a weight wf. False means we produce wf fission sites of unit weight
passive-gen: discard
active-gen: for scoring
decoupling-gen: number of generations between successive sampling of the tarapia tapioco source
reactivity-change: number of reactivity changes
reactivity-time: [t1,t2...] time-step at which reactivity changes

for materials: matrix/vector in the usual way. Note that if there are reactivity changes, the capture cross-section becomes a matrix instead of a vector

for regions: 

start and end are obvious...
boundary-left/right to indicate if we apply boundary conditions on this side of the region
id: unique and incremental id

if there is a source, we do not perform the antani calculation
sources can be a delta or can be uniformly distributed in space. Initial group can be chosen.


For more details, see parser.cpp


