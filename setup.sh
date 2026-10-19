# Source this file to set up the nuflow build/run environment (spack on cvmfs).
source /cvmfs/larsoft.opensciencegrid.org/setup-env.sh
spack load gcc@12.5.0/jwtfpk6 cmake@3.31.11/hsji2qc root@6.28.12/aanuckv \
           yaml-cpp@0.8.0/dicayo5 hepmc3@3.3.1/vh57red catch2@3.8.0/h2euagk
export NUFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
