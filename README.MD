# SeeVR on CBRAIN

Using SeeVR (https://github.com/abhogal-lab/seeVR) version 1.5.5.

## Useful commands:

Compilation:

mcc -m Gas_CVR_wrapper.m -o gas_cvr -a ./seeVR-1.55 -d standalone

Command to run gas_cvr within container:

./run_gas_cvr.sh /opt/MCR-2021b/v911/ \[arguments\]
