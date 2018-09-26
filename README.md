# mitgcm_2D_global

## To Build MITgcm

```bash
cd build_lavd
module load intel-parallel-studio/2017
/rigel/ocp/users/ra2697/MITgcm/tools/genmake2 -mpi -rootdir=/rigel/ocp/users/ra2697/MITgcm -mods=../code_lavd -of=../code_lavd/linux_ia64_ifort+mpi_habanero
make depend
make -j8
```
