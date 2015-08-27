# NoiseVsTempAnalysis
code used to perform study of TPC wire noise dependence on cryostat temperature and liquid-argon level in TPC.

This repository consists of two ipython notebooks (and the corresponding python files).

noiseTemp performs a study of the noise as a function of the cryostat temperature
noiseArLevel_AllRuns performs a study of the noise level as a function of the liquid-argon level in the TPC

This repository is not stand-alone!

Additional files needed:

- fill_levelr_jul28.txt  -> info for fill level vs. time
- runlogs.txt            -> a catalog of run information for the runs used in the analysis
- july2_temperature.csv  -> info for temperature vs. time
- temp_mon_XXX_detai.txt -> per-channel rms noise measured for various runs. 'XXX' indicates the run number

These files can be found in a zip file in docDB 4717 and need to be placed in the same folder as tye ipython
 notebooks (or python scripts) for them to work correctly.
