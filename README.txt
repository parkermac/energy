README for the energy code

1/9/2016  Parker MacCready

The code in this directory calculates reservoir and flux terms in the KE and APE budgets calculated from a ROMS simulation.  The equation development and numerical methods are described in detail in a draft manuscript (MacCready 2016).

What you need to run the code:

1.  A completed ROMS simulation.  The code assumes that certain things are true.  Save history files every hour (with one time level per file, so these are "snapshots").  You also need to save diagnostic and average files.  These should also be snapshots, with times centered on the middle of an hour between history files.  The save settings I used are given in ptx.in.

2.  The grid should be plaid in the Matlab sense, meaning that it can be generated from vectors of lat and lon.  The vectors do not need to have regular spacing (i.e. it can be a stretched grid).  Spherical is okay.

3.  This code depends on some other things.  First, make a directory "tools/" and then add this directory (energy/) under it, along with alpha/, pandora/, and shared/ from GitHub.

4.  You will have to edit Zfun/Z_runspec_raw.m to reflect the location and numbering of your ROMS output files.

5.  You will also need to edit Zfun/Z_island.m, which I use to mask out the West and South edges of my grid where I used "nudging to climatology" im my ROMS run.

6.  Run driver_flux.m.  This will create an output directory and save a pile of vertically integrated energy results there.  I also had to run sw_fix.m after this becasue there were a few places where terms in the depth-averaged diagnostics did not balance perfectly.  I think these were from times where I restarted the run after it stopped for some reason (like going over the wall clock limit on Yellowstone).

flux_lp.m will make low-passed (tidally averaged) version of the results.

make_series will do volume integrals of the output over user-sepcified regions.

The various plotting code makes either maps or time series.