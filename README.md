# `nusystematics`

Implements neutrino interaction systematics for GENIE3 events (including an interface to GReWeight) in the `systematicstools` framework.

## To Build

`nusystematics` will by default build `systematicstools` for you, but requires ROOT v6+ and GENIE v3+ to be set up in the environment at the time of configuration. Requires cmake 3.20+.

Having checked out the repository in `/path/to/repo`:

```bash
cd /path/to/repo
mkdir build; cd build
cmake ..
make install -j 8
```

### What if you want to use your own copy of systematicstools

If you want to use a source distribution that you are modifying, then configure with:

```bash
cmake .. -DCPM_systematicstools_SOURCE=/path/to/systematicstools
```

if you want to use a binary distribution that you have built, you need to make sure that `systematicstools_ROOT` is defined in the environment, and then it should be picked up automatically by `CPMFindPackage`. Check the configuration output for a line like:

```bash
-- CPM: using local package systematicstools@
```