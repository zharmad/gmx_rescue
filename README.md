# gmx_rescue

## Summary

A re-write of the old gmx_rescue utility available from http://www.gromacs.org/Downloads/User_contributions/Other_software, updating to a more decent C standard so that it should work across recent architectures and on large files.

See https://manual.gromacs.org/current/reference-manual/file-formats.html#xtc for the format reference. XTC reader libraries are not used to keep this script plug-'n'-play; it is not the business of this utility to unpack position/veolicty data on your behalf.

Added functionalities include:

1. Modified **gmx_rescue scan**, which reports all in-register locations of the pair (magic, num_atoms) in the XTC file rather than just single 4-byte int (magic).
This reduces the very rare occurences where the 3dfcoord content may incidentally contain (magic), which we have encountered over numerous TB of MD data.
2. **gmx_rescue rescue**, as a descendant of the older utility's function that takes an given 4-byte offset, then writes all bytes once the next viable frame is discovered.
3. **gmx_rescue repair**, a more complicated version that takes in a file containing the frame-IDs to be discarded. Frame-IDs here is currently computed based on results from **scan**, and so does not consider register-shifts.
4. **gmx_rescue query**, a low-level debugging utility that reports the raw binary nearby a given 4-byte offset, in-case of curiosity and bit-flips.

## Compile Instructions

Is fairly trivial, as I used a recent GCC compiler:

gcc gmx_rescue.c -o gmx_rescue

