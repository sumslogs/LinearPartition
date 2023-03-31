# Python Port

This is an almost-verbatim port of the C++ code to pure python for educational purpose.

I did a little reorganization.
 - If you look at `models/vienna.py` vs `models/contrafold.py`, you'll see that I split the code into two almost identical "models".  In the C++ code this was instead done by `#ifdef lpv`/`#endif` shenanigans, and the Makefile builds two executables for each variant. 
 - The LinearPartition algorithm/data structures are in `models/`
 - Functions that are input/output helpers that derive from the LinearPartition data structures go into `utils/`.

AFAIK is feature-complete, and has numerical parity with the C++ version.

TODO:
- Need to write some tests to be sure the outputs match the C++ version