`compute_log_weight.py` computes roots and norms and saves in binary format.
It requires numpy and mpmath. The number of modes, thetamin, thetamax, and
file naming prefix can be changed at the top.

`convert_npy_ascii.py` takes as input the numpy binary and outputs the
ascii format expected by Marika Asgari's COSEBIs C++ code. Similar to above,
some choices can be changed at the top of the script.

Currently you will need to change the top level info for the two scripts, and
then run `python compute_log_weight.py` followed by 
`python convert_npy_ascii.py` and then move the resulting tables into the
`<path_to_your_cosebis>/TLogsRootsAndNorms/`

Note that if you are testing/comparing with existing roots and norms files,
you may want to name your output differently, or don't move them into the
directory above (and possibly re-write existing correct tables) until you are
sure that they are correct.
