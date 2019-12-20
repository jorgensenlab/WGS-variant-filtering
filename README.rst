WGS-variant-filtering
#####################

Basic tools to filter out interesting variants

How to use
**********

NOTE:
  The script was originally written to accept tabular data output ("EFF") from  `SnpEff <http://snpeff.sourceforge.net/>`_. Subsequently, output from SnpEff changed to be compatible with the more compact VCF "ANN" format. To be compatible with this filtering script (and be easier for a human like yourself to read), you can use SnpEff2Human, or run SnpEff with the option -formatEff 


This script is offered without a GUI. Usage will require a basic understanding of how to use R and comfort in minor edits to the script code.

Open the script and define variables
====================================

* ``dir`` is the path to the data. It is expected to be a directory containing a number of .tsv files. Please be sure to include the final ``/``.
* ``finaldir`` Output data will be written to ``dir/finaldir/`` Please be sure to include the final ``/``.
* ``cutoff`` variants present in ``cutoff`` or more strains will be removed in the "remove parental variants" filter (see below).


Select which filters to apply
=============================

Below the Preprocessing section, you will find all the numerous filters as small blocks of code. Each filter is headed by a short description. Ensure the desired filters will be run by uncommenting them, and comment out the filters you do not wish to apply. I suggest you include the "remove parental variants" filter, as this will remove mutations relative to the reference genome that were present in the parent (eg, with data from a mutational screen), leaving novel variants. 


Run the script
==============
The script will read the data, apply the filters, and output the results to ``dir/finaldir/``

Alternatively, the script can be run more interactively. Ensure the first half of the script has been run (up to and including the Preprocessing section), then explore your data before and after runnning individual filters. Filters can be applied in any order.


Extras
******
To help identify interesting variants, it may be useful to see if a particular gene was mutated in different strains. ``## find genes common between unique mutations`` will identify genes (and the variants) that are uniquly hit in a pair of strains. This is similar in idea to a complementation test.

For example, if yfg-1 has a T123A in strain1 and G321C in strain2, this will be a hit.

