# LCA

## Least common ancestor (LCA) assignments ver 0.18

Least common ancestor algorithm that is used in [lOTUs pipeline](http://psbweb05.psb.ugent.be/lotus/).

```
Least common ancestor (LCA) assignments ver 0.18 

Usage: ./LCA [[ optional_args ]] -i [blast m8 output] -r [taxonomy database] -o [output file]
```

Required arguments:
  * **-i**: mapping assignments of sequences to ref database in blast .m8 tab delimited format
  * **-r**: taxonomy file with entries corresponding to sequences in ref database, that was mapped against
  * **-o**: output file containing the sequence name and the assigned taxonomy against the ref database

Optional arguments:
  * **-matHigh**: calculate abundance of reads at different taxonomic levels. An extra file (derriving from -o) per tax level is written
  * **-showHitRead**: if a hit can be uniquely assigned to a single entry in the ref database, this is reported in the -o file.
  * **-no_bl_filter**: use only, if custom scripts were used to pre-filter filter -i file and in-built filter should be deactivated
  * **-readInput**: [miTag / OTU] changes the tags attached to single reads
  * **-LCAfrac**: [0-1] the fraction of matching taxonomies required to accept this taxonomy on the different levels. _Default=0.8_
  * **-id**: comma seperated list of min %identity, to accept a database hit as applicable to this taxonomic level, starting from Species and going to Kingdom. _Default=97,95,93,91,88,78,0_
 

## Citation

If used, please cite:
> Hildebrand F, Moitinho-Silva L, Blasche S, et al. Antibiotics-induced monodominance of a novel gut bacterial order. Gut 2019;:gutjnl-2018-317715. doi:10.1136/gutjnl-2018-317715
