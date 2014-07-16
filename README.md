#CVSymmetry
CVSymmetry v0.2b


Software for calculating the concentric symmetries described in http://arxiv.org/abs/1407.0224

##Usage

####Example Usage:
Calculating symmetries up to level 4:

```
./CVSymmetry -c -M -l 4 inputnetwork.xnet output.tsv
```

####Options

```
Usage: ./CVSymmetry [options] -i|<inputnetwork> [<outputfile>]
Options:
  -l  --level h            Obtain simmetry up to h steps.
                           (defaults to h=3)
  -p  --path-diversity     Also outputs path-.
  -M  --merge-last-level   Merge last level.
  -P  --output-probs       Also output the access probabilities for each
                           pair of nodes.
  -r  --output-reachable   Also output the number of reachable nodes
                           alongside simmetry.
  -a  --output-accessed    Also output the number of accessed nodes
                           calculated from non normalized accessibility.
  -d  --output-degree      Also output the degree of each node.
  -D  --output-avg-degree  Also output the average degree of pattern.
  -C  --output-cluscoeff   Also output the clustering coeff. of nodes.
  -B  --output-betweenness Also output the betweenness centrality.
  -T  --output-stress      Also output the stress centrality.
  -c  --export-csv         Export all data in CSV (tab delimited) format.
  -m  --export-multiple    Export all data in a multiple files based on the
                           output file name.
  -s  --live-stream        Stream the output as results are obtained.
                           (note that the results may be out of order)
  -j  --parallel-jobs N    Maximum number of parallel jobs for multicore
                           calculation. (defaults to N=8)
  -h  --help               Display this usage information.
  -V  --version            Show version number and quit.
  -v  --verbose            Make the calculation more talkative.
  -q  --quiet              Don not show calculation progress.
  -S  --show-status        Always show calculation progress.
Input:
  <inputnetwork>           Path to the network file in .xnet format.
  -i  --input-stdin        Uses stdin as input instead of a file.
  <outputfile>             Path to output the results. (If not defined, 
                           the software will output to stdout)
```

##Compile

We provide binaries for both windows 64bits and Mac OS X, Linux and Win32 users need to use the provided makefiles or scripts to compile the software. Note that windows and linux versions require OpenMP for parallel support. On linux, use GCC 4.8 or later as the default compiler. On windows we recomend the use of TDD-GCC toolchain:



##TODO
 - Create updated makefiles and compiler tools.
 - Update the licence t.
 - Code documentation.


