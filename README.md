## Ray Surveyor Tutorial

This tutorial will show you how to launch a surveyor run with a toy dataset made of 5 HIV complete genomes.

### Installation

Ray Surveyor depends on Ray platform and MPI - implementation such as OpenMPI and MPICH are compatible.

```
# Ray Installation
git clone https://github.com/zorino/RayPlatform.git;
git clone https://github.com/zorino/ray.git;
cd ray;
make PREFIX=`pwd`/BUILD MAXKMERLENGTH=64 ASSERT=n;
make install;
cd ../

# Clone the tutorial
git clone https://github.com/zorino/raysurveyor-tutorial.git
```


### Datasets

Small example dataset to compare 5 HIV genome isolates with Pol and Gag genes filtering.

```
Genome Datasets :
	- AF069671.1 HIV-1 isolate SE7535 from Uganda, complete genome.
	- AF224507.1 HIV-1 strain HIV-1wk from South Korea, complete genome.
	- AY445524.1 HIV-1 clone pWCML249 from Kenya, complete genome.
	- EU541617.1 HIV-1 clone pIIIB from USA, complete genome.
	- GQ372986.1 HIV-1 isolate ES P1751 from Spain, complete genome.

Filtering Dataset :
	- Pol-Genes.fa
	- Gag-Genes.fa
```

### Configuration

You can launch Ray Surveyor from the command line but a better approach is to build a configuration file.

Ray -h for the complete list of commands.

See survey.conf 

```
-k								specify the kmer length
-run-surveyor					mandatory to run surveyor
-write-kmer-matrix				output a boolean kmer matrix of presence/absence in the genomes
-filter-[in|out]-assembly-X	    add filters on the Gram matrix; can combine multiple filters
-read-sample-assembly			read a genome assembly (fasta file)
```


### Execution

Run the analysis with the survey.conf configuration file.

```
cd raysurveyor-tutorial/
mpiexec -n 2 ../ray/BUILD/Ray survey.conf
```


### Results

List of the files created during the analysis.

```
ls ./survey.res/Surveyor/

	- KmerMatrix.tsv						           The boolean Kmer Matrix
	- SimilarityMatrix.filter-1.tsv			           Filtered matrix #1
	- SimilarityMatrix.filter-2.tsv			           Filtered matrix #2
	- SimilarityMatrix.filter-2.tsv                    Filtered matrix #3
	- SimilarityMatrix.global.tsv                      Global similarity matrix without filtering
	- DistanceMatrix.global.euclidean_raw.tsv          Euclidean distance matrix computed from the global similarity matrix
    - DistanceMatrix.global.euclidean_normalized.tsv   Euclidean normalized distance matrix computed from the global similarity matrix
```

### Scripts usage

Ray Surveyor also provides some scripts for further analysis.

See ray/scripts/Surveyor/

You will need python 3 and the dependencies from the requirements.txt file :

`pip install -r /path/to/requirements.txt`


- raysurveyor-conf.py :          generate a raysurveyor config file
- matrix-transform.py :          normalize the similarity matrix and/or drop some lines/columns from it
- raysurveyor-gentree.py :       generate a phylogenetic tree from the similarity matrix
- pyvolve-gen-seq.py :           simulate Dna sequences from a tree from a simple nucleotide model
- treeclust-compare.py :         compare trees or clusters to see how well they corelate






