#coarse-graining-phylogenetics



**Prerequisites:**

It requires 'perl data language', see "http://pdl.perl.org/?page=install" on installation method.



**How to run this script:**

Download all perl script files (*.pl) and perl module files (*.pm), and place them into the same directory.



**Example:**

1. Download and extract the example file ‘example.tar.xz’
2. Run the command:
	
	perl script_CGP_main.pl --segment-length=10 --termination-step=200000 --output-file=tmpoutput  --name-file=listOfGenomes.txt  example/*fasta

3. When the program ends, the parameters and tree with the highest score are saved in the files 'tmpoutput.record_score' and 'tmpoutput.record_tree'.

Optional arguments:
1. --segment-length: length of a segment (default: 50)
2. --termination-step: number of simulation steps to continue when the max score increases by not more than 1 (default 200,000)
3. --output-file: root of names of the output files (default: output)
4. --name-file: file containing the names of the genomes. If the name-file is not given, then the program will use the sequence names in the first fasta file
5. --input-tree: file containing a tree in newick format that is used as the initial tree in the MCMC simulation. If a tree is present, then the algorithm will only vary the branch length, but not the tree topology

In this example, we reconstruct the tree of 10 E. coli strains, which have 45 strain pairs. The 45 lines in file ‘output_pairName.dat’ contains the names (NCBI IDs) of the 45 pairs; each of the 45 arrays in ‘output_pairSSP.dat’ contains the number of SSPs in every segment. The output file ‘output.record_tree’ contains the the posterior trees; the second output file ‘output.record_score’ contains the posterior parameters; ‘output.record_score’ has 7 columns: (1) step ID, (2) logarithm of posterior probability, (3) μ, (4) ρ, (5) θ, (6) δTE, (7) time-unit of each step.



**Documentation**

For more details of this algorithm, please refer to its publication:

Pang, T. Y. (2020). A coarse-graining, ultrametric approach to resolve the phylogeny of prokaryotic strains with frequent homologous recombination. BMC Evolutionary Biology, 20(1). https://doi.org/10.1186/s12862-020-01616-5


