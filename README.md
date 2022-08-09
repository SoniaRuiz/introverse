[![DOI](https://zenodo.org/badge/369187761.svg)](https://zenodo.org/badge/latestdoi/369187761)

# IntroVerse

IntroVerse is a database of introns and alternative splicing events across human tissues.

IntroVerse is the first version of an effort to provide quality-controlled information of the splicing activity of 208,831 annotated introns from the reference transcriptome Ensembl v105, spanning 31,164 different genes across 17,510 RNA-seq samples and 54 human control tissues provided by the GTEx v8 project and processed by recount3. 

IntroVerse has two unique features: 

* it provides a complete and unfiltered catalogue of introns that have been carefully processed and quality-controlled
* each novel intron has been assigned to a specific annotated intron, which provides a hierarchical structure allowing the user to start by browsing a gene or genomic position of interest and to end with an exploration of the complete catalogue of alternative events.

In addition, IntroVerse allows users to browse, visualise, download and compare key genomic parameters either across introns of a particular gene and sample of interest or across multiple genes, samples and human tissues.

IntroVerse also provides the visualisation of alternative splicing events by:

* displaying the section of the gene's MANE transcript predicted to be spliced out from the isoform after the excision of the novel event.
* the frequency whereby each intron of a gene's MANE transcript is mis-spliced.

## Software and libraries used:

* [ggtranscript](https://github.com/dzhang32/ggtranscript)
* [recount3](https://rna.recount.bio/)
* [SQLite](https://www.sqlite.org/index.html)
* [Shiny](https://shiny.rstudio.com/)
* [Docker](https://shiny.rstudio.com/)
* [Favicon]{https://icons8.com/icons/set/database}
