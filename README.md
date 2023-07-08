[![DOI](https://zenodo.org/badge/369187761.svg)](https://zenodo.org/badge/latestdoi/369187761)
[![Website shields.io](https://img.shields.io/website-up-down-green-red/http/shields.io.svg)](http://shields.io/)


# IntroVerse

## Publication

**IntroVerse: a comprehensive database of introns across human tissues**

Sonia García-Ruiz, Emil K Gustavsson, David Zhang, Regina H Reynolds, Zhongbo Chen, Aine Fairbrother-Browne, Ana Luisa Gil-Martínez, Juan A Botia, Leonardo Collado-Torres, Mina Ryten

**Nucleic Acids Research**, *Volume 51, Issue D1, 6 January 2023, Pages D167–D178,* [https://doi.org/10.1093/nar/gkac1056](https://doi.org/10.1093/nar/gkac1056)

## Website

IntroVerse is freely available for use and data download at: [IntroVerse](https://rytenlab.com/browser/app/introverse)

## Description

IntroVerse is a database of introns and alternative splicing events across human tissues.

IntroVerse is the first version of an effort to provide quality-controlled information on the splicing activity of 332,571 annotated introns (Ensembl v105), spanning 32,669 different genes across 17,510 human control RNA samples from 54 tissues provided by the GTEx v8 project and processed by recount3. 

IntroVerse has two unique features: 

* it provides a complete and unfiltered catalogue of introns that have been carefully processed and quality-controlled
* each novel intron has been assigned to a specific annotated intron, which provides a hierarchical structure allowing the user to start by browsing a gene or genomic position of interest and to end with an exploration of the complete catalogue of alternative events.

In addition, IntroVerse allows users to browse, visualise, download and compare key genomic parameters either across introns of a particular gene and sample of interest or across multiple genes, samples and human tissues.

IntroVerse also provides the visualisation of alternative splicing events by:

* displaying the section of the gene's MANE transcript predicted to be spliced out from the isoform after the excision of the novel event.
* the frequency whereby each intron of a gene's MANE transcript is mis-spliced.



## Software and libraries used to build up the UI:

* [ggtranscript](https://github.com/dzhang32/ggtranscript)
* [recount3](https://rna.recount.bio/)
* [SQLite](https://www.sqlite.org/index.html)
* [Shiny](https://shiny.rstudio.com/)
* [Docker](https://shiny.rstudio.com/)
* [Favicon](https://icons8.com/icons/set/database)
