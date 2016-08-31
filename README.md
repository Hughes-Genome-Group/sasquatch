# Sasquatch

### Description

Sasquatch uses DNase-seq data to pile-up average DNase I cut profiles over all possible short k-mers in the context of open-chromatin. By identifying and quantifying footprints in these average profiles, Sasquatch infers a k-mers potential to bind transcription factors in a tissue-specific manner. Based on comparative analysis, Sasquatch predicts the damaging potential of sequence variations considering the tissue context of interest and independent of a specific genotype. Furthermore, *in silico* mutation analysis profiles larger sequences for transcription factor binding sites actively bound in the tissue of interest.
With a large repository of preprocessed, publicly available DNase-seq data, complemented by our own, deep DNase-seq in human primary erythroid cells, Sasquatch provides a powerful tool for rapid transcription factor profiling and for prioritising and interpreting non-coding sequence variations.

### Getting Started 

1. **Clone** the repository.

2. **Get data**. The repository only comes with a minimal dummy of example data to run the example scripts. To start with real data, visit our [webtool](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi "webtool") site and download your data of interest from our [repository](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi?rm=mode_71 "Repository") of preprocessed DNase-seq data. Details about the different samples are avaialble [here](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi?rm=mode_75 "Details"). Extract them into *./sasquatch/data/human/DNase* or *./sasquatch/data/mouse/DNase*. Every tissue data comes in a directory which should form the subdirectory in your *sasquatch/data/organism/DNase* local repository. (When extracting into any other directory just make sure that you link Sasquatch to your right personal repository and keep every tissue dataset in a separate subdirectory.)

3. In your **R-script**, source the the R-functions *./R_utility/functions_sasq_r_utility.R*, point to your local data repository and set basic parameters as described in the beginning of the [Vignette](./R_utility/sasq_R_utility_vignette.pdf "Vignette") and [Example R-script](./R_utility/example_script_sasq_r_utility.R "Example").

4. Also make sure to check out our [webtool](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi "webtool") for further documentation about Sasquatch and average DNase-footprints analysis.


### Documentation

[Vignette](./R_utility/sasq_R_utility_vignette.pdf "Vignette") running through the basic Sasquatch analysis steps.

[Reference Manual](./R_utility/sasq_R_reference_manual.pdf "Manual") for all implemented functions.

[Example R-script](./R_utility/example_script_sasq_r_utility.R "Example") running through the most of Sasquatchs functions.

Also check out our general introduction into the [Sasquatch](http://apps.molbiol.ox.ac.uk/sasquatch/docs/Manual.pdf) approach and into [average DNase I footprints](http://apps.molbiol.ox.ac.uk/sasquatch/docs/Introduction_to_Average_Footprints.pdf "footprints") from our [webtool](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi "webtool") site. 

### Links

Here, you can reach our [Sasquatch webtool](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi "webtool") implementation as well as the [repository](http://apps.molbiol.ox.ac.uk/sasquatch/cgi-bin/foot.cgi?rm=mode_71 "Repository") of preprocessed DNase-seq data to download for local use. 

### License
Sasquatch is pulished under [GPLv3](./LICENSE.txt "License") or later.
