
# Building a Custom Apptainer Container for Use with Nextflow

These steps need to be done on a computer with root access.  It is a good practice to build apptainer images on a local computer and transfer them to your HPC. Find yml files to build all apptainer images for this pipeline [here](https://github.com/AAFC-Bioinfo-AAC/metagenomics-nf/tree/main/docs/apptainer_images)

We will use [Apptainer Definition files](https://apptainer.org/docs/user/latest/definition_files.html#definition-files) provided by the [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) team and in the example provided here, we will build a kraken2 image that can be used with the metagenomic_nf workflow­. 

## Step 1 - Modify the environment.yml file

In the dependencies section, add programs and their version needed in the image. It is recommended to add a python version as well.

The name of the default environment (base) cannot be changed, but any relevant channel, e.g. the bioconda channel, can be added.

For a kraken2 image, the yml file should look like this :

```shell
name: base
channels:
  - conda-forge
  - defaults
  - bioconda
dependencies:
  - python=3.10
  - kraken2=2.1.2
  - bracken=2.9
```

## Step 2 - Build the apptainer image using the definition file

With this approach this is as simple as :

    apptainer build kraken2.sif micromamba.def
 
	 
## Step 3 - Test the image

    apptainer exec -c kraken2.sif kraken2 -v; apptainer exec -c kraken2.sif bracken -v

## Step 4 - Modify an existing image 

You realized that you have forgotten to add some important scripts in your image.

Regenerate the sandox

    sudo apptainer build --sandbox wheezy.dir kraken2.sif


Enter the image using the shell and bind the host ../misc folder that contains your scripts to an empty repository that exists in the image (/tmp)

    sudo apptainer shell -c -B ../misc:/tmp --writable wheezy.dir/

    Apptainer> ls /tmp/
    combine_mpa.py	kreport2mpa.py	produce_bracken_nf.sh
    cp -v /tmp/* /usr/local/bin/
    exit

Make sure that it is working :

    apptainer exec -c wheezy.dir/ combine_mpa.py
    usage: combine_mpa.py [-h] -i IN_FILES [IN_FILES ...] -o O_FILE
    combine_mpa.py: error: the following arguments are required: -i/--input, -o/--output

Rebuild the image

    # keeping the older .sif image in a safe place
    mkdir tmp; mv kraken2.sif tmp
    
    sudo apptainer build kraken2.sif wheezy.dir

Now you are done!

## Note 

To use apptainer images in Nextflow, these lines should be present in the profile :

```shell
laptop {
        
    apptainer.enabled = true

    // scope process
   container = "$baseDir/containers/kraken2.sif"

}
```

You can build multiple container images at once in a loop, e.g.;

```bash
for folder in `find apptainer_images/* -type d`; do
    pushd "$folder"
    singularity build "${folder}.sif" micromamba.def
    popd
done
```