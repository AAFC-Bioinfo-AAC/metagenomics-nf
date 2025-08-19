
# Building a Custom Apptainer Container for Use with Nextflow

These steps need to be done on a computer with root access.  It is a good practice to build apptainer images on a local computer and transfer them to your HPC. Find yml files to build all apptainer images for this pipeline [here](../docs/apptainer_images.md)

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



## Step 3 - Convert the .sif file to sandbox in order to modify the image

As is, the image won't work proprely in Nextflow pipelines. Minor changes need to be done, but as a sif file is read-only, it needs to be converted to a sandbox format to modify it.

    sudo apptainer build --sandbox wheezy.dir kraken2.sif

Then the sif image is moved away because we will create another image with the same name later.

    mkdir tmp; mv kraken2.sif tmp

## Step 4 - Modify the image by entering the shell

Note that we always use -c option.

-c, --contain use minimal /dev and empty other directories (e.g. /tmp and $HOME) instead of sharing filesystems from your host

	    sudo apptainer shell --writable -c wheezy.dir/
        whoami

You are logged in as root.

### Change the permissions of the /tmp

The first thing to do is simply to change the permissions of the /tmp folder. Otherwise, apt-get install won't work. 

    chmod 1777 /tmp
    
### Add the UNIX ps executable that is required by Nextflow engine

    apt update
    apt-get install procps

### Make the conda executable part of a folder that is in the PATH

At this stage the image will work properly if using Apptainer run [...]. You can try it by quitting the shell of the container :

    exit
    apptainer run -c wheezy.dir/ kraken2 -v
  
**Output:**

    Kraken version 2.1.3
    Copyright 2013-2023, Derrick Wood (dwood@cs.jhu.edu)

    
Now if we try the same thing with apptainer exec, which is the commmand used by Nextflow, it is not working! :

    apptainer exec -c wheezy.dir/ kraken2 -v

**output:**

    FATAL:   "kraken2": executable file not found in $PATH

Is is because "apptainer run" allows to run runscripts that are present inside the containers. Here these scripts clearly activate the micromamba base environment contrary to the plain "apptainer exec" command.

A workaround is simply to put all the conda execuables, which are in /opt/conda/bin/, in a folder that is in the PATH variable. Rather than copying these executables, we will make symbolic links in /usr/local/bin.
	
    # Entering the container with the possibility to made changes
    sudo apptainer shell --writable -c wheezy.dir/

	# making sure to be in /usr/local/bin before creating links
	cd /usr/local/bin; for i in `ls /opt/conda/bin/`; do ln -s /opt/conda/bin/$i; done
    kraken2 -v
	exit

## Step 5 - Convert the modified sandbox image to a read-only .sif image

	 sudo apptainer build kraken2.sif wheezy.dir

	 
	 
## Step 6 - Test the image

    apptainer exec -c kraken2.sif kraken2 -v; apptainer exec -c kraken2.sif bracken -v
    

## Step 7 - Clean up

    sudo rm -rf wheezy.dir
    sudo rm -rf tmp

Now you are done!

## Step 8 - Modify an existing image 

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

