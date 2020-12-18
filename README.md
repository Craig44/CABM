# Agent/Individual Based Model (ABM)

This repository is for a generalised ABM that I am developing during my PhD. This apps purpose is as an operating model for testing and validating ecological 
models. The ABM tracks agents which can represent many individuals that have identical attributes over a spatial explicit domain. The ABM allows users to specify 
different life-histories, agent dynamics and interactions with fisheries. The current repository is a mixture of  ideas that need to be acknowledged. 
Firstly to the Casal2 team found [here](https://github.com/NIWAFisheriesModelling/CASAL2),
which formulated the basis of the core code for error handling, parameter structure, configuration syntax and more. I also learnt about ABM's 
via the model found [here](https://github.com/trophia/sna1), where some idea's have been brought across. And finally to the SPM  team found [here](https://github.com/NIWAFisheriesModelling/SPM),
where most of the spatial code structure comes from.

This is an open source project and if anyone is interested in this project please get in touch the more input the better. Because I will be using this as an operating model during my PhD, I would really appreciate
if people who use this tool get in touch to minimize over lap and also to just share ideas. I am interested in adding features that are being used in recent applications of ABM's such as; [energetic functionality](https://www.sciencedirect.com/science/article/pii/S0304380018301327) 
habitat based movement, via currents and active gradient searches like in [this paper](https://www.sciencedirect.com/science/article/pii/S0079661117302896)
Also make it a full life cycle model so that we can do [end to end model](https://www-sciencedirect-com.ezproxy.auckland.ac.nz/science/article/pii/S0079661115000233)
It would also need to simulate data that could easily be used in a stock assessment model for management strategy evaluation.

IBM's are notorious for being limited by CPU, so a big emphasis is to make it as modular and thread safe as possible. It is currently
coded to be used on a desktop (as apposed to HPC) as I belief that is where most users will apply it, although because we use the OpenMP C++ 
library in theory you could make it highly spatial and it should be work fine on a HPC or cluster.

This IBM has a forced spatial structure, where users must define at least one spatial area. The reason I have gone down this road
is mainly because I am interested in spatial characteristics and I believe they are fundamental to IBM's. This means for simple 
spatial models (single area models) the model could be a bit lousy, But when you start having high spatial resolution and spatial
processes such as fishing. The spatial memory management will benefit massively in efficiency's it will also allow easier threading.


### Current Status of Repo
Currently users can have any spatial resolution they desire. I have an example model in found in this [directory](https://github.com/Craig44/IBM/tree/master/Example/SpatialModel) that is 5x5 spatial model,
with a single fishery. This model can generate a range of age, length compositional data and biomass indicies. The parallelization is not working as well as I had hoped
so if you run the model use max_threads 1.
Currently there is no executable to download and use because it is still in development, but if you are interested in the code base
and want to be a collaborator on the project I recommend you fork the repository and compile the code which is easy to do as shown in the
following steps.

### Compiling the code

Since the MSE mode has been added the ABM relies on R being in the system path, to check this open a terminal and type in R. As well as R being in the path, compilation requires the following packages to be installed.
```
install.packages(c("RInside","Rcpp"));
```
Once those dependencies have been addressed open a terminal in the BuildSystem directory, to check you have the right dependencies on your program run the following command

```
doBuild.bat check
```

if you have windows 10 and you are building from the power shell you may have to put a ```.\``` in front of the command, and for linux

```
.\doBuild.sh check
```

This will look for g++, git, R, and g++ version maybe some others.

If the check worked okay, run the command to build the within repository libraries Python, cmake and Boost

```
doBuild.bat thirdparty
```

Once that is done you should (in theory) be able to compile the code and build an executable using the following
call in the command.

```
doBuild.bat release
```

This will put the executable in the following directory ```BuildSystem\bin\'system'\release```, where ```'system'``` can 
be ```windows``` or ```linux``` depending on your OS. An alternative is to build it in a IDE such as eclipse.

For more information on build parameters look at the python scripts or run the command

```
doBuild.bat help
```

### Tips for compiling on windows
I use the cross gcc compiler found [here](https://jmeubank.github.io/tdm-gcc/). Make sure when you install, that openMP and I think fortran tool sets are available in the gcc compiler. 
If you use the compiler that is linked above, this will involve manually selecting extra settings in the install process.


Avoid white space in the absolute path name of the repository, in the past we have found issues when white space is present in the absolute path of the repository.


### Running the Program

There is a partial user manual found in the Documentation folder, but a quick fire run through. This IBM program runs through the command line
(no gui), using the executable name and run parameters e.g. ```ibm.exe -r --loglevel finest > output.log 2> err.log``` This command will run (```-r```) the ibm, which by
default looks for a file named ```config.ibm``` to look for specific model tasks. The main output defined in the configuration files will be printed in
 the file labeled ```output.log```. I have also added an additional print statement (```--loglevel```) which will print extra 
information at run time into the file ```err.log```. There are 4 levels of ```loglevel```, ```medium```,  ```fine```, ```finest``` and ```trace```. If anyone is interested you can 
look in the code base for statements such as ```LOG_FINEST()``` which will print information if the parameters is defined at run time.

Godspeed