# Individual Based Model (IBM)

This repository is for a generalised IBM that I am trying to generate during my PhD. I abuse the term Individual a bit I actually mean agent or super individual, where one entity can represent
many entities that have identical attributes. The purpose of this is for scaling up results to population levels. It will be a stochastic simulator used 
for investigating different stock assessment methods under hypothetical scenarios. The current repository is a mixture of  ideas that need to be acknowledged. 
Firstly to the Casal2 team found [here](https://github.com/NIWAFisheriesModelling/CASAL2),
which formulated the basis of the core code for error handling, parameter structure, configuration syntax and more. I also learnt about IBM's 
via the model found [here](https://github.com/trophia/sna1), where some idea's have been brought across. And finally to the SPM  team found [here](https://github.com/NIWAFisheriesModelling/SPM),
where most of the spatial code structure comes from.

This is an open source project and if anyone is interested in this project please get in touch the more input the better. I am interested in adding features that are being used in recent applications of IBM's such as; [energetic functionality](https://www.sciencedirect.com/science/article/pii/S0304380018301327) 
habitat based movement, via currents and active gradient searches like in [this paper](https://www.sciencedirect.com/science/article/pii/S0079661117302896)
Also make it a full life cycle model so that we can do [end to end model](https://www-sciencedirect-com.ezproxy.auckland.ac.nz/science/article/pii/S0079661115000233)
It would also need to simulate data that could easily be assess in a stock assessment model for management strategy evaluation.

IBM's are notorious for being limited by CPU, so a big emphasis is to make it as modular and thread safe as possible. It is currently
coded to be used on a desktop (as apposed to HPC) as I belief that is where most users will apply it, although because we use the OpenMP C++ 
library in theory you could make it highly spatial and it should be work fine on a HPC or cluster.

This IBM has a forced spatial structure, where users must define at least one spatial area. The reason I have gone down this road
is mainly because I am interested in spatial characteristics and I believe they are fundamental to IBM's. This means for simple 
spatial models (single area models) the model could be a bit lousy, But when you start having high spatial resolution and spatial
processes such as fishing. The spatial memory management will benefit massively in efficiency's it also allows easier threading.


Currently there is no executable to download and use because it is still in development, but if you are interested in the code base
and want to be a collaborator on the project I recommend you fork the repository and compile the code which is easy to do as shown in the
following steps. There will be 

### Compiling the code

Open a terminal in the BuildSystem folder, to check you have the right dependencies on your program run the following command

```
doBuild.bat check
```

if you have windows 10 you may have to put a ```.\``` in front of the command, and for linux

```
.\doBuild.sh check
```

This will look for g++, git, and g++ version maybe some others.

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

### Running the Program

There is a partial user manual found in the Documentation folder, but a quick fire run through. This IBM program runs through the command line
(no gui), using the executable name and run parameters e.g. ```ibm.exe -r --loglevel finest > output.log 2> err.log``` This command will run (```-r```) the ibm, which by
default looks for a file named ```config.ibm``` to look for specific model tasks. The main output defined in the configuration files will be printed in
 the file labeled ```output.log```. I have also added an additional print statement (```--loglevel```) which will print extra 
information at run time into the file ```err.log```. There are 4 levels of ```loglevel```, ```medium```,  ```fine```, ```finest``` and ```trace```. If anyone is interested you can 
look in the code base for statements such as ```LOG_FINEST()``` which will print information if the parameters is defined at run time.

God speed