# Individual Based Model (IBM)

This repository is for a generalised IBM that I am trying to generate during my PhD. I abuse the Individual term a bit I actually mean agent or super individual, where one entity 'can' represent
many entities that have identical attributes. The purpose of this is for scaling up results to population levels. It will be a stochastic simulator used 
for investigating different stock assessment methods under hypothetical scenerios. The current repository is a mixture of  ideas that need to be acknowledged. 
Firstly to the Casal2 team found [here](https://github.com/NIWAFisheriesModelling/CASAL2),
which formulated the basis of the core code for error handling, parameter structure, configuration syntax and more. I also learnt about IBM's 
via the model found [here](https://github.com/trophia/sna1) where some idea's have been brought across.

This is an open source project and if anyone is interested in this project please get in touch the more input the better. Although there is a lot
code that has been pushed to the repository it is still in the design phase currently, have a look at this document for my current thoughts on code 
structure model design see this [document](https://github.com/Craig44/IBM/blob/master/initial%20thoughts.docx).

I am interested in adding features that are being used in recent applications of IBM's such as; [energetic functionality](https://www.sciencedirect.com/science/article/pii/S0304380018301327) 
habitat based movement, via currents and active gradient searches like in [this paper](https://www.sciencedirect.com/science/article/pii/S0079661117302896)
Also make it a full life cycle model so that we can do [end to end model](https://www-sciencedirect-com.ezproxy.auckland.ac.nz/science/article/pii/S0079661115000233)
It would also need to simulate data that could easily be assess in a stock assessment model for management strategy evaluation.

IBM's are notorious for being limited by CPU, so a big emphasis is to make it as modular and thread safe as possible. It is currently
coded to be used on a desktop (as apposed to HPC) as I belief that is where most users will apply it.

This IBM has a forced spatial structure, where users must define at least one spatial area. The reason I have gone down this road
is mainly because I am interested in spatial characteristics and I believe they are fundamental to IBM's. This means for simple 
spatial models (single area models) the model could be a bit lousy, But when you start having high spatial resolution and spatial
processes such as fishing. The spatial memory management will benefit massively in efficiency's it also allows easier threading.


Currently there is no executable to download and use because it is still in development, but if you are interested in the code base
and want to be a collaborator on the project I recommend you fork the repository and compile the code which is easy to do as shown in the
following steps.

Open a terminal in the BuildSystem folder, to check you have the right dependencies on your program run the following command

```
doBuild.bat check
```

if you have windows 10 you may have to put a ```.\``` in front of the command, and for linux

```
doBuild.sh check
```

This will look for g++, git, gfortran (I don't think this is neccesary though will check)

If the check worked okay, run the command to build the within repository libraries Python and Boost

```
doBuild.bat thirdparty
```

Once that is done you should (in theory) be able to compile the code and build an executable using the following
call in the command.

```
doBuild.bat release
```

This will put the executable in the following directory BuildSystem\bin\'system'\release, where 'system' can 
be windows or linux depending on your OS.

God speed