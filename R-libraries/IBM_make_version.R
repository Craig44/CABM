# Get IBM Version
VERSION = ""
this_dir = getwd();
if (grepl("unix", .Platform$OS.type)) {
  setwd("../BuildSystem/bin/linux/release")
  VERSION<-system("./ibm -v",intern=TRUE)
  setwd(this_dir)
} else {
  VERSION<-system("ibm.exe -v",intern=TRUE)
  setwd("../BuildSystem/bin/windows/release")
  VERSION<-system("ibm -v",intern=TRUE)
  setwd(this_dir)
}


version.date<-substring(VERSION,1,regexpr(" ",VERSION)-1)
version.ibm_number<- substring(VERSION,regexpr("rev. ",VERSION) + 5,nchar(VERSION) - 1)
version.package<- "0.1"


# Build DESCRIPTION file
filename<-"ibm/DESCRIPTION"
cat("Package: ibm\nTitle: ibm extract package\nVersion: ",file=filename)
cat(version.package,file=filename,append=TRUE)
cat("\nDate: ",file=filename,append=TRUE)
cat(version.date,file=filename,append=TRUE)
cat("\n",file=filename,append=TRUE)
cat("Author: C. Marsh\n",file=filename,append=TRUE)
cat("Description: A set of R functions for extracting and plotting from ibm output files.\n",file=filename,append=TRUE)
cat("Maintainer: Craig Marsh <Craig.Marsh@bigbluedata.co.nz>\n",file=filename,append=TRUE)
cat("License: GPL v 2. See the IBM User Manual for license details.\n",file=filename,append=TRUE)
cat("LazyData: true\n",file=filename,append=TRUE)
cat("Import:\n",file=filename,append=TRUE)
cat("\t reshape2,\n",file=filename,append=TRUE)
cat("\t ggplot2,\n",file=filename,append=TRUE)
cat("\t RColorBrewer\n",file=filename,append=TRUE)


cat("BugReports: https://github.com/Craig44/IBM/issues\n",file=filename,append=TRUE)


# Create R function to return version number
filename<-"ibm/R/IBM.binary.version.R"
cat("\"IBM.binary.version\"<-\n",file=filename)
cat("function() {\n",file=filename,append=T)
cat(paste("return(\"",VERSION,"\")\n",sep=""),file=filename,append=T)
cat("}\n",file=filename,append=T)

# Create R function to return R library version number
filename<-"ibm/R/Version.R"
cat("\"Version\"<-\n",file=filename)
cat("function() {\n",file=filename,append=T)
cat(paste("return(\"",version.date, "(", version.ibm_number, ")\")\n",sep=""),file=filename,append=T)
cat("}\n",file=filename,append=T)

# Write a .html file to report version number for the Wiki
cat(paste("Version",VERSION),file="ibm.html")

# Exit
q()

