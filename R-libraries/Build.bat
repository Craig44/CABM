R --vanilla < run-roxygen.R
R --vanilla < IBM_make_version.R

R CMD build --force ibm
R CMD INSTALL --build ibm

del IBM.html




