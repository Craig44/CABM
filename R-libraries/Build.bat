R --vanilla < CABM.binary.version
R --vanilla < run-roxygen.R

R CMD build --force ibm
R CMD INSTALL --build ibm

del IBM.html
rm ibm.html



