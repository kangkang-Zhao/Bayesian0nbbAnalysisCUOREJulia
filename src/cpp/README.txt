Instructions for Julia calls cpp code
#1 compile cpp program to executable file
g++ readHist.cxx -o readHist `root-config --cflags --libs`
# run julia code
