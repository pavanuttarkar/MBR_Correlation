#rm *.o
#rm *.so
python setup.py build_ext --inplace
mkdir PROC-CPP
mv *.cpp PROC-CPP/.
