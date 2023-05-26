#rm *.o
#rm *.so
python3 setup.py build_ext --inplace
mkdir PROC-CPP
mv *.cpp PROC-CPP/.
