#rm *.o
#rm *.so
#python setup_header_Fring.py build_ext --inplace
#python setup_header_read.py build_ext --inplace
#python setup_header_geometric.py build_ext --inplace
rm call_to_read_dev.cpp
rm header_Fringe_Search.c
rm geometric.c
python2 setup.py build_ext --inplace
