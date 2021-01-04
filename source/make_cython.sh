#rm *.o
#rm *.so
python setup_header_Fring.py build_ext --inplace
python setup_header_read.py build_ext --inplace
python setup_header_geometric.py build_ext --inplace
