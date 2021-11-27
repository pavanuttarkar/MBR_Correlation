.PHONY: build

build:
	python2.7 setup.py build_ext --inplace

.PHONY: docs

docs:
	make html --directory docs

.PHONY: clean

clean:
	rm -rf build/*
	rm -rf api/*.so api/*.cpp
