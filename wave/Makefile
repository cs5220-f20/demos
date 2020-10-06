waves.so: _waves.pyx waves.c
	python setup.py build_ext -i

wave_demo.x: wave_demo.c waves.c
	gcc wave_demo.c waves.c -o wave_demo.x

clean:
	rm -f _waves.c *.so
