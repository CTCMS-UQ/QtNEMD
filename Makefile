all:
	# F2Py module
	f2py3 driver.f Walls_SSGK.f -m TTCF -h TTCF.pyf --overwrite-signature 
	f2py3 -c TTCF.pyf Walls_SSGK.f driver.f --fcompiler=gfortran \
	--f77flags='-Ofast -march=native -cpp' \
	--f90flags='-Ofast -march-native -cpp' -lgomp
	# GUI elements
	pyuic5 GUI-resources/mainwindow.ui > GUI-resources/ui_mainwindow.py
	pyuic5 GUI-resources/floating_plot.ui > GUI-resources/ui_floating_plot.py
	# Bare exe
	gfortran -march=native -g -O3 -cpp -fopenmp -o Walls_SSGK Walls_SSGK.f

driver: 
	# F2Py module
	f2py3 driver.f Walls_SSGK.f -m TTCF -h TTCF.pyf --overwrite-signature 
	f2py3 -c TTCF.pyf Walls_SSGK.f driver.f --fcompiler=gfortran \
	--f77flags='-Ofast -march=native -cpp -fopenmp' \
	--f90flags='-Ofast -march-native -cpp -fopenmp' -lgomp
gui:
	pyuic5 GUI-resources/mainwindow.ui > GUI-resources/ui_mainwindow.py
	pyuic5 GUI-resources/floating_plot.ui > GUI-resources/ui_floating_plot.py

bare:
	gfortran -march=native -g -ffpe-trap=invalid,zero,overflow -O0 -cpp -fopenmp -o Walls_SSGK Walls_SSGK.f

.PHONY: clean
.PHONY: veryclean
clean:
	rm TTCF.pyf TTCF.cpython* Walls_SSGK
veryclean:	
	rm TTCF.pyf TTCF.cpython* 
	rm GUI-resources/ui_mainwindow.py GUI-resources/ui_floating_plot.py
