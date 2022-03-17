all:
	# F2Py module
	f2py3  shear.f -m TTCF -h TTCF.pyf --overwrite-signature 
	f2py3 -c TTCF.pyf shear.f  --fcompiler=gfortran \
	--f77flags='-Ofast -march=native -cpp' \
	--f90flags='-Ofast -march-native -cpp'
	# GUI elements
	pyuic5 GUI-resources/mainwindow.ui > GUI-resources/ui_mainwindow.py
	pyuic5 GUI-resources/floating_plot.ui > GUI-resources/ui_floating_plot.py

driver: 
	# F2Py module
	f2py3  shear.f -m TTCF -h TTCF.pyf --overwrite-signature 
	f2py3 -c TTCF.pyf shear.f  --fcompiler=gfortran \
	--f77flags='-Ofast -march=native -cpp' \
	--f90flags='-Ofast -march-native -cpp'
gui:
	pyuic5 GUI-resources/mainwindow.ui > GUI-resources/ui_mainwindow.py
	pyuic5 GUI-resources/floating_plot.ui > GUI-resources/ui_floating_plot.py

.PHONY: clean
.PHONY: veryclean
clean:
	rm TTCF.pyf TTCF.cpython* 
veryclean:	
	rm TTCF.pyf TTCF.cpython* 
	rm GUI-resources/ui_mainwindow.py GUI-resources/ui_floating_plot.py
