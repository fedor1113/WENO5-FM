TEMPLATE = app
CONFIG += console c++2a
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -ltbb
LIBS += -ltbb

INCLUDEPATH += Eigen/

SOURCES += \
        _vector4.cpp \
		main.cpp \
    miegruneisen.cpp \
    valarray_extra_operations.cpp

HEADERS += \
        _vector4.h \
		arithmeticwith.h \
		eno3.h \
		euler1d.h \
    eulerforward.h \
    exactsolver.h \
		extrapolation.h \
    inviscidburgers1d.h \
		inviscidburgers1d.h \
    kfr1d.h \
		kfr1d.h \
    lf_flux.h \
    miegruneisen.h \
    rk4_4.h \
    rk6_5.h \
    ssprk10_4.h \
		ssprk33.h \
    tdrk3_5.h \
    twoderrk3_5.h \
		weno5.h \
		weno5coefs.h

DISTFILES += \
    .gitignore \
	plot.gnuplot \
    plot.py
