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
    ebdf5.h \
		eno3.h \
    eos.h \
		euler1d.h \
    eulerforward.h \
    exactsolver.h \
		extrapolation.h \
    hllc_solver.h \
    inviscidburgers1d.h \
		inviscidburgers1d.h \
    kfr1d.h \
		kfr1d.h \
    lf_flux.h \
    lssperk10_9.h \
    lssperk12_11.h \
    miegruneisen.h \
    ossperk4_3.h \
    rk4_4.h \
    rk6_5.h \
    roe.h \
    ssprk10_4.h \
		ssprk33.h \
    ssptserk12_8.h \
    ssptserk8_5.h \
    tdrk3_5.h \
    twoderrk3_5.h \
		weno5.h \
		weno5coefs.h

DISTFILES += \
    .gitignore \
	plot.gnuplot \
    plot.py
