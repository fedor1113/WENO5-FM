TEMPLATE = app
CONFIG += console c++2a
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -ltbb

INCLUDEPATH += Eigen/

SOURCES += \
        _vector4.cpp \
		main.cpp \
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
		kfr1d.h \
		ssprk33.h \
		weno5.h \
		weno5coefs.h
