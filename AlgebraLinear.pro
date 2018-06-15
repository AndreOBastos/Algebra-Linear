TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    LUDec.cpp \
    method.cpp \
    qrdec.cpp \
    eigen.cpp \
    svdec.cpp \
    pca.cpp

HEADERS += \
    LUDec.h \
    Mat.h \
    method.h \
    qrdec.h \
    eigen.h \
    svdec.h \
    pca.h
