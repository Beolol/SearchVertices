TEMPLATE = app
QT += core
CONFIG += console
CONFIG -= app_bundle
CONFIG += c++11

#QMAKE_CXXFLAGS += -std=c++0x #mingw/gcc only!
SOURCES += main.cpp \
    Method.cpp \
    FilterDate.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    Method.h \
    FilterDate.h

