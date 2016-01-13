#-------------------------------------------------
#
# Project created by QtCreator 2014-12-17T14:39:51
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

TARGET = tst_vacc_estest
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += \
    evostrat.cpp \
    main.cpp \
    VaccSEIR_EvoStrat.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    evostrat.h \
    VaccSEIR_EvoStrat.h
