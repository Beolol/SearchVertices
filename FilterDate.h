#ifndef FILTERDATE_H
#define FILTERDATE_H
#include <cmath>
#include <iomanip>
#include <math.h>
#include <QVector>
#include <QPair>
#include <QScopedPointer>
#include "Method.h"

class CFilterDate
{
public:
    CFilterDate(QVector<qreal> &date, const size_t Ndegree
                                            , const size_t SizeRunningWindow = 16);
    QVector< qreal >                        GetAproxPoly();
    void                                    SavePolyToFile( const QString & filename );
    QVector< QPair< qreal, qreal > >        GetVertices( qreal dInputA = 0, qreal dInputB = 127 );
    QVector< QPair< qreal, qreal > >        GetMinVertices( qreal dInputA = 0, qreal dInputB = 127 );
    qreal                                   GetDerivative( const qreal &x );
    qreal                                   GetSecondDrivative( const qreal &x );
    qreal                                   GetCurveRadius( const qreal &x );
    qreal                                   GetCurvative( const qreal &x );
    static QPair< QVector< qreal >, QVector< qreal > >
                                            GetDebPoly
    (
        QVector< qreal > &                  sData,
        const size_t                        ndeg,
        const size_t                        from
    );
private:
    void                                    XTrAlgorithm( const size_t & lBorder, const size_t & rBorder );
    void                                    PolySmoorthM(QVector<qreal> &sData, const size_t ndeg, const size_t index);
    void                                    ConstructionPolynomial
                                            (
                                                QVector<qreal> &date, const size_t Ndegree,
                                                const size_t from, const size_t to
                                             );
    void                                    Filter (QVector<qreal> &date, const size_t Ndegree,const size_t SizeRunningWindow
                                                    );
    QVector< qreal >                        curvative;
    QVector< qreal >                        curveRadius;
    QVector< qreal >                        value;
    QVector< qreal >                        valueDev;
    QVector< qreal >                        valueDev2;
    QVector< qreal >                        valueDerivative;
    QVector< qreal >                        valueSecondDerivative;
    QVector< QPair< qreal, qreal > >        aproxPoly;
    QVector< QPair< qreal, qreal > >        vertices;
    QVector< QPair< qreal, qreal > >        minVertices;
};


#endif // FILTERDATE_H
