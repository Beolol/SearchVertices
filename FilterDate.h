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
    QVector< qreal >                        m_colCurvative;
    QVector< qreal >                        m_colCurveRadius;
    QVector< qreal >                        m_colValue;
    QVector< qreal >                        m_colValueDev;
    QVector< qreal >                        m_colValueDev2;
    QVector< qreal >                        m_colValueDerivative;
    QVector< qreal >                        m_colValueSecondDerivative;
    QVector< QPair< qreal, qreal > >        m_colAproxPoly;
    QVector< QPair< qreal, qreal > >        m_colVertices;
    QVector< QPair< qreal, qreal > >        m_colMinVertices;
};


#endif // FILTERDATE_H
