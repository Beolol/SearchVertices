#ifndef METHOD_H
#define METHOD_H
#include <cmath>
#include <iomanip>
#include <math.h>
#include <QVector>
#include <QPair>

class CPolinom
{
public:
    CPolinom( size_t NSize, size_t Ndegree, const QVector< QPair< qreal, qreal > > & function );
    void                                    LeastSquareMethod();
    void                                    PrintMassive();
    void                                    Gauss();
    void                                    CalcDerivative();
    void                                    SearchRoots( const qreal & lBorder, const qreal & section );
    void                                    Frebenius();
    qreal                                   ValueFunction(const qreal &x );
    QVector< QPair< qreal, qreal > >        ReturnFunction();
    QVector< QPair< qreal, qreal > >        ReturnMaxVertices();
    QVector< QPair< qreal, qreal > >        ReturnMinVertices();
    qreal                                   ValueDerFunction( const qreal & x );
    qreal                                   ValueSecondDerFunction( const qreal & x );
private:
    qreal                                   MethodSecants( qreal & x1, qreal & x2, const qreal & e );
    qreal                                   Newton( qreal & x1, qreal & x2, const qreal & e );

private:

    QVector< QPair< qreal, qreal > >        m_colMinVertices;
    QVector< QPair< qreal, qreal > >        m_colMaxVertices;
    QVector< QVector< qreal > >             m_colMassiv;
    QVector< qreal >                        m_colInputX;
    QVector< qreal >                        m_colInputY;
    QVector< qreal >                        m_colVectorC;
    QVector< qreal >                        m_colDerivative;
    size_t                                  m_iSize;
    size_t                                  m_iDegree;

};

#endif // METHOD_H
