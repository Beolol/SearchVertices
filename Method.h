#ifndef METHOD_H
#define METHOD_H
#include <cmath>
#include <iomanip>
#include <math.h>
#include <QVector>
#include <QPair>

class Polinom
{
public:
    Polinom( size_t NSize, size_t Ndegree, const QVector< QPair< qreal, qreal > > & function );
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

    QVector< QPair< qreal, qreal > >        minVertices;
    QVector< QPair< qreal, qreal > >        maxVertices;
    QVector< QVector< qreal > >             massiv;
    QVector< qreal >                        inputX;
    QVector< qreal >                        inputY;
    QVector< qreal >                        vectorC;
    QVector< qreal >                        derivative;
    size_t                                  size;
    size_t                                  degree;

};

#endif // METHOD_H
